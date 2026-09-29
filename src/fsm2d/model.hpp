// The 2D FSM model (Charlesworth & Turner, PNAS 2019, continuous visual state: SI section 2), shared by the CPU and
// CUDA simulators: parameters, initial state, the decision rule and the motion update. docs/NOTES.md explains the
// model and the score.
//
// Every step each agent considers the five moves (turn by -dtheta or +dtheta at speed v0; keep its heading at speed
// v0 - dv, v0 or v0 + dv), imagines every sequence of tau moves starting with each, assumes everyone else flies
// straight on at v0, discards futures that collide, and picks the first move whose surviving futures have the most
// varied visual states. All agents decide from the same snapshot, then move together.
#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../common/cli.hpp"
#include "../common/json.hpp"
#include "../common/rng.hpp"

namespace fsm2d {

constexpr int NM = 5;  // moves: 0 turn -dtheta, 1 turn +dtheta, 2/3/4 straight on at v0 - dv / v0 / v0 + dv
constexpr double TWO_PI = 2.0 * M_PI;
struct Move { int dh, ds; };  // heading change (in dtheta) and speed change (in dv)
constexpr std::array<Move, NM> MOVES{{{-1, 0}, {1, 0}, {0, -1}, {0, 0}, {0, 1}}};

struct Params {
    int N = 50, tau = 4, steps = 1000;
    double v0 = 10, dv = 2, dtheta_deg = 15, R = 1, dt = 1;
    bool half_asin = false;   // half-angle subtended by a disc: atan(R/d) (as in the paper) or asin(R/d) (exact)
    double noise = 0;         // sd of Gaussian noise added to heading and speed after each move
    double box = 100;         // initial positions uniform in [0, box]^2 ...
    double heading_sd_deg = -1;  // ... headings ~ N(0, sd) (default: dtheta)
    uint64_t seed = 1;
    std::string out = "run", init;  // init: text file with "x y theta [v]" per agent (replaces N, box, heading-sd)
    int log_every = 1;
    bool resume = false;

    double dth() const { return dtheta_deg * M_PI / 180.0; }

    // Handles the flags common to both simulators; false if `flag` is not one of them.
    bool parse(const std::string& flag, fsm::Args& a) {
        if (flag == "--N") N = a.integer();
        else if (flag == "--tau") tau = a.integer();
        else if (flag == "--steps") steps = a.integer();
        else if (flag == "--v0") v0 = a.real();
        else if (flag == "--dv") dv = a.real();
        else if (flag == "--dtheta") dtheta_deg = a.real();
        else if (flag == "--R") R = a.real();
        else if (flag == "--dt") dt = a.real();
        else if (flag == "--half-angle") half_asin = a.choice({"atan", "asin"}) == "asin";
        else if (flag == "--noise") noise = a.real();
        else if (flag == "--box") box = a.real();
        else if (flag == "--heading-sd") heading_sd_deg = a.real();
        else if (flag == "--seed") seed = a.u64();
        else if (flag == "--out") out = a.str();
        else if (flag == "--init") init = a.str();
        else if (flag == "--log-every") log_every = a.integer();
        else if (flag == "--resume") resume = true;
        else return false;
        return true;
    }
    void finish() {
        if (heading_sd_deg < 0) heading_sd_deg = dtheta_deg;
        if (N < 2 || tau < 1 || steps < 0) throw std::runtime_error("need --N >= 2, --tau >= 1, --steps >= 0");
    }

    // The model's part of PREFIX.json.
    void describe(fsm::Json& j, int n_agents) const {
        j.num("dims", 2).num("N", n_agents).num("tau", tau)
            .num("v0", v0).num("dv", dv).num("dtheta_deg", dtheta_deg).num("R", R).num("dt", dt)
            .str("half_angle", half_asin ? "asin" : "atan").num("noise", noise)
            .num("box", box).num("heading_sd_deg", heading_sd_deg).num("seed", (unsigned long long)seed).str("init", init);
    }
};

// The simulation state: position, heading and speed of every agent.
struct State {
    std::vector<double> x, y, th, v;
    int size() const { return (int)x.size(); }

    // .traj frame layout: [N][4] doubles (x, y, theta, v)
    std::vector<char> frame() const {
        std::vector<char> buf(sizeof(double) * 4 * x.size());
        double* d = (double*)buf.data();
        for (size_t i = 0; i < x.size(); i++) { d[4 * i] = x[i]; d[4 * i + 1] = y[i]; d[4 * i + 2] = th[i]; d[4 * i + 3] = v[i]; }
        return buf;
    }
    void load_frame(const std::vector<char>& buf) {
        const double* d = (const double*)buf.data();
        for (size_t i = 0; i < x.size(); i++) { x[i] = d[4 * i]; y[i] = d[4 * i + 1]; th[i] = d[4 * i + 2]; v[i] = d[4 * i + 3]; }
    }
};

// Random positions in [0, box]^2 and headings ~ N(0, heading_sd), or the --init file.
inline State initial_state(Params& p) {
    State s;
    if (!p.init.empty()) {
        std::ifstream f(p.init);
        if (!f) throw std::runtime_error("cannot open " + p.init);
        std::string line;
        while (std::getline(f, line)) {
            std::istringstream ls(line);
            double x, y, th, v = p.v0;
            if (!(ls >> x >> y >> th)) continue;
            ls >> v;
            s.x.push_back(x); s.y.push_back(y); s.th.push_back(th); s.v.push_back(v);
        }
        p.N = s.size();
        return s;
    }
    auto rng = fsm::make_rng(p.seed, fsm::kInitStream, 0);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    std::normal_distribution<double> G(0.0, 1.0);
    const double sd = p.heading_sd_deg * M_PI / 180.0;
    for (int i = 0; i < p.N; i++) {
        s.x.push_back(p.box * U(rng));
        s.y.push_back(p.box * U(rng));
        s.th.push_back(sd * G(rng));
        s.v.push_back(p.v0);
    }
    return s;
}

// The move with the highest score; ties (within 1e-9 relative) are broken uniformly at random. Moves whose every
// future collides score -1; if all do, the agent picks among all of them.
template <class S>
int choose(const S* score, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> U(0.0, 1.0);
    double mx = score[0];
    for (int b = 1; b < NM; b++) mx = std::max(mx, (double)score[b]);
    const double tol = 1e-9 * std::max(1.0, std::fabs(mx));
    int cand[NM], nc = 0;
    for (int b = 0; b < NM; b++)
        if (score[b] >= mx - tol) cand[nc++] = b;
    return cand[std::min(nc - 1, (int)(U(rng) * nc))];
}

// Agent i takes move m (plus optional motor noise) and moves.
inline void apply_move(const Params& p, State& s, int i, int m, int step) {
    double th = s.th[i] + MOVES[m].dh * p.dth(), v = p.v0 + MOVES[m].ds * p.dv;
    if (p.noise > 0) {
        auto rng = fsm::make_rng(p.seed ^ fsm::kNoiseSeedXor, step, i);
        std::normal_distribution<double> G(0.0, 1.0);
        th += p.noise * G(rng);
        v += p.noise * G(rng);
    }
    s.th[i] = th;
    s.v[i] = v;
    s.x[i] += v * p.dt * std::cos(th);
    s.y[i] += v * p.dt * std::sin(th);
}

}  // namespace fsm2d
