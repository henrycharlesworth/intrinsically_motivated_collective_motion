// The 3D FSM model (docs/NOTES.md, "3D"): parameters, initial state, the decision rule and the motion update.
// python/fsmtools/reference3d.py is the brute-force reference implementation.
//
// Agents are spheres of radius R with a heading given by an azimuth theta and an elevation phi in a gravity-aligned
// frame (no roll). Seven moves: yaw by -+dtheta, pitch by -+dtheta, speed v0 -+ dv, or straight on at v0; a pitch past
// +-phi_max is unavailable. The visual field is the sphere of directions, seen in each future's yaw frame and sampled
// on a latitude-longitude grid; the score is the 2D one with solid-angle weights. Elevations stay on the lattice of
// multiples of dtheta, so the future tree depends only on the elevation index.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../common/cli.hpp"
#include "../common/json.hpp"
#include "../common/rng.hpp"

namespace fsm3d {

constexpr int NM = 7;  // moves: yaw -, yaw +, pitch -, pitch +, slower, straight on, faster
constexpr int STAY = 5;
constexpr double TWO_PI = 2.0 * M_PI;
struct Move { int da, de, ds; };  // azimuth, elevation (in dtheta) and speed (in dv) changes
constexpr std::array<Move, NM> MOVES{{{-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 0}, {0, 0, 1}}};

// A predator: a lethal sphere flying at the flock (docs/NOTES.md, "Predators"). The agents predict it straight on;
// their futures that come within R + radius + grow * k of it at lookahead level k die like collisions. It is not
// part of their visual fields.
struct PredatorParams {
    std::string mode = "off";   // off | line (straight) | com (turns towards the flock centre) | nearest | intercept
    int start = 200;            // step at which it appears ...
    double dist = 300;          // ... this far from the flock's centre ...
    double az = 0, el = 0;      // ... in this direction relative to the flock's heading (0, 0: ahead, degrees)
    double speed = 15, radius = 3, turn = 10, grow = 0;
    bool vision = false;        // only agents that see it know of it (a rear blind cone, and a clear line of sight)
    double blind = 60, memory = 30, range = 1e30;
    std::string mode2;          // a second predator (--pred2 MODE,AZ,EL); it shares the other parameters
    double az2 = 180, el2 = 0;
    std::vector<double> state;  // --pred-state x,y,z,vx,vy,vz[,...]: explicit state at `start` (or for --score)

    int count() const { return mode == "off" ? 0 : mode2.empty() ? 1 : 2; }
};

struct Params {
    int N = 100, tau = 3, steps = 1000;
    double v0 = 10, dv = 2, dtheta_deg = 15, R = 1, dt = 1, phi_max_deg = 75;
    bool half_asin = false;
    double box = 50, heading_sd_deg = -1;  // initial positions uniform in [0, box]^3, headings ~ N(0, sd) around +x
    uint64_t seed = 1;
    std::string out = "run3d", init;       // init: text file with "x y z theta phi v" per agent
    int log_every = 10;
    bool resume = false;
    int nlon = 384, nlat = 192;            // visual-field grid
    double tile = 24, merge = 0;           // collision tile size; position merging grid (0: exact)
    PredatorParams pred;

    double dth() const { return dtheta_deg * M_PI / 180.0; }
    int emax() const { return (int)std::floor(phi_max_deg / dtheta_deg + 1e-9); }  // largest elevation index

    bool parse(const std::string& f, fsm::Args& a) {
        if (f == "--N") N = a.integer();
        else if (f == "--tau") tau = a.integer();
        else if (f == "--steps") steps = a.integer();
        else if (f == "--v0") v0 = a.real();
        else if (f == "--dv") dv = a.real();
        else if (f == "--dtheta") dtheta_deg = a.real();
        else if (f == "--R") R = a.real();
        else if (f == "--dt") dt = a.real();
        else if (f == "--phi-max") phi_max_deg = a.real();
        else if (f == "--half-angle") half_asin = a.choice({"atan", "asin"}) == "asin";
        else if (f == "--box") box = a.real();
        else if (f == "--heading-sd") heading_sd_deg = a.real();
        else if (f == "--seed") seed = a.u64();
        else if (f == "--out") out = a.str();
        else if (f == "--init") init = a.str();
        else if (f == "--log-every") log_every = a.integer();
        else if (f == "--resume") resume = true;
        else if (f == "--grid") {
            const auto g = fsm::Args::split(a.str(), 'x');
            if (g.size() != 2) throw std::runtime_error("--grid: NLONxNLAT, e.g. 192x96");
            nlon = std::stoi(g[0]); nlat = std::stoi(g[1]);
        }
        else if (f == "--tile") tile = a.real();
        else if (f == "--merge") merge = a.real();
        else if (f == "--predator") pred.mode = a.choice({"off", "line", "com", "nearest", "intercept"});
        else if (f == "--pred-start") pred.start = a.integer();
        else if (f == "--pred-dist") pred.dist = a.real();
        else if (f == "--pred-from") {
            const auto v = a.reals();
            pred.az = v.at(0); pred.el = v.size() > 1 ? v[1] : 0;
        }
        else if (f == "--pred-speed") pred.speed = a.real();
        else if (f == "--pred-radius") pred.radius = a.real();
        else if (f == "--pred-turn") pred.turn = a.real();
        else if (f == "--pred-grow") pred.grow = a.real();
        else if (f == "--pred-vision") pred.vision = a.choice({"on", "off"}) == "on";
        else if (f == "--pred-blind") pred.blind = a.real();
        else if (f == "--pred-memory") pred.memory = a.real();
        else if (f == "--pred-range") pred.range = a.real();
        else if (f == "--pred2") {
            const auto v = fsm::Args::split(a.str(), ',');
            if (v.size() != 3) throw std::runtime_error("--pred2: MODE,AZ,EL");
            pred.mode2 = v[0]; pred.az2 = std::stod(v[1]); pred.el2 = std::stod(v[2]);
        }
        else if (f == "--pred-state") {
            pred.state = a.reals();
            if (pred.state.size() != 6 && pred.state.size() != 12) throw std::runtime_error("--pred-state: x,y,z,vx,vy,vz[,x2,..,vz2]");
        }
        else return false;
        return true;
    }

    void finish() {
        if (heading_sd_deg < 0) heading_sd_deg = dtheta_deg;
        if (!pred.state.empty() && pred.mode == "off") pred.mode = "line";            // given explicitly
        if (pred.state.size() == 12 && pred.mode2.empty()) pred.mode2 = "line";
        const std::string m2 = pred.mode2.empty() ? "line" : pred.mode2;
        if (m2 != "line" && m2 != "com" && m2 != "nearest" && m2 != "intercept") throw std::runtime_error("--pred2: line, com, nearest or intercept");
        if (!pred.mode2.empty() && pred.mode == "off") throw std::runtime_error("--pred2 needs --predator");
        const double per = nlon * dtheta_deg / 360.0;
        if (std::fabs(per - std::round(per)) > 1e-9) throw std::runtime_error("--grid: nlon must be a multiple of 360 / dtheta");
        if (nlon % 32) throw std::runtime_error("--grid: nlon must be a multiple of 32");
        if (N < 2 || tau < 1 || steps < 0) throw std::runtime_error("need --N >= 2, --tau >= 1, --steps >= 0");
    }

    void describe(fsm::Json& j, int n_agents) const {
        j.num("dims", 3).num("N", n_agents).num("tau", tau)
            .num("v0", v0).num("dv", dv).num("dtheta_deg", dtheta_deg).num("R", R).num("dt", dt).num("phi_max_deg", phi_max_deg)
            .str("half_angle", half_asin ? "asin" : "atan").num("box", box).num("heading_sd_deg", heading_sd_deg)
            .num("seed", (unsigned long long)seed).str("init", init)
            .nums("grid", {(double)nlon, (double)nlat}).num("tile", tile).num("merge", merge);
        if (pred.mode == "off") return;
        fsm::Json q;
        q.str("mode", pred.mode).num("start", pred.start).num("speed", pred.speed).num("radius", pred.radius).num("grow", pred.grow)
            .boolean("vision", pred.vision).num("blind_deg", pred.blind).num("memory", pred.memory)
            .num("range", pred.range > 1e29 ? -1.0 : pred.range).num("dist", pred.dist).num("turn_deg", pred.turn)
            .nums("from_deg", {pred.az, pred.el}).num("n", pred.count());
        if (!pred.mode2.empty()) {
            fsm::Json s;
            s.str("mode", pred.mode2).nums("from_deg", {pred.az2, pred.el2});
            q.obj("second", s);
        }
        q.str("file", "float64 [steps+1, n, 7] (x, y, z, vx, vy, vz, caught so far; NaN before start)");
        if (pred.vision) q.str("seen", "uint8 [steps, N]: per predator r a base-3 digit (x 3^r): 0 unaware, 1 sees it, 2 remembers it");
        j.obj("predator", q);
    }
};

// The state: [N][6] doubles (x, y, z, theta, phi, v), world frame; also the .traj frame layout.
using State = std::vector<double>;

inline State read_state(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    State X;
    std::string line;
    while (std::getline(f, line)) {
        std::istringstream ls(line);
        double v[6];
        if (!(ls >> v[0] >> v[1] >> v[2] >> v[3] >> v[4] >> v[5])) continue;
        X.insert(X.end(), v, v + 6);
    }
    return X;
}

// Random positions in a cube, headings around +x (elevations rounded to the lattice), or the --init file.
inline State initial_state(const Params& p) {
    if (!p.init.empty()) return read_state(p.init);
    auto rng = fsm::make_rng(p.seed, fsm::kInitStream, 0);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    std::normal_distribution<double> G(0.0, 1.0);
    const double sd = p.heading_sd_deg * M_PI / 180.0, dth = p.dth();
    const int EM = p.emax();
    State X(6 * (size_t)p.N);
    for (int i = 0; i < p.N; i++) {
        X[6 * i] = p.box * U(rng); X[6 * i + 1] = p.box * U(rng); X[6 * i + 2] = p.box * U(rng);
        X[6 * i + 3] = sd * G(rng);
        X[6 * i + 4] = dth * std::clamp((double)std::lround(sd * G(rng) / dth), (double)-EM, (double)EM);
        X[6 * i + 5] = p.v0;
    }
    return X;
}

// The best move, ties broken uniformly at random (as in 2D). When every move dies, stay (a pitch move could be past
// the elevation cap).
inline int choose(const float* sc, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> U(0.0, 1.0);
    double mx = sc[0];
    for (int b = 1; b < NM; b++) mx = std::max(mx, (double)sc[b]);
    if (mx < 0) return STAY;
    const double tol = 1e-9 * std::max(1.0, std::fabs(mx));
    int cand[NM], nc = 0;
    for (int b = 0; b < NM; b++)
        if (sc[b] >= mx - tol) cand[nc++] = b;
    return cand[std::min(nc - 1, (int)(U(rng) * nc))];
}

// Agent i takes move m; returns its new unit heading.
inline std::array<double, 3> apply_move(const Params& p, State& X, int i, int m) {
    double* x = &X[6 * (size_t)i];
    const double dth = p.dth();
    x[3] += MOVES[m].da * dth;
    x[4] = std::lround(x[4] / dth + MOVES[m].de) * dth;  // stays exactly on the lattice
    x[5] = p.v0 + MOVES[m].ds * p.dv;
    const std::array<double, 3> u{std::cos(x[4]) * std::cos(x[3]), std::cos(x[4]) * std::sin(x[3]), std::sin(x[4])};
    x[0] += x[5] * p.dt * u[0]; x[1] += x[5] * p.dt * u[1]; x[2] += x[5] * p.dt * u[2];
    return u;
}

}  // namespace fsm3d
