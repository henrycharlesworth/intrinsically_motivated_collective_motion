// Predators (3D): their motion, what the agents know of them, and catches (docs/NOTES.md, "Predators").
//
// Modes: "line" keeps its initial velocity; "com" / "nearest" turn it by at most --pred-turn degrees per step towards
// the flock's centre of mass / the nearest agent; "intercept" towards the point where the flock will be when the
// predator could get there. The agents predict it straight on at its current velocity; for their futures it is
// lethal within R + radius + grow * k at lookahead level k (a growing zone for their uncertainty about where it
// goes). Catches are counted (the agents are not removed) with R + radius.
//
// With --pred-vision on, an agent knows of a predator only if it sees it (outside its rear blind cone, within range,
// with a clear line of sight) and remembers the last sighting for --pred-memory steps, extrapolated straight on, with
// the danger zone grown by grow per step since. Unaware agents do plain FSM.
#pragma once

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "model.hpp"

namespace fsm3d {

// A predator appears at step --pred-start, --pred-dist from the flock's centre of mass in the direction (az, el)
// relative to the flock's mean heading (0, 0: ahead of the flock, flying at it), heading for the centre of mass at
// --pred-speed.
inline std::vector<double> predator_start(const Params& p, const std::vector<double>& X, int N, double az, double el) {
    double c[3] = {0, 0, 0}, h[3] = {0, 0, 0};
    for (int i = 0; i < N; i++) {
        const double* x = &X[6 * (size_t)i];
        for (int d = 0; d < 3; d++) c[d] += x[d] / N;
        h[0] += std::cos(x[4]) * std::cos(x[3]); h[1] += std::cos(x[4]) * std::sin(x[3]); h[2] += std::sin(x[4]);
    }
    const double th = std::atan2(h[1], h[0]) + az * M_PI / 180.0;
    const double ph = std::atan2(h[2], std::hypot(h[0], h[1])) + el * M_PI / 180.0;
    const double e[3] = {std::cos(ph) * std::cos(th), std::cos(ph) * std::sin(th), std::sin(ph)};
    return {c[0] + p.pred.dist * e[0], c[1] + p.pred.dist * e[1], c[2] + p.pred.dist * e[2],
            -p.pred.speed * e[0], -p.pred.speed * e[1], -p.pred.speed * e[2]};
}

// Which agents see the predator P (state X): the direction to it is outside the agent's rear blind cone and within
// range, and the sight line to its centre passes no other agent (sphere of radius R) before it. O(N^2), threaded.
inline void predator_seen(const Params& p, const std::vector<double>& X, int N, const double* P, std::vector<uint8_t>& seen) {
    seen.assign(N, 0);
    const double cblind = std::cos(M_PI - 0.5 * p.pred.blind * M_PI / 180.0), R2 = p.R * p.R;
    auto work = [&](int i0, int i1) {
        for (int i = i0; i < i1; i++) {
            const double* x = &X[6 * (size_t)i];
            const double d[3] = {P[0] - x[0], P[1] - x[1], P[2] - x[2]};
            const double L = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
            if (L > p.pred.range || L <= 0) { seen[i] = L <= 0; continue; }
            const double u[3] = {std::cos(x[4]) * std::cos(x[3]), std::cos(x[4]) * std::sin(x[3]), std::sin(x[4])};
            if ((u[0] * d[0] + u[1] * d[1] + u[2] * d[2]) / L < cblind) continue;  // in the blind cone behind
            const double e[3] = {d[0] / L, d[1] / L, d[2] / L};
            bool blocked = false;
            for (int j = 0; j < N && !blocked; j++) {
                if (j == i) continue;
                const double r[3] = {X[6 * (size_t)j] - x[0], X[6 * (size_t)j + 1] - x[1], X[6 * (size_t)j + 2] - x[2]};
                const double t = r[0] * e[0] + r[1] * e[1] + r[2] * e[2];
                if (t <= 0 || t >= L) continue;
                blocked = r[0] * r[0] + r[1] * r[1] + r[2] * r[2] - t * t < R2;
            }
            seen[i] = !blocked;
        }
    };
    const int nt = std::max(1, std::min(16, (int)std::thread::hardware_concurrency()));
    std::vector<std::thread> th;
    for (int t = 0; t < nt; t++) th.emplace_back(work, (int)((long long)N * t / nt), (int)((long long)N * (t + 1) / nt));
    for (auto& t : th) t.join();
}

// Each agent's knowledge of the predator from sightings: know [N][7] (state projected to `step`, age; -1: none);
// last [N]: step of the last sighting (-1: never), with the predator frames PF(step) giving the state then.
template <class F>
inline void predator_knowledge(const Params& p, int N, int step, const std::vector<int>& last, F frame, std::vector<double>& know) {
    know.assign(7 * (size_t)N, 0.0);
    for (int i = 0; i < N; i++) {
        double* q = &know[7 * (size_t)i];
        const int age = last[i] < 0 ? -1 : step - last[i];
        if (age < 0 || age > p.pred.memory) { q[6] = -1; continue; }
        const double* s = frame(last[i]);
        for (int d = 0; d < 3; d++) { q[d] = s[d] + age * p.dt * s[3 + d]; q[3 + d] = s[3 + d]; }
        q[6] = age;
    }
}

// Turn the predator's velocity towards its target (agents' positions X, before this step's move) and move it.
// intercept: the flock's centre where it will be when the predator could reach it (moving straight on at the
// flock's mean velocity; the earliest time at which the predator, at pred_speed, can be there).
inline void predator_step(const Params& p, const std::vector<double>& X, int N, const std::string& mode, double* P) {
    if (mode == "com" || mode == "nearest" || mode == "intercept") {
        double t[3] = {0, 0, 0};
        if (mode == "com" || mode == "intercept") {
            for (int i = 0; i < N; i++) for (int d = 0; d < 3; d++) t[d] += X[6 * (size_t)i + d] / N;
            if (mode == "intercept") {
                double v[3] = {0, 0, 0};
                for (int i = 0; i < N; i++) {
                    const double* x = &X[6 * (size_t)i];
                    v[0] += x[5] * std::cos(x[4]) * std::cos(x[3]) / N; v[1] += x[5] * std::cos(x[4]) * std::sin(x[3]) / N; v[2] += x[5] * std::sin(x[4]) / N;
                }
                // |t + v s - P| = pred_speed s: a s^2 + b s + c = 0 with a = |v|^2 - vp^2 < 0 (the predator is faster)
                const double r[3] = {t[0] - P[0], t[1] - P[1], t[2] - P[2]};
                const double a = v[0] * v[0] + v[1] * v[1] + v[2] * v[2] - p.pred.speed * p.pred.speed;
                const double b = 2 * (r[0] * v[0] + r[1] * v[1] + r[2] * v[2]), c = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
                const double disc = b * b - 4 * a * c;
                if (a < 0 && disc >= 0) {
                    const double sT = (-b - std::sqrt(disc)) / (2 * a);
                    if (sT > 0) for (int d = 0; d < 3; d++) t[d] += v[d] * sT;
                }
            }
        } else {
            double best = 1e300;
            for (int i = 0; i < N; i++) {
                const double* x = &X[6 * (size_t)i];
                const double d2 = (x[0] - P[0]) * (x[0] - P[0]) + (x[1] - P[1]) * (x[1] - P[1]) + (x[2] - P[2]) * (x[2] - P[2]);
                if (d2 < best) { best = d2; t[0] = x[0]; t[1] = x[1]; t[2] = x[2]; }
            }
        }
        double u[3] = {P[3], P[4], P[5]}, g[3] = {t[0] - P[0], t[1] - P[1], t[2] - P[2]};
        const double nu = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]), ng = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
        if (nu > 0 && ng > 0) {
            for (int d = 0; d < 3; d++) { u[d] /= nu; g[d] /= ng; }
            const double cg = std::clamp(u[0] * g[0] + u[1] * g[1] + u[2] * g[2], -1.0, 1.0), ang = std::acos(cg);
            const double turn = std::min(ang, p.pred.turn * M_PI / 180.0);
            double w[3] = {g[0] - cg * u[0], g[1] - cg * u[1], g[2] - cg * u[2]};
            double nw = std::sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
            if (nw < 1e-12) { w[0] = -u[1]; w[1] = u[0]; w[2] = 0; nw = std::hypot(w[0], w[1]); if (nw < 1e-12) { w[0] = 1; nw = 1; } }
            for (int d = 0; d < 3; d++) P[3 + d] = p.pred.speed * (std::cos(turn) * u[d] + std::sin(turn) * w[d] / nw);
        }
    }
    for (int d = 0; d < 3; d++) P[d] += p.dt * P[3 + d];
}

// Agents that came within R + pred_radius of the predator during the step (both moving linearly from X0 / P0).
inline int predator_catches(const Params& p, const std::vector<double>& X0, const std::vector<double>& X1,
                            const double* P0, const double* P1, int N) {
    const double D = p.R + p.pred.radius;
    int n = 0;
    for (int i = 0; i < N; i++) {
        double r0[3], dr[3];
        for (int d = 0; d < 3; d++) { r0[d] = X0[6 * (size_t)i + d] - P0[d]; dr[d] = (X1[6 * (size_t)i + d] - P1[d]) - r0[d]; }
        const double dd = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
        const double s = dd > 0 ? std::clamp(-(r0[0] * dr[0] + r0[1] * dr[1] + r0[2] * dr[2]) / dd, 0.0, 1.0) : 0.0;
        double m2 = 0;
        for (int d = 0; d < 3; d++) m2 += (r0[d] + s * dr[d]) * (r0[d] + s * dr[d]);
        n += m2 < D * D;
    }
    return n;
}

// The predators of a run and the agents' knowledge of them, step by step. Output: PREFIX.pred (one record per
// frame: [n][7] x, y, z, vx, vy, vz, agents caught so far; NaN before the start) and, with vision, PREFIX.seen (one
// record per step: per agent the base-3 digits of what it knows of each predator).
class Predators {
public:
    std::vector<double> state;  // [n][6], world frame; empty: not there yet
    std::vector<double> know;   // [n][N][7] each agent's knowledge (state projected to now, age; -1: none); empty: all know
    std::vector<uint8_t> flags; // [N] this step's sightings (base-3 digits)
    long long caught = 0;

    Predators(const Params& p, int N) : p_(p), n_(p.pred.count()), N_(N), flags(N, 0) {
        if (p_.pred.vision) {
            hist_.assign(6 * (size_t)n_ * (p_.steps + 1), std::nan(""));
            last_.assign((size_t)std::max(1, n_) * N_, -1);
        }
    }
    int count() const { return n_; }
    bool present() const { return !state.empty(); }
    size_t record_bytes() const { return sizeof(double) * 7 * n_; }

    // The .pred record of the current state.
    std::vector<double> record() const {
        std::vector<double> r(7 * (size_t)n_, std::nan(""));
        if (present())
            for (int q = 0; q < n_; q++) {
                for (int d = 0; d < 6; d++) r[7 * q + d] = state[6 * q + d];
                r[7 * q + 6] = (double)caught;
            }
        return r;
    }

    // Resume at `start`: the state from the .pred record of that frame, and with vision the sighting history.
    void restore(const std::string& prefix, int start) {
        std::ifstream f(prefix + ".pred", std::ios::binary);
        std::vector<double> r(7 * (size_t)n_);
        for (int q = 0; q <= start && f.read((char*)r.data(), (std::streamsize)record_bytes()); q++) {
            if (p_.pred.vision)
                for (int k = 0; k < n_; k++) for (int d = 0; d < 6; d++) hist_[6 * ((size_t)q * n_ + k) + d] = r[7 * k + d];
            if (q == start && !std::isnan(r[0])) {
                state.clear();
                for (int k = 0; k < n_; k++) state.insert(state.end(), r.begin() + 7 * k, r.begin() + 7 * k + 6);
                caught = (long long)r[6];
            }
        }
        if (!p_.pred.vision) return;
        std::ifstream fs(prefix + ".seen", std::ios::binary);
        std::vector<uint8_t> row(N_);
        for (int q = 0; q < start && fs.read((char*)row.data(), N_); q++)
            for (int i = 0; i < N_; i++)
                for (int k = 0, v = row[i]; k < n_; k++, v /= 3) if (v % 3 == 1) last_[(size_t)k * N_ + i] = q;
    }

    // The predators appear at frame --pred-start (placed from the flock's state X at that frame). Called for every
    // frame before its .pred record is written, so the first frame with the predators is recorded too.
    void appear(int frame, const State& X) {
        if (present() || frame < p_.pred.start) return;
        const auto& given = p_.pred.state;
        if (!given.empty()) state.assign(given.begin(), given.begin() + std::min<size_t>(given.size(), 6 * (size_t)n_));
        else state = predator_start(p_, X, N_, p_.pred.az, p_.pred.el);
        if (n_ == 2 && state.size() < 12) {
            const auto q = predator_start(p_, X, N_, p_.pred.az2, p_.pred.el2);
            state.insert(state.end(), q.begin(), q.end());
        }
    }

    // Before the agents decide at `step`: with vision, who sees the predators now and what each agent knows of them.
    void before_step(int step, const State& X) {
        if (!p_.pred.vision) return;
        std::fill(flags.begin(), flags.end(), 0);
        if (!present()) return;
        know.clear();
        std::vector<int> last(N_);
        std::vector<double> kr;
        std::vector<uint8_t> seen;
        for (int k = 0, pw = 1; k < n_; k++, pw *= 3) {
            for (int d = 0; d < 6; d++) hist_[6 * ((size_t)step * n_ + k) + d] = state[6 * k + d];
            predator_seen(p_, X, N_, &state[6 * k], seen);
            for (int i = 0; i < N_; i++) {
                if (seen[i]) last_[(size_t)k * N_ + i] = step;
                last[i] = last_[(size_t)k * N_ + i];
            }
            predator_knowledge(p_, N_, step, last, [&](int q) { return &hist_[6 * ((size_t)q * n_ + k)]; }, kr);
            know.insert(know.end(), kr.begin(), kr.end());
            for (int i = 0; i < N_; i++) flags[i] += (uint8_t)(pw * (seen[i] ? 1 : kr[7 * (size_t)i + 6] >= 0 ? 2 : 0));
        }
    }

    // After the agents moved from X0 to X: the predators move; returns the agents caught during the step.
    int after_step(const State& X0, const State& X) {
        if (!present()) return 0;
        const std::vector<double> P0 = state;
        const std::string modes[2] = {p_.pred.mode, p_.pred.mode2};
        int c = 0;
        for (int k = 0; k < n_; k++) {
            predator_step(p_, X0, N_, modes[k], &state[6 * k]);
            c += predator_catches(p_, X0, X, &P0[6 * k], &state[6 * k], N_);
        }
        caught += c;
        return c;
    }

    // For the progress line: agents that see any predator / know of one (with vision), distance of the first from
    // the flock's centre.
    std::string summary(const State& X, int caught_step) const {
        if (!present()) return "";
        double c[3] = {0, 0, 0};
        for (int i = 0; i < N_; i++) for (int d = 0; d < 3; d++) c[d] += X[6 * (size_t)i + d] / N_;
        char buf[160];
        std::snprintf(buf, sizeof buf, "  predator %.0f from the centre, caught %d (%lld)",
                      std::sqrt((c[0] - state[0]) * (c[0] - state[0]) + (c[1] - state[1]) * (c[1] - state[1]) + (c[2] - state[2]) * (c[2] - state[2])),
                      caught_step, caught);
        std::string s = buf;
        if (p_.pred.vision) {
            int ns = 0, nk = 0;
            for (int i = 0; i < N_; i++) {
                bool any = false;
                for (int k = 0, v = flags[i]; k < n_; k++, v /= 3) any |= v % 3 == 1;
                ns += any; nk += flags[i] != 0;
            }
            s += ", sees " + std::to_string(ns) + ", knows " + std::to_string(nk);
        }
        return s;
    }

private:
    const Params& p_;
    int n_, N_;
    std::vector<double> hist_;  // [step][n][6] states, for the agents' memories
    std::vector<int> last_;     // [n][N] step of each agent's last sighting (-1: never)
};

}  // namespace fsm3d
