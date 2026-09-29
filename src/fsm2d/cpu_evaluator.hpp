// Exact CPU evaluation of one agent's decision (docs/OPTIMIZATIONS.md, "CPU").
//
// For agent i: place the DAG in its frame; prune the edges that collide with the other agents' ballistic futures,
// counting for each node and first move b the collision-free paths into it; build the visual state of every
// reachable node (the union of the angular intervals subtended by the other agents, relative to the node's
// heading); and score every first move b with
//     W_b = (1 / n_b) sum_{nodes u, v} m_b(u) m_b(v) d(u, v),   d = measure(view_u XOR view_v) / 2 pi
// (m_b: path counts, n_b = sum_u m_b(u)), computed in O(n log n) instead of O(n^2) as
//     W_b = (1 / (pi n_b)) integral S_b(theta) (n_b - S_b(theta)) dtheta,   S_b(theta) = sum_{u sees theta} m_b(u)
// (docs/NOTES.md derives this). Everything is in double precision and exact up to rounding.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>
#include <vector>

#include "dag.hpp"
#include "model.hpp"

namespace fsm2d {

// Monotone map from a direction (x, y) to [0, 4), ordered like atan2 mapped to [0, 2 pi): sorting by it needs no
// trigonometry. pseudo_to_angle inverts it (used only on the merged interval endpoints).
inline double pseudo_angle(double x, double y) {
    const double a = x / (std::fabs(x) + std::fabs(y));
    return (y >= 0) ? 1.0 - a : 3.0 + a;
}
inline double pseudo_to_angle(double p) {
    if (p <= 0) return 0;
    if (p >= 4) return TWO_PI;
    if (p <= 2) {
        const double a = 1 - p;
        return std::atan2(1 - std::fabs(a), a);
    }
    const double a = p - 3;
    return std::atan2(-(1 - std::fabs(a)), a) + TWO_PI;
}

// Exact sort by a key in [0, kmax): a counting sort into n buckets, then an insertion-sort pass (only elements of the
// same bucket can be out of order), so O(n) for keys that are reasonably spread out.
template <class T, class K>
void bucket_sort(std::vector<T>& a, std::vector<T>& tmp, std::vector<int>& bkt, std::vector<int>& cnt, double kmax, K key) {
    const int n = (int)a.size();
    if (n < 32) {
        std::sort(a.begin(), a.end(), [&](const T& x, const T& y) { return key(x) < key(y); });
        return;
    }
    const double scale = n / kmax;
    cnt.assign(n + 1, 0);
    bkt.resize(n);
    for (int q = 0; q < n; q++) {
        const int b = std::min(n - 1, std::max(0, (int)(key(a[q]) * scale)));
        bkt[q] = b;
        cnt[b + 1]++;
    }
    for (int b = 0; b < n; b++) cnt[b + 1] += cnt[b];
    tmp.resize(n);
    for (int q = 0; q < n; q++) tmp[cnt[bkt[q]]++] = a[q];
    for (int q = 1; q < n; q++) {
        const T x = tmp[q];
        const double kx = key(x);
        int r = q;
        while (r > 0 && key(tmp[r - 1]) > kx) {
            tmp[r] = tmp[r - 1];
            r--;
        }
        tmp[r] = x;
    }
    a.swap(tmp);
}

using Scores = std::array<double, NM>;  // -1: every future of that first move collides

// Per-thread scratch space, reused across agents.
struct Workspace {
    struct Event {
        double th;
        int node, sign;
    };
    std::vector<double> px, py, ch, sh;
    std::vector<std::array<int32_t, NM>> mult;  // [node][first move]: collision-free paths into the node
    std::vector<std::vector<int>> cand;         // [level]: collision candidates
    std::vector<std::pair<double, double>> ivs, ivs_tmp;
    std::vector<int> bkt, cnt;
    std::vector<double> vis;                    // merged intervals (true angles, node frame) of all nodes, flat
    std::vector<int> voff, vcnt;
    std::vector<Event> ev, ev_tmp;
    uint64_t nodes = 0, edges = 0;              // work counters (since last read)
};

class CpuEvaluator {
public:
    CpuEvaluator(const Params& p, bool merge_dag = true) : p_(p), dag_(build_dag(p, merge_dag)), D2_(4 * p.R * p.R) {}

    const Dag& dag() const { return dag_; }

    // The other agents' ballistic futures (straight on at v0), once per step for all agents: level k of agent j at
    // [k * N + j].
    void set_state(const State& s) {
        s_ = &s;
        const int N = s.size();
        QX_.assign((size_t)(p_.tau + 1) * N, 0);
        QY_.assign(QX_.size(), 0);
        for (int j = 0; j < N; j++) {
            const double cx = std::cos(s.th[j]) * p_.v0 * p_.dt, cy = std::sin(s.th[j]) * p_.v0 * p_.dt;
            for (int k = 0; k <= p_.tau; k++) {
                QX_[(size_t)k * N + j] = s.x[j] + k * cx;
                QY_[(size_t)k * N + j] = s.y[j] + k * cy;
            }
        }
    }

    Scores evaluate(Workspace& ws, int i) const {
        const Dag& t = dag_;
        const int N = s_->size(), T = p_.tau, nn = t.size();
        const double x0 = s_->x[i], y0 = s_->y[i], th0 = s_->th[i], dth = p_.dth();
        const double c0 = std::cos(th0), s0 = std::sin(th0);
        ws.ch.resize(2 * T + 1);
        ws.sh.resize(2 * T + 1);
        for (int hh = -T; hh <= T; hh++) {
            ws.ch[hh + T] = std::cos(th0 + hh * dth);
            ws.sh[hh + T] = std::sin(th0 + hh * dth);
        }
        ws.px.resize(nn);
        ws.py.resize(nn);
        for (int n = 0; n < nn; n++) {
            ws.px[n] = x0 + c0 * t.x[n] - s0 * t.y[n];
            ws.py[n] = y0 + s0 * t.x[n] + c0 * t.y[n];
        }
        // Collision candidates per level: other agents whose step segment's bounding box meets that of the nodes at
        // levels k-1 and k (inflated by 2R).
        ws.cand.resize(T + 1);
        for (int k = 1; k <= T; k++) {
            double bx0 = 1e300, bx1 = -1e300, by0 = 1e300, by1 = -1e300;
            for (int n = t.lstart[k - 1]; n < t.lstart[k + 1]; n++) {
                bx0 = std::min(bx0, ws.px[n]); bx1 = std::max(bx1, ws.px[n]);
                by0 = std::min(by0, ws.py[n]); by1 = std::max(by1, ws.py[n]);
            }
            const double m = 2 * p_.R;
            bx0 -= m; bx1 += m; by0 -= m; by1 += m;
            auto& cl = ws.cand[k];
            cl.clear();
            const double* ax = &QX_[(size_t)(k - 1) * N];
            const double* ay = &QY_[(size_t)(k - 1) * N];
            const double* bx = &QX_[(size_t)k * N];
            const double* by = &QY_[(size_t)k * N];
            for (int j = 0; j < N; j++) {
                if (j == i) continue;
                if (std::max(ax[j], bx[j]) < bx0 || std::min(ax[j], bx[j]) > bx1) continue;
                if (std::max(ay[j], by[j]) < by0 || std::min(ay[j], by[j]) > by1) continue;
                cl.push_back(j);
            }
        }
        // Collision-free path counts per first move.
        ws.mult.assign(nn, {0, 0, 0, 0, 0});
        for (int e = t.estart[1]; e < t.estart[2]; e++) {
            const int c = t.echild[e];
            ws.edges++;
            if (!edge_collides(ws, 1, x0, y0, ws.px[c], ws.py[c])) ws.mult[c][t.emove[e]] = 1;
        }
        for (int k = 2; k <= T; k++)
            for (int e = t.estart[k]; e < t.estart[k + 1]; e++) {
                const int a = t.eparent[e], c = t.echild[e];
                const auto& ma = ws.mult[a];
                if ((ma[0] | ma[1] | ma[2] | ma[3] | ma[4]) == 0) continue;
                ws.edges++;
                if (edge_collides(ws, k, ws.px[a], ws.py[a], ws.px[c], ws.py[c])) continue;
                auto& mc = ws.mult[c];
                for (int b = 0; b < NM; b++) mc[b] += ma[b];
            }
        // Visual states of the reachable nodes.
        ws.vis.clear();
        ws.voff.assign(nn, 0);
        ws.vcnt.assign(nn, 0);
        std::array<double, NM> nb{};
        for (int k = 1; k <= T; k++)
            for (int n = t.lstart[k]; n < t.lstart[k + 1]; n++) {
                const auto& m = ws.mult[n];
                if ((m[0] | m[1] | m[2] | m[3] | m[4]) == 0) continue;
                for (int b = 0; b < NM; b++) nb[b] += m[b];
                ws.voff[n] = (int)ws.vis.size();
                const int hh = t.h[n] + T;
                visual_state(ws, i, k, ws.px[n], ws.py[n], ws.ch[hh], ws.sh[hh], ws.vis);
                ws.vcnt[n] = ((int)ws.vis.size() - ws.voff[n]) / 2;
                ws.nodes++;
            }
        // One sweep over all interval endpoints: the integral of S_b (n_b - S_b) for all first moves at once.
        Scores score;
        for (int b = 0; b < NM; b++) score[b] = nb[b] > 0 ? 0.0 : -1.0;
        auto& ev = ws.ev;
        ev.clear();
        for (int n = t.lstart[1]; n < nn; n++)
            for (int q = 0; q < ws.vcnt[n]; q++) {
                const double* iv = &ws.vis[ws.voff[n] + 2 * q];
                ev.push_back({iv[0], n, +1});
                ev.push_back({iv[1], n, -1});
            }
        bucket_sort(ev, ws.ev_tmp, ws.bkt, ws.cnt, TWO_PI, [](const Workspace::Event& e) { return e.th; });
        std::array<double, NM> S{}, acc{};
        double prev = 0;
        for (const auto& e : ev) {
            const double d = e.th - prev;
            if (d > 0) {
                for (int b = 0; b < NM; b++) acc[b] += S[b] * (nb[b] - S[b]) * d;
                prev = e.th;
            }
            const auto& m = ws.mult[e.node];
            for (int b = 0; b < NM; b++) S[b] += e.sign * m[b];
        }
        for (int b = 0; b < NM; b++)
            if (nb[b] > 0) score[b] = acc[b] / (M_PI * nb[b]);
        return score;
    }

    // The same scores by brute force, from the definition (all pairs of futures, O(n^2)), using the workspace left by
    // evaluate(). For --verify only.
    Scores brute_force(const Workspace& ws) const {
        auto overlap = [&](int a, int b) {
            const double* A = &ws.vis[ws.voff[a]];
            const double* B = &ws.vis[ws.voff[b]];
            int ia = 0, ib = 0;
            double tot = 0;
            while (ia < ws.vcnt[a] && ib < ws.vcnt[b]) {
                const double lo = std::max(A[2 * ia], B[2 * ib]), hi = std::min(A[2 * ia + 1], B[2 * ib + 1]);
                if (hi > lo) tot += hi - lo;
                if (A[2 * ia + 1] < B[2 * ib + 1]) ia++; else ib++;
            }
            return tot;
        };
        auto measure = [&](int a) {
            double tot = 0;
            for (int q = 0; q < ws.vcnt[a]; q++) tot += ws.vis[ws.voff[a] + 2 * q + 1] - ws.vis[ws.voff[a] + 2 * q];
            return tot;
        };
        const int nn = dag_.size();
        Scores out;
        for (int b = 0; b < NM; b++) {
            double tot = 0, n = 0;
            for (int a = dag_.lstart[1]; a < nn; a++) {
                if (ws.mult[a][b] == 0) continue;
                n += ws.mult[a][b];
                for (int c = dag_.lstart[1]; c < nn; c++) {
                    if (ws.mult[c][b] == 0) continue;
                    const double dij = (measure(a) + measure(c) - 2 * overlap(a, c)) / TWO_PI;
                    tot += (double)ws.mult[a][b] * ws.mult[c][b] * dij;
                }
            }
            out[b] = n > 0 ? tot / n : -1;
        }
        return out;
    }

private:
    // Does the agent, moving p0 -> p1 during step k (level k-1 -> k), collide with any candidate (discs of radius R
    // overlapping at the end of the step, or passing through each other during it)?
    bool edge_collides(const Workspace& ws, int k, double p0x, double p0y, double p1x, double p1y) const {
        const int N = s_->size();
        const double* q0x = &QX_[(size_t)(k - 1) * N];
        const double* q0y = &QY_[(size_t)(k - 1) * N];
        const double* q1x = &QX_[(size_t)k * N];
        const double* q1y = &QY_[(size_t)k * N];
        for (int j : ws.cand[k]) {
            const double r1x = p1x - q1x[j], r1y = p1y - q1y[j];
            if (r1x * r1x + r1y * r1y < D2_) return true;
            const double r0x = p0x - q0x[j], r0y = p0y - q0y[j];
            if (r0x * r0x + r0y * r0y < D2_) continue;  // already overlapping at the start: ignore the pair
            const double dx = r1x - r0x, dy = r1y - r0y, dd = dx * dx + dy * dy;
            if (dd <= 0) continue;
            const double s = -(r0x * dx + r0y * dy) / dd;  // closest approach along the step
            if (s > 0 && s < 1) {
                const double mx = r0x + s * dx, my = r0y + s * dy;
                if (mx * mx + my * my < D2_) return true;
            }
        }
        return false;
    }

    // The merged occupied angular intervals (true angles in [0, 2 pi], relative to the heading (c, s)) seen from
    // (px, py) at level k, appended to out as (lo, hi) pairs.
    void visual_state(Workspace& ws, int i, int k, double px, double py, double c, double s, std::vector<double>& out) const {
        const int N = s_->size();
        const double* qx = &QX_[(size_t)k * N];
        const double* qy = &QY_[(size_t)k * N];
        const double R = p_.R, R2 = R * R;
        auto& ivs = ws.ivs;
        ivs.clear();
        for (int j = 0; j < N; j++) {
            if (j == i) continue;
            const double rx = qx[j] - px, ry = qy[j] - py;
            const double u = c * rx + s * ry, v = -s * rx + c * ry;  // node frame: u ahead, v to the left
            const double d2 = u * u + v * v;
            // The tangent directions are (u, v) rotated by -+alpha; scaled by d^2 / R they are (u a -+ v R, +-u R + v a)
            // with a = sqrt(d^2 - R^2) (asin rule) or a = d (atan rule, as in the paper).
            double a;
            if (p_.half_asin) {
                if (d2 <= R2) {  // inside another disc: the whole field is occupied
                    out.push_back(0);
                    out.push_back(TWO_PI);
                    return;
                }
                a = std::sqrt(d2 - R2);
            } else {
                if (d2 <= 0) continue;
                a = std::sqrt(d2);
            }
            const double lo = pseudo_angle(u * a + v * R, -u * R + v * a);
            const double hi = pseudo_angle(u * a - v * R, u * R + v * a);
            if (lo <= hi) {
                ivs.emplace_back(lo, hi);
            } else {  // wraps through angle 0
                ivs.emplace_back(lo, 4.0);
                ivs.emplace_back(0.0, hi);
            }
        }
        if (ivs.empty()) return;
        bucket_sort(ivs, ws.ivs_tmp, ws.bkt, ws.cnt, 4.0, [](const std::pair<double, double>& x) { return x.first; });
        double clo = ivs[0].first, chi = ivs[0].second;
        for (size_t q = 1; q < ivs.size(); q++) {
            if (ivs[q].first <= chi) {
                if (ivs[q].second > chi) chi = ivs[q].second;
            } else {
                out.push_back(pseudo_to_angle(clo));
                out.push_back(pseudo_to_angle(chi));
                clo = ivs[q].first;
                chi = ivs[q].second;
            }
        }
        out.push_back(pseudo_to_angle(clo));
        out.push_back(pseudo_to_angle(chi));
    }

    Params p_;
    Dag dag_;
    double D2_;  // (2R)^2
    const State* s_ = nullptr;
    std::vector<double> QX_, QY_;
};

}  // namespace fsm2d
