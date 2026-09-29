// The 3D future-state DAG, one template per elevation index e0 (docs/OPTIMIZATIONS.md, "GPU 3D"). In the agent's yaw
// frame (origin at the agent) a node's position depends only on the headings visited, so as in 2D, move sequences
// with the same (azimuth, elevation) cells, step counts and speed-change sums merge into one node.
//   node      (level, azimuth index a, elevation index e, exact position): path counts live here
//   position  nodes at the same spot (or within --merge): one visual mask
//   unit      the nodes of a position with the same azimuth: the same view up to a longitude shift of a * per pixels
//   tile      spatial cell of positions at one level; collision candidates are filtered per tile
// Nodes are numbered level, tile, position, azimuth-major (contiguous ranges); node 0 is the root.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <thread>
#include <unordered_map>
#include <vector>

#include "model.hpp"

namespace fsm3d {

struct Template3 {
    int tau = 0, e0 = 0, M = 0, P = 0, U = 0, NT = 0;
    size_t tree_nodes = 0;
    std::vector<int> ltile, lpos;                  // [tau+2] tile / position ranges per level
    std::vector<float> node_x, node_y, node_z;
    std::vector<int> node_unit;
    std::vector<int> unit_sh;                      // azimuth shift of a unit, in longitude pixels [0, nlon)
    std::vector<int> pos_u0;                       // [P+1] unit range of each position
    std::vector<int> pos_n0;                       // [P+1] node range of each position
    std::vector<float> pos_x, pos_y, pos_z;
    std::vector<int> tile_estart;                  // [NT+1] incoming edges of each tile
    std::vector<float> tile_bb;                    // [6 NT] bbox of each tile's nodes and their parents
    std::vector<int> eparent, echild, emove;
    std::vector<float> lvl_bb;                     // [6 (tau+1)] bbox of the nodes at levels k-1 and k
};

// 3D Morton (Z-order) key of non-negative cell coordinates (< 2^21): consecutive keys are spatially compact, so
// batches of consecutive positions are small in space (tight culling bounds).
inline uint64_t morton3(uint64_t x, uint64_t y, uint64_t z) {
    auto spread = [](uint64_t v) {
        v &= 0x1fffff;
        v = (v | v << 32) & 0x1f00000000ffffULL;
        v = (v | v << 16) & 0x1f0000ff0000ffULL;
        v = (v | v << 8) & 0x100f00f00f00f00fULL;
        v = (v | v << 4) & 0x10c30c30c30c30c3ULL;
        v = (v | v << 2) & 0x1249249249249249ULL;
        return v;
    };
    return spread(x) | spread(y) << 1 | spread(z) << 2;
}

struct VecHash {
    size_t operator()(const std::vector<int>& k) const {
        uint64_t h = 1469598103934665603ULL;
        for (int v : k) h = (h ^ (uint64_t)(uint32_t)v) * 1099511628211ULL;
        return (size_t)h;
    }
};

inline Template3 build_template(const Params& p, int e0) {
    const int T = p.tau, EM = p.emax();
    const double dth = p.dth();
    const int per = (int)std::lround(p.nlon * p.dtheta_deg / 360.0);
    // DAG key: current (a, e), then the visited heading cells with their step counts and speed-change sums
    // (the position is the sum over steps of (v0 + ds dv) dt u(a, e), independent of the order).
    using Key = std::vector<int>;
    auto cell = [&](int a, int e) { return (a + T) * (2 * T + 1) + (e - e0 + T); };
    auto pos_of = [&](const Key& k) {
        std::array<double, 3> q{0, 0, 0};
        for (size_t s = 2; s < k.size(); s += 3) {
            const int c = k[s], a = c / (2 * T + 1) - T, e = c % (2 * T + 1) - T + e0;
            const double d = (k[s + 1] * p.v0 + k[s + 2] * p.dv) * p.dt;
            const double ca = std::cos(a * dth), sa = std::sin(a * dth), ce = std::cos(e * dth), se = std::sin(e * dth);
            q[0] += d * ce * ca; q[1] += d * ce * sa; q[2] += d * se;
        }
        return q;
    };
    struct Raw { int level, a, e; double x, y, z; };
    std::vector<Raw> raw{{0, 0, e0, 0, 0, 0}};
    std::vector<std::array<int, 3>> redges;
    std::vector<Key> cur{Key{0, e0}};
    std::vector<int> cur_ids{0};
    Template3 t;
    t.tau = T; t.e0 = e0;
    size_t paths = 1;
    for (int k = 1; k <= T; k++) {
        std::unordered_map<Key, int, VecHash> idx;
        std::vector<Key> nk;
        std::vector<int> ni;
        for (size_t q = 0; q < cur.size(); q++)
            for (int m = 0; m < NM; m++) {
                const int a = cur[q][0] + MOVES[m].da, e = cur[q][1] + MOVES[m].de;
                if (std::abs(e) > EM) continue;  // past the elevation cap: unavailable
                Key c = cur[q];
                c[0] = a; c[1] = e;
                const int cc = cell(a, e);
                size_t s = 2;
                while (s < c.size() && c[s] < cc) s += 3;
                if (s < c.size() && c[s] == cc) { c[s + 1] += 1; c[s + 2] += MOVES[m].ds; }
                else c.insert(c.begin() + s, {cc, 1, MOVES[m].ds});
                auto it = idx.find(c);
                int cid;
                if (it == idx.end()) {
                    const auto xyz = pos_of(c);
                    cid = (int)raw.size();
                    raw.push_back({k, a, e, xyz[0], xyz[1], xyz[2]});
                    idx.emplace(c, cid);
                    nk.push_back(std::move(c));
                    ni.push_back(cid);
                } else cid = it->second;
                redges.push_back({cur_ids[q], cid, m});
            }
        cur = std::move(nk);
        cur_ids = std::move(ni);
        paths *= NM;
        t.tree_nodes += paths;
    }
    // positions (merged on a grid), units, tiles, per level
    const double eps = p.merge > 0 ? p.merge : 1e-6;
    struct PosRec { std::array<long long, 3> tk; uint64_t tm = 0, pm = 0; double x, y, z; std::vector<int> nodes; };
    std::vector<int> new_id(raw.size(), -1), node_tile(1, -1);
    new_id[0] = 0;
    t.node_x = {0.f}; t.node_y = {0.f}; t.node_z = {0.f}; t.node_unit = {-1};
    t.ltile = {0, 0}; t.lpos = {0, 0};
    t.pos_u0 = {0}; t.pos_n0 = {1};
    std::vector<std::vector<int>> by_level(T + 1);
    for (int n = 1; n < (int)raw.size(); n++) by_level[raw[n].level].push_back(n);
    int ntile = 0;
    for (int k = 1; k <= T; k++) {
        std::map<std::array<long long, 3>, int> pidx;
        std::vector<PosRec> prs;
        for (int n : by_level[k]) {
            const std::array<long long, 3> key{std::llround(raw[n].x / eps), std::llround(raw[n].y / eps), std::llround(raw[n].z / eps)};
            auto it = pidx.find(key);
            int q;
            if (it == pidx.end()) { q = (int)prs.size(); pidx.emplace(key, q); prs.push_back({{0, 0, 0}, 0, 0, 0, {}}); }
            else q = it->second;
            prs[q].nodes.push_back(n);
        }
        for (auto& pr : prs) {
            for (int n : pr.nodes) { pr.x += raw[n].x; pr.y += raw[n].y; pr.z += raw[n].z; }
            pr.x /= pr.nodes.size(); pr.y /= pr.nodes.size(); pr.z /= pr.nodes.size();
            pr.tk = {(long long)std::floor(pr.x / p.tile), (long long)std::floor(pr.y / p.tile), (long long)std::floor(pr.z / p.tile)};
            std::stable_sort(pr.nodes.begin(), pr.nodes.end(), [&](int a, int b) { return raw[a].a < raw[b].a; });
        }
        // tiles, and positions within a tile, in Morton order (offsets make the coordinates non-negative)
        const long long off = (long long)std::ceil((T + 1) * (p.v0 + p.dv) * p.dt / std::min(p.tile, 1.0)) + 2;
        for (auto& pr : prs) {
            pr.tm = morton3(pr.tk[0] + off, pr.tk[1] + off, pr.tk[2] + off);
            pr.pm = morton3((uint64_t)std::llround(pr.x * 4) + 4 * off, (uint64_t)std::llround(pr.y * 4) + 4 * off, (uint64_t)std::llround(pr.z * 4) + 4 * off);
        }
        std::stable_sort(prs.begin(), prs.end(), [](const PosRec& a, const PosRec& b) { return a.tm != b.tm ? a.tm < b.tm : a.pm < b.pm; });
        for (size_t q = 0; q < prs.size(); q++) {
            if (q == 0 || prs[q].tm != prs[q - 1].tm) ntile++;
            const auto& nodes = prs[q].nodes;
            t.pos_x.push_back((float)prs[q].x); t.pos_y.push_back((float)prs[q].y); t.pos_z.push_back((float)prs[q].z);
            for (size_t z = 0; z < nodes.size(); z++) {
                const int n = nodes[z];
                if (z == 0 || raw[n].a != raw[nodes[z - 1]].a) {
                    const long long sh = ((long long)raw[n].a * per) % p.nlon;
                    t.unit_sh.push_back((int)((sh + p.nlon) % p.nlon));
                }
                new_id[n] = (int)t.node_x.size();
                t.node_x.push_back((float)raw[n].x); t.node_y.push_back((float)raw[n].y); t.node_z.push_back((float)raw[n].z);
                t.node_unit.push_back((int)t.unit_sh.size() - 1);
                node_tile.push_back(ntile - 1);
            }
            t.pos_u0.push_back((int)t.unit_sh.size());
            t.pos_n0.push_back((int)t.node_x.size());
        }
        t.ltile.push_back(ntile);
        t.lpos.push_back((int)t.pos_x.size());
    }
    t.M = (int)t.node_x.size(); t.P = (int)t.pos_x.size(); t.U = (int)t.unit_sh.size(); t.NT = ntile;
    // edges sorted by child tile; per-tile and per-level collision bboxes
    std::vector<std::array<int, 3>> es;
    for (auto& e : redges) es.push_back({new_id[e[0]], new_id[e[1]], e[2]});
    std::stable_sort(es.begin(), es.end(), [&](const auto& a, const auto& b) { return node_tile[a[1]] < node_tile[b[1]]; });
    t.tile_estart.assign(t.NT + 1, 0);
    t.tile_bb.assign(6 * (size_t)t.NT, 0.f);
    for (int q = 0; q < t.NT; q++)
        for (int d = 0; d < 3; d++) { t.tile_bb[6 * q + 2 * d] = 1e30f; t.tile_bb[6 * q + 2 * d + 1] = -1e30f; }
    for (auto& e : es) {
        t.eparent.push_back(e[0]); t.echild.push_back(e[1]); t.emove.push_back(e[2]);
        const int q = node_tile[e[1]];
        t.tile_estart[q + 1]++;
        for (int n : {e[0], e[1]}) {
            const float c[3] = {t.node_x[n], t.node_y[n], t.node_z[n]};
            for (int d = 0; d < 3; d++) {
                t.tile_bb[6 * q + 2 * d] = std::min(t.tile_bb[6 * q + 2 * d], c[d]);
                t.tile_bb[6 * q + 2 * d + 1] = std::max(t.tile_bb[6 * q + 2 * d + 1], c[d]);
            }
        }
    }
    for (int q = 0; q < t.NT; q++) t.tile_estart[q + 1] += t.tile_estart[q];
    t.lvl_bb.assign(6 * (size_t)(T + 1), 0.f);
    for (int k = 1; k <= T; k++) {
        float* bb = &t.lvl_bb[6 * k];
        for (int d = 0; d < 3; d++) { bb[2 * d] = 1e30f; bb[2 * d + 1] = -1e30f; }
        for (int q = t.ltile[k]; q < t.ltile[k + 1]; q++)
            for (int d = 0; d < 3; d++) {
                bb[2 * d] = std::min(bb[2 * d], t.tile_bb[6 * q + 2 * d]);
                bb[2 * d + 1] = std::max(bb[2 * d + 1], t.tile_bb[6 * q + 2 * d + 1]);
            }
    }
    return t;
}

// The templates of every elevation index (built once, in parallel; shared by all GPUs).
inline std::vector<Template3> build_templates(const Params& p) {
    const int EM = p.emax();
    std::vector<Template3> tpls(2 * EM + 1);
    std::vector<std::thread> th;
    for (int e = -EM; e <= EM; e++) th.emplace_back([&, e] { tpls[e + EM] = build_template(p, e); });
    for (auto& t : th) t.join();
    return tpls;
}

}  // namespace fsm3d
