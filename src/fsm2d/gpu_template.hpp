// The DAG (dag.hpp) regrouped for the GPU kernel (docs/OPTIMIZATIONS.md, "GPU 2D"). Agent independent; offsets in
// the agent's frame.
//   node      a DAG node (level, heading, exact position): path counts live here
//   position  a projection point: nodes at the same spot (or within --merge of each other) share one visual mask
//   unit      the nodes of a position with the same heading: identical for scoring, so their path counts are summed
//             and accumulated once
//   group     <= G positions / <= GU units, processed together (painted as one batch)
//   tile      spatial cell of positions at one level; collision candidates are filtered per tile
// Nodes are numbered level, tile, position, heading-major, so every tile/group/position/unit owns a contiguous node
// range. Node 0 is the root.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <vector>

#include "dag.hpp"

namespace fsm2d {

constexpr int G = 32;   // positions per group
constexpr int GU = 64;  // units per group (a position has at most 2 tau + 1)

struct GpuTemplate {
    int tau = 0, M = 0, P = 0, U = 0, NT = 0, NG = 0;
    size_t tree_nodes = 0;
    std::vector<int> lstart, ltile;                 // [tau+2]: node / tile ranges per level
    std::vector<float> node_x, node_y;              // exact node offsets (collision geometry)
    std::vector<int> node_unit;
    std::vector<int> unit_sh, unit_pos;             // unit: heading as a bin shift in [0, B), its position
    std::vector<float> pos_x, pos_y;                // projection point of each position
    std::vector<int> tile_gstart, tile_estart;      // groups / incoming edges per tile
    std::vector<float> tile_bb;                     // bbox of each tile's nodes and their parents (4 per tile)
    std::vector<int> grp_p0, grp_u0, grp_u1, grp_n0, grp_n1;
    std::vector<int> eparent, echild, emove;        // edges sorted by child tile
    std::vector<float> lvl_bb;                      // [4 (tau+1)]: bbox of the nodes at levels k-1..k
    bool exact_shift = true;                        // dtheta is a whole number of bins
};

// B: angular bins of the visual field; merge: position grid (0: exact); tile: tile size.
inline GpuTemplate build_gpu_template(const Params& p, int B, double merge, double tile_size) {
    const Dag d = build_dag(p);
    const int T = p.tau;
    GpuTemplate t;
    t.tau = T;
    t.tree_nodes = d.tree_nodes;

    // Positions (nodes merged on a grid), units, tiles and groups, per level.
    const double eps = merge > 0 ? merge : 1e-6;
    const double tile = tile_size > 0 ? tile_size : 1e30;
    struct PosRec { long long tx, ty; double x, y; std::vector<int> nodes; };
    std::vector<int> new_id(d.size(), -1), node_tile(1, -1);
    t.lstart = {0, 1};
    t.ltile = {0, 0};
    t.node_x = {0.f}; t.node_y = {0.f}; t.node_unit = {-1};
    new_id[0] = 0;
    const double bins_per_dth = p.dth() / (TWO_PI / B);
    t.exact_shift = std::fabs(bins_per_dth - std::round(bins_per_dth)) < 1e-9;
    for (int k = 1; k <= T; k++) {
        std::map<std::pair<long long, long long>, int> pidx;
        std::vector<PosRec> prs;
        for (int n = d.lstart[k]; n < d.lstart[k + 1]; n++) {
            const auto key = std::make_pair(std::llround(d.x[n] / eps), std::llround(d.y[n] / eps));
            auto it = pidx.find(key);
            int q;
            if (it == pidx.end()) { q = (int)prs.size(); pidx.emplace(key, q); prs.push_back({0, 0, 0, 0, {}}); }
            else q = it->second;
            prs[q].nodes.push_back(n);
        }
        for (auto& pr : prs) {  // projection point = mean of the merged nodes; nodes ordered by heading
            for (int n : pr.nodes) { pr.x += d.x[n]; pr.y += d.y[n]; }
            pr.x /= pr.nodes.size(); pr.y /= pr.nodes.size();
            pr.tx = (long long)std::floor(pr.x / tile);
            pr.ty = (long long)std::floor(pr.y / tile);
            std::stable_sort(pr.nodes.begin(), pr.nodes.end(), [&](int a, int b) { return d.h[a] < d.h[b]; });
        }
        std::stable_sort(prs.begin(), prs.end(), [](const PosRec& a, const PosRec& b) { return a.tx != b.tx ? a.tx < b.tx : a.ty < b.ty; });
        for (size_t q = 0; q < prs.size();) {  // one tile per grid cell
            size_t q1 = q;
            while (q1 < prs.size() && prs[q1].tx == prs[q].tx && prs[q1].ty == prs[q].ty) q1++;
            const int tid = (int)t.tile_gstart.size();
            t.tile_gstart.push_back((int)t.grp_p0.size());
            int gp = 0, gu = 0;
            for (size_t r = q; r < q1; r++) {
                const auto& nodes = prs[r].nodes;
                int nu = 0;
                for (size_t z = 0; z < nodes.size(); z++) nu += (z == 0 || d.h[nodes[z]] != d.h[nodes[z - 1]]);
                if (gp == 0 || gp == G || gu + nu > GU) {  // new group
                    t.grp_p0.push_back((int)t.pos_x.size());
                    t.grp_u0.push_back((int)t.unit_sh.size());
                    t.grp_n0.push_back((int)t.node_x.size());
                    gp = 0; gu = 0;
                }
                const int pid = (int)t.pos_x.size();
                t.pos_x.push_back((float)prs[r].x); t.pos_y.push_back((float)prs[r].y);
                for (size_t z = 0; z < nodes.size(); z++) {
                    const int n = nodes[z];
                    if (z == 0 || d.h[n] != d.h[nodes[z - 1]]) {
                        const long long sh = std::llround(d.h[n] * bins_per_dth) % B;
                        t.unit_sh.push_back((int)((sh + B) % B));
                        t.unit_pos.push_back(pid);
                    }
                    new_id[n] = (int)t.node_x.size();
                    t.node_x.push_back((float)d.x[n]); t.node_y.push_back((float)d.y[n]);
                    t.node_unit.push_back((int)t.unit_sh.size() - 1);
                    node_tile.push_back(tid);
                }
                gp++; gu += nu;
                t.grp_u1.resize(t.grp_p0.size()); t.grp_n1.resize(t.grp_p0.size());
                t.grp_u1.back() = (int)t.unit_sh.size();
                t.grp_n1.back() = (int)t.node_x.size();
            }
            q = q1;
        }
        t.lstart.push_back((int)t.node_x.size());
        t.ltile.push_back((int)t.tile_gstart.size());
    }
    t.M = (int)t.node_x.size(); t.P = (int)t.pos_x.size(); t.U = (int)t.unit_sh.size();
    t.NT = (int)t.tile_gstart.size(); t.NG = (int)t.grp_p0.size();
    t.tile_gstart.push_back(t.NG);

    // Edges sorted by child tile; per-tile collision bbox (children and their parents); per-level bbox.
    std::vector<std::array<int, 3>> es;
    for (size_t e = 0; e < d.eparent.size(); e++) es.push_back({new_id[d.eparent[e]], new_id[d.echild[e]], d.emove[e]});
    std::stable_sort(es.begin(), es.end(), [&](const auto& a, const auto& b) { return node_tile[a[1]] < node_tile[b[1]]; });
    t.tile_estart.assign(t.NT + 1, 0);
    t.tile_bb.assign(4 * (size_t)t.NT, 0.f);
    for (int q = 0; q < t.NT; q++) { t.tile_bb[4 * q] = 1e30f; t.tile_bb[4 * q + 1] = -1e30f; t.tile_bb[4 * q + 2] = 1e30f; t.tile_bb[4 * q + 3] = -1e30f; }
    for (auto& e : es) {
        t.eparent.push_back(e[0]); t.echild.push_back(e[1]); t.emove.push_back(e[2]);
        const int q = node_tile[e[1]];
        t.tile_estart[q + 1]++;
        for (int n : {e[0], e[1]}) {
            float* bb = &t.tile_bb[4 * q];
            bb[0] = std::min(bb[0], t.node_x[n]); bb[1] = std::max(bb[1], t.node_x[n]);
            bb[2] = std::min(bb[2], t.node_y[n]); bb[3] = std::max(bb[3], t.node_y[n]);
        }
    }
    for (int q = 0; q < t.NT; q++) t.tile_estart[q + 1] += t.tile_estart[q];
    t.lvl_bb.assign(4 * (size_t)(T + 1), 0.f);
    for (int k = 1; k <= T; k++) {
        float* bb = &t.lvl_bb[4 * k];
        bb[0] = 1e30f; bb[1] = -1e30f; bb[2] = 1e30f; bb[3] = -1e30f;
        for (int q = t.ltile[k]; q < t.ltile[k + 1]; q++) {
            bb[0] = std::min(bb[0], t.tile_bb[4 * q]); bb[1] = std::max(bb[1], t.tile_bb[4 * q + 1]);
            bb[2] = std::min(bb[2], t.tile_bb[4 * q + 2]); bb[3] = std::max(bb[3], t.tile_bb[4 * q + 3]);
        }
    }
    return t;
}

}  // namespace fsm2d
