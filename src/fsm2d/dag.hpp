// The future-state DAG of one agent (docs/OPTIMIZATIONS.md, "The future-state DAG").
//
// The naive decision tree has 5^k nodes at depth k. But a node's position is the sum of its steps, and that sum does
// not depend on their order: it is fixed by the current heading index h and, for each heading visited, the number of
// steps taken at it and their summed speed changes. Move sequences with the same such key reach the same point with
// the same heading, so their subtrees are identical. Merging them turns the tree into a DAG (4,594 nodes instead of
// 19,530 at tau = 6). Pruning collisions becomes counting, for every node and first move, the collision-free paths
// into it, which gives exactly the scores of the full tree.
//
// The DAG is agent independent: node offsets are in the agent's own frame (x along its current heading) and are
// rotated and translated per agent. Nodes are numbered level by level (node 0 is the root); edges are listed level by
// level in parent order.
#pragma once

#include <cmath>
#include <map>
#include <vector>

#include "model.hpp"

namespace fsm2d {

struct Dag {
    int tau = 0;
    std::vector<int> level, h;       // per node: depth and heading index in [-tau, tau]
    std::vector<double> x, y;        // per node: offset in the agent's frame
    std::vector<int> lstart;         // nodes of level k: [lstart[k], lstart[k+1])
    std::vector<int> eparent, echild, emove;
    std::vector<int> estart;         // edges into level k: [estart[k], estart[k+1]), k = 1..tau
    size_t tree_nodes = 0;           // nodes of the equivalent full tree (levels 1..tau)
    int size() const { return (int)h.size(); }
};

// merge = false keeps the full tree (every node a distinct path): only for checking that the DAG is exact.
inline Dag build_dag(const Params& p, bool merge = true) {
    const int T = p.tau, H = 2 * T + 1;
    const double dth = p.dth();
    Dag d;
    d.tau = T;
    using Key = std::vector<int>;  // [h, n_0, s_0, ..., n_{H-1}, s_{H-1}]: steps and speed-change sum per heading
    auto add = [&](const Key& k, int level) {
        double x = 0, y = 0;
        for (int slot = 0; slot < H; slot++) {
            const int n = k[1 + 2 * slot], s = k[2 + 2 * slot];
            if (n == 0) continue;
            const double dist = (n * p.v0 + s * p.dv) * p.dt, ang = (slot - T) * dth;
            x += dist * std::cos(ang);
            y += dist * std::sin(ang);
        }
        d.level.push_back(level);
        d.h.push_back(k[0]);
        d.x.push_back(x);
        d.y.push_back(y);
        return d.size() - 1;
    };
    std::vector<Key> cur{Key(1 + 2 * H, 0)};
    std::vector<int> cur_ids{add(cur[0], 0)};
    d.lstart = {0, 1};
    d.estart.assign(T + 2, 0);
    size_t paths = 1;
    for (int k = 1; k <= T; k++) {
        d.estart[k] = (int)d.eparent.size();
        std::map<Key, int> index;
        std::vector<Key> next;
        std::vector<int> next_ids;
        for (size_t q = 0; q < cur.size(); q++)
            for (int m = 0; m < NM; m++) {
                Key c = cur[q];
                const int hh = c[0] + MOVES[m].dh;
                c[0] = hh;
                c[1 + 2 * (hh + T)] += 1;
                c[2 + 2 * (hh + T)] += MOVES[m].ds;
                auto it = merge ? index.find(c) : index.end();
                int id;
                if (it == index.end()) {
                    id = add(c, k);
                    if (merge) index.emplace(c, id);
                    next.push_back(std::move(c));
                    next_ids.push_back(id);
                } else {
                    id = it->second;
                }
                d.eparent.push_back(cur_ids[q]);
                d.echild.push_back(id);
                d.emove.push_back(m);
            }
        d.lstart.push_back(d.size());
        cur = std::move(next);
        cur_ids = std::move(next_ids);
        paths *= NM;
        d.tree_nodes += paths;
    }
    d.estart[T + 1] = (int)d.eparent.size();
    return d;
}

}  // namespace fsm2d
