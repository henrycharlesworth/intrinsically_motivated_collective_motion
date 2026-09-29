// fsm2d: the 2D FSM simulator on the CPU (exact; OpenMP over agents). The reference implementation.
//
//   bin/fsm2d --N 500 --tau 6 --box 250 --steps 1000 --out runs/n500_t6
//   bin/fsm2d --N 500 --tau 6 --box 250 --steps 2000 --out runs/n500_t6 --resume   # continue or extend a run
//
// Model flags (both 2D simulators): --N --tau --steps --v0 --dv --dtheta --R --dt --half-angle atan|asin --noise
// --box --heading-sd --seed --init FILE --out PREFIX --log-every --resume. CPU only: --threads (default: all cores);
// --verify (step 0: check every score against the O(n^2) definition); --full-tree (use the full decision tree instead
// of the DAG; the results must be identical).
#include <cstdio>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "cpu_evaluator.hpp"
#include "simulate.hpp"

namespace fsm2d {

class CpuBackend {
public:
    using Score = double;

    CpuBackend(const Params& p, int threads, bool merge_dag, bool verify) : eval_(p, merge_dag), verify_(verify) {
#ifdef _OPENMP
        if (threads > 0) omp_set_num_threads(threads);
        nthreads_ = omp_get_max_threads();
#else
        (void)threads;
#endif
        ws_.resize(nthreads_);
        const Dag& d = eval_.dag();
        std::fprintf(stderr, "fsm2d (CPU, %d threads): tau %d, DAG %d nodes (full tree %zu)\n", nthreads_, p.tau, d.size() - 1, d.tree_nodes);
    }

    void score(const State& s, int step, std::vector<Score>& out) {
        const auto t0 = std::chrono::steady_clock::now();
        eval_.set_state(s);
        const int N = s.size();
        out.resize((size_t)NM * N);
        const bool check = verify_ && step == 0;
        double err = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(max : err)
        for (int i = 0; i < N; i++) {
#ifdef _OPENMP
            Workspace& ws = ws_[omp_get_thread_num()];
#else
            Workspace& ws = ws_[0];
#endif
            const Scores sc = eval_.evaluate(ws, i);
            for (int b = 0; b < NM; b++) out[(size_t)NM * i + b] = sc[b];
            if (check) {
                const Scores bf = eval_.brute_force(ws);
                for (int b = 0; b < NM; b++) err = std::max(err, std::fabs(bf[b] - sc[b]) / std::max(1e-12, std::fabs(bf[b])));
            }
        }
        if (check) std::fprintf(stderr, "verify: max relative |W - W_brute_force| = %.3e\n", err);
        busy_ += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }

    void describe(fsm::Json& j) const {
        const Dag& d = eval_.dag();
        j.str("backend", "cpu").num("dag_nodes", d.size() - 1).num("tree_nodes", d.tree_nodes);
    }

    std::string log(int steps) {
        uint64_t nodes = 0, edges = 0;
        for (auto& ws : ws_) { nodes += ws.nodes; edges += ws.edges; ws.nodes = ws.edges = 0; }
        const double per = (double)std::max(1, steps) * std::max<size_t>(1, n_agents_);
        char buf[96];
        std::snprintf(buf, sizeof buf, "  nodes/agent %.0f  edges/agent %.0f", nodes / per, edges / per);
        return buf;
    }
    double busy_seconds() const { return busy_; }
    void set_agents(size_t n) { n_agents_ = n; }

private:
    CpuEvaluator eval_;
    bool verify_;
    int nthreads_ = 1;
    std::vector<Workspace> ws_;
    double busy_ = 0;
    size_t n_agents_ = 1;
};

}  // namespace fsm2d

int main(int argc, char** argv) {
    using namespace fsm2d;
    try {
        Params p;
        int threads = 0;
        bool verify = false, merge_dag = true;
        fsm::Args a(argc, argv);
        while (!a.done()) {
            const std::string f = a.next();
            if (p.parse(f, a)) continue;
            if (f == "--threads") threads = a.integer();
            else if (f == "--verify") verify = true;
            else if (f == "--full-tree") merge_dag = false;
            else a.unknown(f);
        }
        p.finish();
        State s = initial_state(p);
        CpuBackend be(p, threads, merge_dag, verify);
        be.set_agents(s.size());
        return simulate(p, s, be);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
