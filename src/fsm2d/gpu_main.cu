// fsm2d_gpu: the 2D FSM simulator on NVIDIA GPUs. Same model, random streams and output files as the CPU simulator
// (bin/fsm2d); scores are computed on the GPU (float32 geometry, the visual field sampled on --bins angular bins,
// optionally merged future positions), decisions and motion on the host in double precision. The GPU chooses the
// same move as the exact CPU code for >= 99.8% of agents; the rest are near-ties (docs/OPTIMIZATIONS.md).
//
//   bin/fsm2d_gpu --N 500 --tau 6 --box 250 --steps 1000 --out runs/n500_t6
//   bin/fsm2d_gpu --N 50000 --tau 6 --box 2500 --steps 10000 --gpus all --out runs/n50000_t6    # split over GPUs
//   bin/fsm2d_gpu --score runs/n500_t6 --frames 200,600 --dump s.f32     # re-score saved frames (for validation)
//   bin/fsm2d_gpu --benchmark 2000 --tau 7                               # time one step on a synthetic flock
//
// Model flags as bin/fsm2d. GPU flags: --bins 1536|3072|6144, --merge EPS (default: 0 for tau <= 7, 1.5 above),
// --device D, --gpus 0,1,..|all, and the tuning knobs --threads, --d-shared, --tile, --cand-cap, --edge-cap, --retune.
#include <chrono>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "gpu_scorer.cuh"
#include "simulate.hpp"

namespace fsm2d {

class GpuBackend {
public:
    using Score = float;

    GpuBackend(const Params& p, const GpuOptions& o, int N) : o_(o), ms_(p, o, N) {}

    void score(const State& s, int step, std::vector<Score>& out) {
        if (first_ || (o_.retune > 0 && step % o_.retune == 0)) ms_.tune(s, out);
        else ms_.score(s, out);
        first_ = false;
        busy_ms_ += ms_.last_ms;
    }
    void describe(fsm::Json& j) const {
        const GpuTemplate& t = ms_.t();
        j.str("backend", "cuda").num("bins", o_.bins).num("merge", o_.merge).num("tile", o_.tile)
            .num("dag_nodes", t.M - 1).num("tree_nodes", t.tree_nodes).num("positions", t.P);
    }
    std::string log(int) const {
        char buf[96];
        if (ms_.gpus() > 1) std::snprintf(buf, sizeof buf, "  (GPU %.3fs, imbalance %.2f)", ms_.last_ms / 1e3, ms_.imbalance);
        else std::snprintf(buf, sizeof buf, "  (GPU %.3fs)", ms_.last_ms / 1e3);
        return buf;
    }
    double busy_seconds() const { return busy_ms_ / 1e3; }

private:
    GpuOptions o_;
    MultiScorer ms_;
    bool first_ = true;
    double busy_ms_ = 0;
};

// --score PREFIX: frames of a saved run (its model parameters are read from PREFIX.json).
// The run's tau is used unless --tau is given (to score its frames with another horizon).
std::vector<State> load_frames(Params& p, const std::string& prefix, const std::string& spec, bool keep_tau, std::vector<int>& frames) {
    const std::string js = fsm::read_text(prefix + ".json");
    p.N = (int)fsm::json_number(js, "N", p.N);
    if (!keep_tau) p.tau = (int)fsm::json_number(js, "tau", p.tau);
    p.v0 = fsm::json_number(js, "v0", p.v0); p.dv = fsm::json_number(js, "dv", p.dv);
    p.dtheta_deg = fsm::json_number(js, "dtheta_deg", p.dtheta_deg);
    p.R = fsm::json_number(js, "R", p.R); p.dt = fsm::json_number(js, "dt", p.dt);
    p.half_asin = js.find("\"half_angle\": \"asin\"") != std::string::npos;
    std::ifstream f(prefix + ".traj", std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + prefix + ".traj");
    f.seekg(0, std::ios::end);
    const long nframes = (long)(f.tellg() / (sizeof(double) * 4 * p.N));
    for (const auto& item : fsm::Args::split(spec.empty() ? std::to_string(nframes - 2) : spec, ',')) {
        const auto parts = fsm::Args::split(item, ':');  // FRAME or FROM:TO[:STEP]
        if (parts.size() == 1) frames.push_back(std::stoi(parts[0]));
        else for (int fr = std::stoi(parts[0]); fr < std::stoi(parts[1]); fr += parts.size() > 2 ? std::stoi(parts[2]) : 1) frames.push_back(fr);
    }
    std::vector<State> out;
    std::vector<char> buf(sizeof(double) * 4 * p.N);
    for (int fr : frames) {
        if (fr < 0 || fr >= nframes) throw std::runtime_error("frame " + std::to_string(fr) + " out of range");
        f.seekg((std::streamoff)fr * buf.size());
        f.read(buf.data(), (std::streamsize)buf.size());
        State s;
        s.x.resize(p.N); s.y.resize(p.N); s.th.resize(p.N); s.v.resize(p.N);
        s.load_frame(buf);
        out.push_back(s);
    }
    return out;
}

// --benchmark N: an aligned synthetic flock (jittered lattice in a disc, nearest-neighbour spacing ~ `spacing`).
State synthetic_flock(int N, double spacing, double v0) {
    std::mt19937_64 g(1);
    std::normal_distribution<double> Gd(0, 1);
    const int side = (int)std::ceil(std::sqrt(N * 4 / M_PI)) + 2;
    std::vector<std::array<double, 2>> pts;
    for (int a = -side; a <= side; a++)
        for (int b = -side; b <= side; b++) pts.push_back({a * spacing + 0.25 * spacing * Gd(g), b * spacing + 0.25 * spacing * Gd(g)});
    std::sort(pts.begin(), pts.end(), [](auto& A, auto& B) { return A[0] * A[0] + A[1] * A[1] < B[0] * B[0] + B[1] * B[1]; });
    State s;
    for (int j = 0; j < N; j++) { s.x.push_back(pts[j][0]); s.y.push_back(pts[j][1]); s.th.push_back(0.15 * Gd(g)); s.v.push_back(v0); }
    return s;
}

// Score the given snapshots on one GPU; print the kernel time and optionally dump the scores (float32 [F][N][5]).
int score_frames(Params& p, GpuOptions& o, std::vector<State> states, const std::vector<int>& frames, int reps, const std::string& dump) {
    p.N = states[0].size();
    o.resolve(p.tau);
    o.gpus.clear();
    GpuScorer gs(p, o, p.N, o.device);
    std::vector<float> sc, all;
    gs.tune(states[0], sc);
    double tot = 0, best = 1e30;
    int runs = 0;
    for (auto& s : states) {
        gs.score(s, sc);  // warm-up for this frame
        for (int r = 0; r < reps; r++) {
            gs.score(s, sc);
            tot += gs.last_ms; best = std::min(best, (double)gs.last_ms); runs++;
        }
        all.insert(all.end(), sc.begin(), sc.end());
    }
    std::printf("N=%d tau=%d bins=%d merge=%g | %zu frame(s) | kernel %.2f ms mean, %.2f ms best\n", p.N, p.tau, o.bins, o.merge,
                states.size(), runs ? tot / runs : 0.0, runs ? best : 0.0);
    if (!dump.empty()) {
        std::ofstream fo(dump, std::ios::binary);
        fo.write((const char*)all.data(), (std::streamsize)(all.size() * sizeof(float)));
        std::string fl;
        for (size_t q = 0; q < frames.size(); q++) fl += (q ? "," : "") + std::to_string(frames[q]);
        fsm::write_text(dump + ".json", "{\"N\": " + std::to_string(p.N) + ", \"n_moves\": 5, \"frames\": [" + fl + "]}\n");
    }
    return 0;
}

}  // namespace fsm2d

int main(int argc, char** argv) {
    using namespace fsm2d;
    try {
        Params p;
        GpuOptions o;
        std::string score_prefix, frames, dump;
        int benchmark = 0, reps = 3;
        double spacing = 7.8;
        bool tau_given = false;
        fsm::Args a(argc, argv);
        while (!a.done()) {
            const std::string f = a.next();
            tau_given |= f == "--tau";
            if (p.parse(f, a)) continue;
            if (f == "--bins") o.bins = a.integer();
            else if (f == "--merge") o.merge = a.real();
            else if (f == "--tile") o.tile = a.real();
            else if (f == "--threads") o.threads = a.integer();
            else if (f == "--d-shared") o.d_shared = a.integer();
            else if (f == "--cand-cap") o.cand_cap = a.integer();
            else if (f == "--edge-cap") o.edge_cap = a.integer();
            else if (f == "--retune") o.retune = a.integer();
            else if (f == "--device") o.device = a.integer();
            else if (f == "--gpus") o.gpus = fsm::parse_gpus(a.str());
            else if (f == "--quiet") o.quiet = true;
            else if (f == "--score") score_prefix = a.str();
            else if (f == "--frames") frames = a.str();
            else if (f == "--dump") dump = a.str();
            else if (f == "--benchmark") benchmark = a.integer();
            else if (f == "--spacing") spacing = a.real();
            else if (f == "--reps") reps = a.integer();
            else a.unknown(f);
        }
        p.finish();
        CK(cudaSetDevice(o.device));
        if (!score_prefix.empty()) {
            std::vector<int> fr;
            auto states = load_frames(p, score_prefix, frames, tau_given, fr);
            return score_frames(p, o, states, fr, reps, dump);
        }
        if (benchmark > 0) return score_frames(p, o, {synthetic_flock(benchmark, spacing, p.v0)}, {0}, reps, dump);
        o.resolve(p.tau);
        State s = initial_state(p);
        GpuBackend be(p, o, s.size());
        return simulate(p, s, be);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
