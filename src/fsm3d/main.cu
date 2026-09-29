// fsm3d_gpu: the 3D FSM simulator on NVIDIA GPUs (docs/NOTES.md, "3D"; python/fsmtools/reference3d.py is the
// brute-force reference). Output files as the 2D simulators, with "dims": 3 and frames of (x, y, z, theta, phi, v).
//
//   bin/fsm3d_gpu --N 1000 --tau 6 --box 110 --grid 192x96 --merge 1.5 --steps 2000 --out runs/n1000_t6
//   bin/fsm3d_gpu ... --predator com --pred-start 100 --pred-grow 4 --pred2 intercept,0,70 --pred-vision on
//   bin/fsm3d_gpu --score STATE.txt --tau 4 --dump s.f32      # score one state (x y z theta phi v per line)
//   bin/fsm3d_gpu --benchmark 1000 --tau 5 --grid 192x96        # time one step on a synthetic flock
//
// Model flags: --N --tau --steps --v0 --dv --dtheta --R --dt --phi-max --half-angle --box --heading-sd --seed --init
// --out --log-every --resume; --grid NLONxNLAT (visual field; 192x96 is validated, nlon a multiple of 32 and of
// 360 / dtheta), --merge EPS (merge future positions within EPS: 1.5 recommended for tau >= 6), --tile.
// Predators: --predator line|com|nearest|intercept, --pred-start, --pred-dist, --pred-from AZ,EL, --pred-speed,
// --pred-radius, --pred-turn, --pred-grow, --pred-vision on|off, --pred-blind, --pred-memory, --pred-range,
// --pred2 MODE,AZ,EL, --pred-state x,y,z,vx,vy,vz[,...]. GPU: --device, --gpus 0,1,..|all, --batch, --band-config.
#include <array>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "../common/run_files.hpp"
#include "predator.hpp"
#include "scorer.cuh"

namespace fsm3d {

std::string params_json(const Params& p, const MultiScorer& s, int N, int steps_done, double seconds) {
    fsm::Json j;
    p.describe(j, N);
    j.str("backend", "cuda").num("dag_nodes", s.dag_nodes).num("tree_nodes", s.tree_nodes).num("positions", s.positions)
        .num("gpus", s.ngpu()).num("steps", steps_done).num("wall_seconds", seconds);
    fsm::Json files;
    files.str("traj", "float64 [steps+1, N, 6] (x, y, z, theta, phi, v)").str("moves", "int8 [steps, N]")
        .str("scores", "float32 [steps, N, 7] (-1: every future of that move dies)");
    j.obj("files", files);
    return j.text();
}

int simulate(Params& p, const GpuOptions& o) {
    State X = initial_state(p);
    const int N = (int)(X.size() / 6);
    p.N = N;
    MultiScorer sc(p, o, N);
    Predators pred(p, N);
    std::vector<fsm::RunFiles::Stream> extra;
    if (pred.count()) extra.push_back({".pred", pred.record_bytes(), true});
    if (pred.count() && p.pred.vision) extra.push_back({".seen", (size_t)N, false});
    const size_t fb = sizeof(double) * 6 * N, sb = sizeof(float) * NM * N;
    fsm::RunFiles files(p.out, fb, N, sb, extra);
    std::vector<char> frame;
    double prev_seconds = 0;
    const int start = files.open(params_json(p, sc, N, 0, 0), p.resume, frame, prev_seconds);
    if (start > 0) {
        std::memcpy(X.data(), frame.data(), fb);
        if (pred.count()) pred.restore(p.out, start);
    } else {
        files.write(0, X.data(), fb);
        if (pred.count()) {
            pred.appear(0, X);
            files.write(3, pred.record().data(), pred.record_bytes());
        }
    }
    std::vector<float> s;
    std::vector<int8_t> mv(N);
    State X0;
    const auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
    double gpu_ms = 0;
    for (int step = start; step < p.steps; step++) {
        const double ts = elapsed();
        if (pred.count()) {
            pred.before_step(step, X);
            sc.pred = pred.state;
            sc.know = pred.know;
            if (p.pred.vision) files.write(4, pred.flags.data(), N);
        }
        if (step == start || step % 250 == 0) sc.tune(X, s); else sc.score(X, s);
        gpu_ms += sc.last_ms;
        X0 = X;
        double ox = 0, oy = 0, oz = 0;
        for (int i = 0; i < N; i++) {
            auto rng = fsm::make_rng(p.seed, step, i);
            const int m = choose(&s[(size_t)i * NM], rng);
            mv[i] = (int8_t)m;
            const auto u = apply_move(p, X, i, m);
            ox += u[0]; oy += u[1]; oz += u[2];
        }
        const int caught = pred.count() ? pred.after_step(X0, X) : 0;
        files.write(1, mv.data(), N);
        files.write(2, s.data(), sb);
        files.write(0, X.data(), fb);
        if (pred.count()) {
            pred.appear(step + 1, X);
            files.write(3, pred.record().data(), pred.record_bytes());
        }
        if (p.log_every > 0 && (step % p.log_every == 0 || step == p.steps - 1))
            std::fprintf(stderr, "step %d  %.3fs (scoring %.3fs, GPU imbalance %.2f)  order %.4f%s\n", step, elapsed() - ts, sc.last_ms / 1e3,
                         sc.imbalance(), std::sqrt(ox * ox + oy * oy + oz * oz) / N, pred.summary(X, caught).c_str());
        if (step % 50 == 0 || step == p.steps - 1) files.write_json(params_json(p, sc, N, step + 1, prev_seconds + elapsed()));
    }
    const double secs = elapsed();
    const int ran = std::max(1, p.steps - start);
    std::fprintf(stderr, "done: %d steps (%d this session) in %.2fs (%.3fs/step, scoring %.3fs/step)\n", std::max(start, p.steps),
                 std::max(0, p.steps - start), secs, secs / ran, gpu_ms / 1e3 / ran);
    return 0;
}

// A jittered cubic lattice at the given spacing, headings near +x.
State synthetic_flock(int n, double spacing, double v0) {
    std::mt19937_64 g(7);
    std::normal_distribution<double> Gd(0, 1);
    const int side = (int)std::ceil(std::cbrt((double)n)) + 1;
    std::vector<std::array<double, 3>> pts;
    for (int a = 0; a < side; a++)
        for (int b = 0; b < side; b++)
            for (int c = 0; c < side; c++)
                pts.push_back({a * spacing + 0.25 * spacing * Gd(g), b * spacing + 0.25 * spacing * Gd(g), c * spacing + 0.25 * spacing * Gd(g)});
    State X;
    for (int j = 0; j < n; j++) X.insert(X.end(), {pts[j][0], pts[j][1], pts[j][2], 0.15 * Gd(g), 0.0, v0});
    return X;
}

// Score one state: print the kernel time, optionally dump the scores (float32 [N][7]) and, with a predator and
// vision, who sees it (DUMP.seen). Without memory: only the agents that see the predator now know of it.
int score_state(Params& p, const GpuOptions& o, State X, int reps, const std::string& dump) {
    const int N = (int)(X.size() / 6);
    p.N = N;
    MultiScorer sc(p, o, N);
    sc.pred = p.pred.state;
    if (!sc.pred.empty() && p.pred.vision) {
        std::vector<uint8_t> seen, all;
        sc.know.clear();
        for (size_t r = 0; r < sc.pred.size() / 6; r++) {
            predator_seen(p, X, N, &sc.pred[6 * r], seen);
            std::vector<int> last(N);
            for (int i = 0; i < N; i++) last[i] = seen[i] ? 0 : -1;
            std::vector<double> kr;
            predator_knowledge(p, N, 0, last, [&](int) { return &sc.pred[6 * r]; }, kr);
            sc.know.insert(sc.know.end(), kr.begin(), kr.end());
            all.insert(all.end(), seen.begin(), seen.end());
        }
        if (!dump.empty()) {
            std::ofstream fo(dump + ".seen", std::ios::binary);
            fo.write((const char*)all.data(), (std::streamsize)all.size());
        }
    }
    std::vector<float> s;
    sc.tune(X, s);  // also the warm-up
    double tot = 0, best = 1e30;
    if (reps <= 0) { tot = best = sc.last_ms; reps = 0; }
    for (int r = 0; r < reps; r++) { sc.score(X, s); tot += sc.last_ms; best = std::min(best, (double)sc.last_ms); }
    std::printf("N=%d tau=%d grid %dx%d merge=%g | kernels %.2f ms mean, %.2f ms best\n", N, p.tau, p.nlon, p.nlat, p.merge,
                tot / std::max(1, reps), best);
    if (!dump.empty()) {
        std::ofstream fo(dump, std::ios::binary);
        fo.write((const char*)s.data(), (std::streamsize)(s.size() * sizeof(float)));
    }
    return 0;
}

}  // namespace fsm3d

int main(int argc, char** argv) {
    using namespace fsm3d;
    try {
        Params p;
        GpuOptions o;
        std::string score, dump;
        int benchmark = 0, reps = 3;
        double spacing = 11;
        fsm::Args a(argc, argv);
        while (!a.done()) {
            const std::string f = a.next();
            if (p.parse(f, a)) continue;
            if (f == "--device") o.device = a.integer();
            else if (f == "--gpus") o.gpus = fsm::parse_gpus(a.str());
            else if (f == "--batch") o.batch = a.integer();
            else if (f == "--band-config") o.band_config = a.integer();
            else if (f == "--quiet") o.quiet = true;
            else if (f == "--score") score = a.str();
            else if (f == "--dump") dump = a.str();
            else if (f == "--benchmark") benchmark = a.integer();
            else if (f == "--spacing") spacing = a.real();
            else if (f == "--reps") reps = a.integer();
            else a.unknown(f);
        }
        p.finish();
        CK(cudaSetDevice(o.device));
        if (!score.empty()) return score_state(p, o, read_state(score), reps, dump);
        if (benchmark > 0) return score_state(p, o, synthetic_flock(benchmark, spacing, p.v0), reps, dump);
        return simulate(p, o);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
