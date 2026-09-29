// The 2D simulation loop, shared by the CPU and CUDA simulators: score every agent's moves from the current
// snapshot, let every agent choose, move everyone together, write the step.
//
// A Backend provides
//   using Score = ...;                                   score type (double on the CPU, float on the GPU)
//   void score(const State&, int step, std::vector<Score>& out)   out: [N][NM] (-1: every future collides)
//   void describe(fsm::Json&) const                      its fields of PREFIX.json
//   std::string log(int steps) const                     extra text for the progress line (and resets counters)
//   double busy_seconds() const                          time spent scoring so far (for the summary line)
#pragma once

#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "../common/run_files.hpp"
#include "model.hpp"

namespace fsm2d {

template <class Backend>
std::string params_json(const Params& p, const Backend& be, int n_agents, int steps_done, double seconds) {
    fsm::Json j;
    p.describe(j, n_agents);
    be.describe(j);
    j.num("steps", steps_done).num("wall_seconds", seconds);
    fsm::Json files;
    files.str("traj", "float64 [steps+1, N, 4] (x, y, theta, v)").str("moves", "int8 [steps, N]")
        .str("scores", "float32 [steps, N, 5] (-1: every future of that move collides)");
    j.obj("files", files);
    return j.text();
}

template <class Backend>
int simulate(Params& p, State& s, Backend& be) {
    const int N = s.size();
    fsm::RunFiles files(p.out, sizeof(double) * 4 * N, N, sizeof(float) * NM * N);
    std::vector<char> frame;
    double prev_seconds = 0;
    const int start = files.open(params_json(p, be, N, 0, 0), p.resume, frame, prev_seconds);
    if (start > 0) {
        s.load_frame(frame);
    } else {
        const auto f = s.frame();
        files.write(0, f.data(), f.size());
    }

    std::vector<typename Backend::Score> sc;
    std::vector<int8_t> moves(N);
    std::vector<float> scores(NM * (size_t)N);
    const auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
    int last_log = start - 1;
    for (int step = start; step < p.steps; step++) {
        const double ts = elapsed();
        be.score(s, step, sc);
        double ox = 0, oy = 0;
        for (int i = 0; i < N; i++) {
            auto rng = fsm::make_rng(p.seed, step, i);
            const int m = choose(&sc[(size_t)NM * i], rng);
            moves[i] = (int8_t)m;
            for (int b = 0; b < NM; b++) scores[(size_t)NM * i + b] = (float)sc[(size_t)NM * i + b];
            apply_move(p, s, i, m, step);
            ox += std::cos(s.th[i]);
            oy += std::sin(s.th[i]);
        }
        files.write(1, moves.data(), moves.size());
        files.write(2, scores.data(), scores.size() * sizeof(float));
        const auto f = s.frame();
        files.write(0, f.data(), f.size());
        if (p.log_every > 0 && (step % p.log_every == 0 || step == p.steps - 1)) {
            std::fprintf(stderr, "step %d  %.3fs  order %.4f%s\n", step, elapsed() - ts, std::hypot(ox, oy) / N,
                         be.log(step - last_log).c_str());
            last_log = step;
        }
        if (step % 50 == 0) files.write_json(params_json(p, be, N, step + 1, prev_seconds + elapsed()));
    }
    const double secs = elapsed();
    files.write_json(params_json(p, be, N, std::max(start, p.steps), prev_seconds + secs));
    const int ran = std::max(1, p.steps - start);
    std::fprintf(stderr, "done: %d steps (%d this session) in %.2fs (%.3fs/step, scoring %.3fs/step)\n", std::max(start, p.steps),
                 std::max(0, p.steps - start), secs, secs / ran, be.busy_seconds() / ran);
    return 0;
}

}  // namespace fsm2d
