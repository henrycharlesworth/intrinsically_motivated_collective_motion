// Host side of the 2D GPU scoring: a GpuScorer owns one GPU's buffers and scores a range of agents from a snapshot;
// MultiScorer splits every step's agents over several GPUs. Scores are per agent and bit-identical whatever the
// launch configuration or the split, so trajectories do not depend on either.
#pragma once

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "gpu_kernel.cuh"

namespace fsm2d {

// Options of the GPU scorer (flags of fsm2d_gpu).
struct GpuOptions {
    int bins = 3072;       // angular bins of the visual field: 1536, 3072 or 6144 (a multiple of 360 / dtheta)
    double merge = -1;     // position merging grid; -1: automatic (exact for tau <= 7, 1.5 for tau >= 8)
    double tile = 24;      // tile size for the collision culling
    int threads = 0;       // threads per block (0: autotuned)
    int d_shared = -1;     // difference arrays in shared memory (1) or global (0); -1: autotuned
    int cand_cap = 4096, edge_cap = 4096;
    int retune = 250;      // re-run the autotuner every this many steps (0: only at the start)
    int device = 0;
    std::vector<int> gpus;  // several GPUs: split every step's agents over them
    bool quiet = false;

    // Positions within 1.5 of each other are projected from their mean at tau >= 8: 3-5x faster, with >= 99.8% of
    // decisions unchanged (docs/OPTIMIZATIONS.md).
    void resolve(int tau) {
        if (merge < 0) merge = tau >= 8 ? 1.5 : 0.0;
    }
};

using KernelFn = void (*)(KernelArgs);
struct KernelChoice { KernelFn fn; int threads; size_t mask_bytes; };

inline KernelChoice pick_kernel(int B, int threads) {
#define KC(BB, TT) \
    if (B == BB && threads == TT) return {score_kernel<BB, TT>, TT, SmemLayout<BB>::mask_bytes};
    KC(1536, 256) KC(3072, 256) KC(6144, 256) KC(1536, 512) KC(3072, 512) KC(6144, 512) KC(3072, 1024) KC(6144, 1024)
#undef KC
    return {nullptr, threads, 0};
}

class GpuScorer {
public:
    // One launch configuration: kernel (threads), where the difference arrays live, shared memory, grid.
    struct Config { KernelChoice kc; int d_shared = 0; size_t smem = 0; int per_sm = 0, grid = 0; float ms = -1; bool competitive = true; };
    GpuTemplate t;
    float last_ms = 0;
    int a0 = 0, a1 = 0;  // agents scored by score() (default all)

    GpuScorer(const Params& p, const GpuOptions& o, int N, int dev) : p_(p), o_(o), N_(N), dev_(dev) {
        CK(cudaSetDevice(dev_));
        a1 = N;
        t = build_gpu_template(p, o.bins, o.merge, o.tile);
        CK(cudaGetDeviceProperties(&prop_, dev_));
        // Even, so the work area stays 16-byte aligned. When the other agents don't all fit they are read from global
        // memory, and nothing is reserved (an unused buffer would only cost occupancy).
        ka_.cand_cap = N - 1 <= o.cand_cap ? N / 2 * 2 : 0;
        ka_.edge_cap = o.edge_cap;
        for (int th : {256, 512, 1024}) {  // candidate launch configurations
            if (o.threads > 0 && th != o.threads) continue;
            const KernelChoice kc = pick_kernel(o.bins, th);
            if (!kc.fn) continue;
            for (int ds : {0, 1}) {
                if (o.d_shared >= 0 && ds != o.d_shared) continue;
                Config c;
                c.kc = kc;
                c.d_shared = ds;
                cudaFuncAttributes fa;
                CK(cudaFuncGetAttributes(&fa, c.kc.fn));
                const size_t smem_max = (size_t)prop_.sharedMemPerBlockOptin - fa.sharedSizeBytes;
                c.smem = (size_t)ka_.cand_cap * 8 + c.kc.mask_bytes + (size_t)ka_.edge_cap * 4 + (ds ? (size_t)NM * o.bins * 4 : 0);
                if (c.smem > smem_max) continue;
                CK(cudaFuncSetAttribute(c.kc.fn, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)c.smem));
                CK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&c.per_sm, c.kc.fn, th, c.smem));
                if (c.per_sm < 1) continue;
                c.grid = std::max(1, std::min(N, c.per_sm * prop_.multiProcessorCount));
                configs_.push_back(c);
            }
        }
        if (configs_.empty())
            throw std::runtime_error("no launch configuration (--bins 1536|3072|6144, --threads 256|512|1024, or lower --cand-cap / --edge-cap)");
        int max_grid = 0;
        for (auto& c : configs_) max_grid = std::max(max_grid, c.grid);
        DevTemplate& d = ka_.t;
        d.tau = t.tau; d.M = t.M;
        d.lstart = upload(t.lstart); d.ltile = upload(t.ltile); d.node_unit = upload(t.node_unit);
        d.unit_sh = upload(t.unit_sh); d.unit_pos = upload(t.unit_pos);
        d.node_x = upload(t.node_x); d.node_y = upload(t.node_y); d.pos_x = upload(t.pos_x); d.pos_y = upload(t.pos_y);
        d.tile_gstart = upload(t.tile_gstart); d.tile_estart = upload(t.tile_estart);
        d.tile_bb = upload(t.tile_bb); d.lvl_bb = upload(t.lvl_bb);
        d.grp_p0 = upload(t.grp_p0); d.grp_u0 = upload(t.grp_u0); d.grp_u1 = upload(t.grp_u1);
        d.grp_n0 = upload(t.grp_n0); d.grp_n1 = upload(t.grp_n1);
        d.eparent = upload(t.eparent); d.echild = upload(t.echild); d.emove = upload(t.emove);
        ka_.N = N;
        ka_.R = (float)p.R; ka_.D2 = (float)(4 * p.R * p.R);
        ka_.half_asin = p.half_asin ? 1 : 0;
        const size_t nq = (size_t)(p.tau + 1) * N, grid = (size_t)max_grid;  // per-block scratch for the largest grid
        dQX_ = alloc<float>(nq); dQY_ = alloc<float>(nq); dTH_ = alloc<float>(N);
        ka_.QX = dQX_; ka_.QY = dQY_; ka_.TH = dTH_;
        ka_.mult = alloc<int>(grid * t.M * NM);
        ka_.lvl = alloc<float2>(grid * N);
        ka_.cl = alloc<float4>(grid * N);
        ka_.ct = alloc<float4>(grid * N);
        ka_.diff = alloc<int>(grid * NM * o.bins);
        ka_.edges = alloc<int>(grid * G * o.bins);
        ka_.next_agent = alloc<int>(1);
        ka_.score = alloc<float>((size_t)N * NM);
        hQX_.resize(nq); hQY_.resize(nq); hTH_.resize(N);
        CK(cudaEventCreate(&e0_)); CK(cudaEventCreate(&e1_));
        if (!o.quiet)
            std::fprintf(stderr, "GPU %d %s (%d SMs): %d bins | DAG %d nodes (full tree %zu), %d positions, %d units, %d tiles | merge %g%s\n",
                         dev_, prop_.name, prop_.multiProcessorCount, o.bins, t.M - 1, t.tree_nodes, t.P, t.U, t.NT, o.merge,
                         t.exact_shift ? "" : " | WARNING: dtheta is not a whole number of bins (headings rounded)");
    }
    ~GpuScorer() {
        cudaSetDevice(dev_);
        for (void* ptr : allocs_) cudaFree(ptr);
    }
    GpuScorer(const GpuScorer&) = delete;
    GpuScorer& operator=(const GpuScorer&) = delete;

    // Scores of agents [a0, a1) into out ([N][NM]; other entries are left alone). Positions in double (world frame).
    void score(const State& s, std::vector<float>& out) {
        CK(cudaSetDevice(dev_));
        const int T = p_.tau, N = N_;
        double cx = 0, cy = 0;  // re-centre: the device works in float32
        for (int j = 0; j < N; j++) { cx += s.x[j]; cy += s.y[j]; }
        cx /= N; cy /= N;
        for (int j = 0; j < N; j++) {
            const double ux = std::cos(s.th[j]) * p_.v0 * p_.dt, uy = std::sin(s.th[j]) * p_.v0 * p_.dt;
            for (int k = 0; k <= T; k++) {
                hQX_[(size_t)k * N + j] = (float)(s.x[j] - cx + k * ux);
                hQY_[(size_t)k * N + j] = (float)(s.y[j] - cy + k * uy);
            }
            hTH_[j] = (float)std::remainder(s.th[j], TWO_PI);
        }
        CK(cudaMemcpy(dQX_, hQX_.data(), hQX_.size() * 4, cudaMemcpyHostToDevice));
        CK(cudaMemcpy(dQY_, hQY_.data(), hQY_.size() * 4, cudaMemcpyHostToDevice));
        CK(cudaMemcpy(dTH_, hTH_.data(), hTH_.size() * 4, cudaMemcpyHostToDevice));
        CK(cudaMemset(ka_.next_agent, 0, 4));
        const Config& c = configs_[cur_];
        ka_.d_shared = c.d_shared;
        ka_.a0 = a0; ka_.a1 = a1;
        CK(cudaEventRecord(e0_));
        c.kc.fn<<<std::max(1, std::min(c.grid, a1 - a0)), c.kc.threads, c.smem>>>(ka_);
        CK(cudaGetLastError());
        CK(cudaEventRecord(e1_));
        out.resize((size_t)N * NM);
        CK(cudaMemcpy(out.data() + (size_t)a0 * NM, ka_.score + (size_t)a0 * NM, (size_t)(a1 - a0) * NM * 4, cudaMemcpyDeviceToHost));
        CK(cudaEventElapsedTime(&last_ms, e0_, e1_));
    }

    // Time the launch configurations on this snapshot and keep the fastest; out gets its scores (identical for every
    // configuration). The first call times all of them (after a warm-up run each); later calls re-time only those
    // within 25% of the best, once each.
    void tune(const State& s, std::vector<float>& out) {
        if (configs_.size() == 1) {
            score(s, out);
            return;
        }
        const bool first = !tuned_;
        std::vector<int> timed;
        for (size_t c = 0; c < configs_.size(); c++) {
            if (!first && !configs_[c].competitive) continue;
            cur_ = (int)c;
            if (first) score(s, out);  // warm-up
            score(s, out);
            configs_[c].ms = last_ms;
            timed.push_back((int)c);
        }
        cur_ = timed[0];
        for (int c : timed) if (configs_[c].ms < configs_[cur_].ms) cur_ = c;
        if (first) {
            for (auto& c : configs_) c.competitive = c.ms <= 1.25f * configs_[cur_].ms;
            tuned_ = true;
        }
        last_ms = configs_[cur_].ms;
        if (!o_.quiet) {
            std::string all;
            for (int c : timed)
                all += " " + std::to_string(configs_[c].kc.threads) + (configs_[c].d_shared ? "s" : "g") + ":" + std::to_string((int)std::lround(configs_[c].ms)) + "ms";
            const Config& c = configs_[cur_];
            std::fprintf(stderr, "autotune GPU %d: %d threads, D in %s memory, %zu KB shared, %d blocks/SM |%s\n", dev_, c.kc.threads,
                         c.d_shared ? "shared" : "global", c.smem / 1024, c.per_sm, all.c_str());
        }
    }

private:
    template <class E> E* alloc(size_t n) {
        void* ptr;
        CK(cudaMalloc(&ptr, std::max<size_t>(n * sizeof(E), 16)));
        allocs_.push_back(ptr);
        return (E*)ptr;
    }
    template <class E> const E* upload(const std::vector<E>& v) {
        E* d = alloc<E>(std::max<size_t>(1, v.size()));
        if (!v.empty()) CK(cudaMemcpy(d, v.data(), v.size() * sizeof(E), cudaMemcpyHostToDevice));
        return d;
    }

    Params p_;
    GpuOptions o_;
    int N_, dev_;
    cudaDeviceProp prop_{};
    std::vector<Config> configs_;
    int cur_ = 0;
    bool tuned_ = false;
    KernelArgs ka_{};
    std::vector<void*> allocs_;
    std::vector<float> hQX_, hQY_, hTH_;
    float *dQX_ = nullptr, *dQY_ = nullptr, *dTH_ = nullptr;
    cudaEvent_t e0_, e1_;
};

// Several GPUs: every GPU holds all positions and scores a contiguous share of the agents, one host thread each. The
// shares follow each GPU's measured speed (smoothed), so the split stays balanced as the flock changes.
class MultiScorer {
public:
    float last_ms = 0, imbalance = 1;

    MultiScorer(const Params& p, const GpuOptions& o, int N) : N_(N) {
        const std::vector<int> devs = o.gpus.empty() ? std::vector<int>{o.device} : o.gpus;
        for (int d : devs) g_.push_back(std::make_unique<GpuScorer>(p, o, N, d));
        rate_.assign(g_.size(), 1.0);
        split();
    }
    const GpuTemplate& t() const { return g_[0]->t; }
    int gpus() const { return (int)g_.size(); }

    void score(const State& s, std::vector<float>& out) { run([&](GpuScorer& g) { g.score(s, out); }, out); }
    void tune(const State& s, std::vector<float>& out) { run([&](GpuScorer& g) { g.tune(s, out); }, out); }

private:
    void split() {  // shares proportional to rate
        double tot = 0;
        for (double r : rate_) tot += r;
        double acc = 0;
        for (size_t k = 0; k < g_.size(); k++) {
            g_[k]->a0 = (int)std::lround(N_ * acc / tot);
            acc += rate_[k];
            g_[k]->a1 = k + 1 == g_.size() ? N_ : (int)std::lround(N_ * acc / tot);
        }
    }
    template <class F> void run(F f, std::vector<float>& out) {
        out.resize((size_t)N_ * NM);
        if (g_.size() == 1) {
            f(*g_[0]);
            last_ms = g_[0]->last_ms;
            return;
        }
        std::vector<std::thread> th;
        std::vector<std::string> err(g_.size());
        for (size_t k = 0; k < g_.size(); k++)
            th.emplace_back([&, k] { try { f(*g_[k]); } catch (const std::exception& e) { err[k] = e.what(); } });
        for (auto& x : th) x.join();
        for (auto& e : err) if (!e.empty()) throw std::runtime_error(e);
        float mx = 0, sum = 0;
        for (size_t k = 0; k < g_.size(); k++) {
            const float ms = std::max(1e-3f, g_[k]->last_ms);
            mx = std::max(mx, ms);
            sum += ms;
            rate_[k] = 0.5 * rate_[k] + 0.5 * (g_[k]->a1 - g_[k]->a0) / ms;
        }
        last_ms = mx;
        imbalance = mx / (sum / g_.size());
        split();
    }

    int N_;
    std::vector<std::unique_ptr<GpuScorer>> g_;
    std::vector<double> rate_;  // agents per ms, per GPU
};

}  // namespace fsm2d
