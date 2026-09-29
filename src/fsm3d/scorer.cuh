// Host side of the 3D scoring. A Scorer owns one GPU's buffers and scores a range of agents; MultiScorer builds the
// templates once, splits every step's agents over the GPUs (contiguous ranges of equal predicted cost: each agent's
// band-kernel time in the previous step), and turns the exact row sums into scores in a fixed order, so scores (and
// trajectories) do not depend on the number of GPUs, the split or the launch configuration.
#pragma once

#include <chrono>
#include <exception>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "band_kernel.cuh"
#include "prune_kernel.cuh"
#include "template.hpp"

namespace fsm3d {

// Options of the GPU scorer (flags of fsm3d_gpu).
struct GpuOptions {
    int device = 0;
    std::vector<int> gpus;  // several GPUs: split every step's agents over them
    int batch = 0;          // agents per kernel pass (0: as many as fit in memory)
    int band_config = -1;   // band kernel configuration (index into band_choices; -1: autotuned)
    bool quiet = false;
};

struct BandChoice { void (*fn)(Args, int); int rb, pb, threads, list; size_t smem; };
template <int NLON, int RB, int PB, int TH, int L>
inline BandChoice band_choice() { return {band_kernel<NLON, RB, PB, TH, L>, RB, PB, TH, L, BandCfg<NLON, RB, PB, TH, L>::bytes}; }
inline std::vector<BandChoice> band_choices(int nlon) {  // candidate configurations, autotuned at start-up
    if (nlon == 384) return {band_choice<384, 4, 32, 512, 2048>(), band_choice<384, 4, 32, 1024, 4096>(), band_choice<384, 2, 32, 512, 2048>()};
    if (nlon == 768) return {band_choice<768, 2, 32, 512, 2048>(), band_choice<768, 2, 16, 512, 2048>()};
    if (nlon == 192) return {band_choice<192, 8, 32, 512, 2048>(), band_choice<192, 8, 32, 1024, 4096>(), band_choice<192, 8, 64, 1024, 4096>(),
                             band_choice<192, 4, 32, 512, 2048>(), band_choice<192, 8, 16, 512, 2048>()};
    if (nlon == 96) return {band_choice<96, 16, 32, 512, 2048>(), band_choice<96, 16, 32, 1024, 4096>(), band_choice<96, 8, 32, 512, 2048>()};
    throw std::runtime_error("--grid: nlon 96, 192, 384 or 768");
}

// One GPU: scores the agents of a range [lo, hi) into the shared host buffers (exact row sums, n_b, and each
// agent's cost in band-kernel block cycles, which balances the next step's split between GPUs).
struct Scorer {
    Params p;
    GpuOptions o;
    int N, EM, nt, dev;
    const std::vector<Template3>& tpls;
    std::vector<DevTemplate> dts;
    DevTemplate* d_tpl = nullptr;
    std::vector<void*> allocs;
    Args ka{};
    int grid_prune = 0, batch = 0;
    struct BandCfgRun { BandChoice c; int grid = 0; float ms = -1; };
    std::vector<BandCfgRun> bands;  // candidate configurations; cur: the fastest
    int cur = 0;
    float *dQX, *dQY, *dQZ, *dTH, *dPX, *dPY, *dPZ, *dPA;
    int* dtix;
    cudaEvent_t e0, e1;
    float last_ms = 0;

    template <class T> T* dalloc(size_t n) { void* ptr; CK(cudaMalloc(&ptr, std::max<size_t>(1, n) * sizeof(T))); allocs.push_back(ptr); return (T*)ptr; }
    template <class T> const T* upload(const std::vector<T>& v) { T* d = dalloc<T>(v.size()); CK(cudaMemcpy(d, v.data(), v.size() * sizeof(T), cudaMemcpyHostToDevice)); return d; }

    Scorer(const Params& pp, const GpuOptions& oo, int N_, int dev_, const std::vector<Template3>& tp) : p(pp), o(oo), N(N_), dev(dev_), tpls(tp) {
        CK(cudaSetDevice(dev));
        EM = p.emax(); nt = 2 * EM + 1;
        int Umax = 0, Mmax = 0;
        for (auto& t : tpls) { Umax = std::max(Umax, t.U); Mmax = std::max(Mmax, t.M); }
        for (auto& t : tpls) {
            DevTemplate d{};
            d.M = t.M; d.P = t.P; d.U = t.U;
            d.ltile = upload(t.ltile); d.lpos = upload(t.lpos); d.node_unit = upload(t.node_unit); d.unit_sh = upload(t.unit_sh);
            d.pos_u0 = upload(t.pos_u0); d.pos_n0 = upload(t.pos_n0); d.tile_estart = upload(t.tile_estart);
            d.eparent = upload(t.eparent); d.echild = upload(t.echild); d.emove = upload(t.emove);
            d.node_x = upload(t.node_x); d.node_y = upload(t.node_y); d.node_z = upload(t.node_z);
            d.pos_x = upload(t.pos_x); d.pos_y = upload(t.pos_y); d.pos_z = upload(t.pos_z);
            d.tile_bb = upload(t.tile_bb); d.lvl_bb = upload(t.lvl_bb);
            dts.push_back(d);
        }
        d_tpl = dalloc<DevTemplate>(dts.size());
        CK(cudaMemcpy(d_tpl, dts.data(), dts.size() * sizeof(DevTemplate), cudaMemcpyHostToDevice));
        cudaDeviceProp prop;
        CK(cudaGetDeviceProperties(&prop, dev));
        int per_sm = 0, grid_band = 0;
        auto choices = band_choices(p.nlon);
        if (o.band_config >= 0) {
            if (o.band_config >= (int)choices.size()) throw std::runtime_error("--band-config: 0.." + std::to_string(choices.size() - 1));
            choices = {choices[o.band_config]};
        }
        for (const auto& c : choices) {
            cudaFuncAttributes fa;
            CK(cudaFuncGetAttributes(&fa, c.fn));
            if (c.smem + fa.sharedSizeBytes > prop.sharedMemPerBlockOptin) continue;  // dynamic + static
            CK(cudaFuncSetAttribute(c.fn, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)c.smem));
            CK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&per_sm, c.fn, c.threads, c.smem));
            if (per_sm < 1) continue;
            bands.push_back({c, per_sm * prop.multiProcessorCount});
            grid_band = std::max(grid_band, bands.back().grid);
        }
        if (bands.empty()) throw std::runtime_error("no band kernel configuration fits");
        CK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&per_sm, prune_kernel<256>, 256, 0));
        grid_prune = std::max(1, per_sm) * prop.multiProcessorCount;
        // agents per pass: bounded by the unit-count buffers (and --batch)
        size_t freeb, totb;
        CK(cudaMemGetInfo(&freeb, &totb));
        const size_t per_agent = (size_t)Umax * NM * 4 + NM * 8 + (size_t)NM * p.nlat * 8 + 8;
        const size_t fixed = (size_t)grid_prune * ((size_t)Mmax * NM * 4 + 2 * (size_t)N * 16 * 2) + (size_t)grid_band * (N * 16 * 3 + OVF * 4 + PAIRS * 2);
        const size_t budget = freeb > fixed + (512ull << 20) ? (freeb - fixed - (512ull << 20)) / 2 : 0;
        batch = (int)std::min<size_t>(N, std::max<size_t>(1, budget / per_agent));
        if (o.batch > 0) batch = std::min(batch, o.batch);
        ka.tpl = d_tpl; ka.tau = p.tau; ka.N = N; ka.R = (float)p.R; ka.D2 = (float)(4 * p.R * p.R); ka.half_asin = p.half_asin;
        ka.Umax = Umax; ka.Mmax = Mmax;
        ka.mult = dalloc<int>((size_t)grid_prune * Mmax * NM);
        ka.cl = dalloc<float4>((size_t)grid_prune * 2 * N);
        ka.ct = dalloc<float4>((size_t)grid_prune * 2 * N);
        ka.um = dalloc<int>((size_t)batch * Umax * NM);
        ka.nb = dalloc<long long>((size_t)batch * NM);
        ka.cost = dalloc<unsigned long long>(batch);
        ka.cand = dalloc<float4>((size_t)grid_band * 3 * N);
        ka.ovf = dalloc<unsigned>((size_t)grid_band * OVF);
        ka.pairs = dalloc<unsigned short>((size_t)grid_band * PAIRS);
        ka.part = dalloc<long long>((size_t)batch * NM * p.nlat);
        ka.next_task = dalloc<int>(1);
        const size_t nq = (size_t)(p.tau + 1) * N;
        dQX = dalloc<float>(nq); dQY = dalloc<float>(nq); dQZ = dalloc<float>(nq);
        dTH = dalloc<float>(N); dtix = dalloc<int>(N);
        ka.QX = dQX; ka.QY = dQY; ka.QZ = dQZ; ka.TH = dTH; ka.tix = dtix;
        const int np = std::max(1, p.pred.count());
        dPX = dalloc<float>(np * nq); dPY = dalloc<float>(np * nq); dPZ = dalloc<float>(np * nq); dPA = dalloc<float>((size_t)np * N);
        ka.PX = dPX; ka.PY = dPY; ka.PZ = dPZ; ka.PA = dPA; ka.npred = 0;
        ka.pm = (float)(p.R + p.pred.radius); ka.pgrow = (float)p.pred.grow;
        CK(cudaEventCreate(&e0)); CK(cudaEventCreate(&e1));
        if (!o.quiet)
            std::fprintf(stderr, "GPU %d (%s): %d agents per pass, %zu band configurations, prune blocks %d\n", dev, prop.name, batch, bands.size(), grid_prune);
    }
    ~Scorer() { cudaSetDevice(dev); for (void* ptr : allocs) cudaFree(ptr); }

    // Agents [lo, hi): futures from the host arrays, results into part / nb / cost (full-size host buffers).
    void run(const float* hQX, const float* hQY, const float* hQZ, const float* hTH, const int* htix, const float* hP, int npred,
             int lo, int hi, long long* part, long long* nb, unsigned long long* cost) {
        CK(cudaSetDevice(dev));
        const size_t nq = (size_t)(p.tau + 1) * N;
        ka.npred = npred;
        if (npred) {  // hP: x, y, z [npred][(tau + 1) N] (each agent's prediction), then ages [npred][N]
            const size_t m = (size_t)npred * nq;
            CK(cudaMemcpy(dPX, hP, m * 4, cudaMemcpyHostToDevice));
            CK(cudaMemcpy(dPY, hP + m, m * 4, cudaMemcpyHostToDevice));
            CK(cudaMemcpy(dPZ, hP + 2 * m, m * 4, cudaMemcpyHostToDevice));
            CK(cudaMemcpy(dPA, hP + 3 * m, (size_t)npred * N * 4, cudaMemcpyHostToDevice));
        }
        CK(cudaMemcpy(dQX, hQX, nq * 4, cudaMemcpyHostToDevice));
        CK(cudaMemcpy(dQY, hQY, nq * 4, cudaMemcpyHostToDevice));
        CK(cudaMemcpy(dQZ, hQZ, nq * 4, cudaMemcpyHostToDevice));
        CK(cudaMemcpy(dTH, hTH, N * 4, cudaMemcpyHostToDevice));
        CK(cudaMemcpy(dtix, htix, N * 4, cudaMemcpyHostToDevice));
        const auto& bc = bands[cur];
        ka.nbands = (p.nlat + bc.c.rb - 1) / bc.c.rb;
        CK(cudaEventRecord(e0));
        for (int a0 = lo; a0 < hi; a0 += batch) {
            ka.a0 = a0; ka.na = std::min(batch, hi - a0);
            CK(cudaMemset(ka.next_task, 0, 4));
            CK(cudaMemset(ka.cost, 0, (size_t)ka.na * 8));
            prune_kernel<256><<<grid_prune, 256>>>(ka);
            CK(cudaGetLastError());
            CK(cudaMemset(ka.next_task, 0, 4));
            bc.c.fn<<<bc.grid, bc.c.threads, bc.c.smem>>>(ka, p.nlat);
            CK(cudaGetLastError());
            CK(cudaMemcpy(part + (size_t)a0 * NM * p.nlat, ka.part, (size_t)ka.na * NM * p.nlat * 8, cudaMemcpyDeviceToHost));
            CK(cudaMemcpy(nb + (size_t)a0 * NM, ka.nb, (size_t)ka.na * NM * 8, cudaMemcpyDeviceToHost));
            CK(cudaMemcpy(cost + a0, ka.cost, (size_t)ka.na * 8, cudaMemcpyDeviceToHost));
        }
        CK(cudaEventRecord(e1));
        CK(cudaEventSynchronize(e1));
        CK(cudaEventElapsedTime(&last_ms, e0, e1));
    }
};

// All GPUs: every step the agents are split into contiguous ranges of equal predicted cost (the previous step's
// per-agent cost), one per GPU, scored in parallel (a host thread per GPU). The scores do not depend on the split.
struct MultiScorer {
    Params p;
    GpuOptions o;
    int N, EM;
    std::vector<Template3> tpls;
    std::vector<std::unique_ptr<Scorer>> gs;
    std::vector<double> omega;      // solid angle of a pixel in each row
    std::vector<float> hQX, hQY, hQZ, hTH, hP;
    std::vector<int> htix;
    std::vector<double> pred;       // predators' states [np][6] (x, y, z, vx, vy, vz), world frame; empty: none
    std::vector<double> know;       // [np][N][7] each agent's knowledge of them: state projected to now, age (-1:
                                    // none); empty: every agent knows the current states (age 0)
    std::vector<long long> part, nb;
    std::vector<unsigned long long> cost;
    std::vector<int> split;         // [G + 1] agent ranges
    std::vector<float> gpu_ms;      // per GPU, last step
    float last_ms = 0;              // wall time of the scoring (all GPUs)
    size_t tree_nodes = 0, positions = 0, dag_nodes = 0;

    MultiScorer(const Params& pp, const GpuOptions& oo, int N_) : p(pp), o(oo), N(N_) {
        EM = p.emax();
        const auto tb = std::chrono::steady_clock::now();
        tpls = build_templates(p);
        const auto& t0 = tpls[EM];
        tree_nodes = t0.tree_nodes; positions = t0.P; dag_nodes = t0.M - 1;
        if (!o.quiet)
            std::fprintf(stderr, "3D templates: %d elevations, DAG %zu nodes (tree %zu), %zu positions (level elevation); grid %dx%d; built in %.1fs\n",
                         2 * EM + 1, dag_nodes, tree_nodes, positions, p.nlon, p.nlat,
                         std::chrono::duration<double>(std::chrono::steady_clock::now() - tb).count());
        std::vector<int> devs = o.gpus.empty() ? std::vector<int>{o.device} : o.gpus;
        gs.resize(devs.size());
        {
            std::vector<std::thread> th;
            for (size_t g = 0; g < devs.size(); g++) th.emplace_back([&, g] { gs[g] = std::make_unique<Scorer>(p, o, N, devs[g], tpls); });
            for (auto& t : th) t.join();
        }
        for (int j = 0; j < p.nlat; j++) {
            const double lo = -0.5 * M_PI + j * (M_PI / p.nlat);
            omega.push_back((TWO_PI / p.nlon) * (std::sin(lo + M_PI / p.nlat) - std::sin(lo)));
        }
        hQX.resize((size_t)(p.tau + 1) * N); hQY = hQX; hQZ = hQX; hTH.resize(N); htix.resize(N);
        part.resize((size_t)N * NM * p.nlat); nb.resize((size_t)N * NM); cost.assign(N, 1);
        gpu_ms.assign(gs.size(), 0.f);
    }

    int ngpu() const { return (int)gs.size(); }

    void prepare(const std::vector<double>& X) {
        const int T = p.tau;
        double cx = 0, cy = 0, cz = 0;  // re-centre (float32 on the device)
        for (int j = 0; j < N; j++) { cx += X[6 * j]; cy += X[6 * j + 1]; cz += X[6 * j + 2]; }
        cx /= N; cy /= N; cz /= N;
        const double dth = p.dtheta_deg * M_PI / 180.0;
        for (int j = 0; j < N; j++) {
            const double th = X[6 * j + 3], ph = X[6 * j + 4];
            const double ux = std::cos(ph) * std::cos(th) * p.v0 * p.dt, uy = std::cos(ph) * std::sin(th) * p.v0 * p.dt, uz = std::sin(ph) * p.v0 * p.dt;
            for (int k = 0; k <= T; k++) {
                hQX[(size_t)k * N + j] = (float)(X[6 * j] - cx + k * ux);
                hQY[(size_t)k * N + j] = (float)(X[6 * j + 1] - cy + k * uy);
                hQZ[(size_t)k * N + j] = (float)(X[6 * j + 2] - cz + k * uz);
            }
            hTH[j] = (float)std::remainder(th, TWO_PI);
            const long e = std::lround(ph / dth);
            if (std::fabs(ph - e * dth) > 1e-6 || std::abs(e) > EM) throw std::runtime_error("elevation off the lattice (agent " + std::to_string(j) + ")");
            htix[j] = (int)e + EM;
        }
        const size_t nq = (size_t)(T + 1) * N;
        const int np = (int)(pred.size() / 6);
        const size_t m = (size_t)np * nq;
        hP.assign(3 * m + (size_t)np * N, 0.f);
        for (int r = 0; r < np; r++)  // predicted straight on at the (known) velocity
            for (int j = 0; j < N; j++) {
                const double* q = know.empty() ? &pred[6 * (size_t)r] : &know[7 * ((size_t)r * N + j)];
                float& age = hP[3 * m + (size_t)r * N + j];
                age = know.empty() ? 0.f : (float)q[6];
                if (age < 0) continue;
                for (int k = 0; k <= T; k++) {
                    const size_t o = (size_t)r * nq + (size_t)k * N + j;
                    hP[o] = (float)(q[0] - cx + k * p.dt * q[3]);
                    hP[m + o] = (float)(q[1] - cy + k * p.dt * q[4]);
                    hP[2 * m + o] = (float)(q[2] - cz + k * p.dt * q[5]);
                }
            }
    }

    // contiguous ranges of equal cost (the last measured costs)
    void make_split() {
        const int G = ngpu();
        split.assign(G + 1, N);
        split[0] = 0;
        double tot = 0;
        for (int j = 0; j < N; j++) tot += (double)cost[j];
        double acc = 0;
        int g = 1;
        for (int j = 0; j < N && g < G; j++) {
            acc += (double)cost[j];
            while (g < G && acc >= tot * g / G) split[g++] = j + 1;
        }
        for (int q = 1; q <= G; q++) split[q] = std::max(split[q], split[q - 1]);
    }

    // scores of agents [lo, hi) from the row sums (fixed order: the same for every configuration and split)
    void reduce(int lo, int hi, std::vector<float>& out) const {
        for (int j = lo; j < hi; j++)
            for (int b = 0; b < NM; b++) {
                const long long n = nb[(size_t)j * NM + b];
                out[(size_t)j * NM + b] = -1.f;
                if (!n) continue;
                double s = 0;
                for (int q = 0; q < p.nlat; q++) s += omega[q] * (double)part[((size_t)j * NM + b) * p.nlat + q];
                out[(size_t)j * NM + b] = (float)(s / (TWO_PI * (double)n));
            }
    }

    void run_all(bool tune, std::vector<float>& out) {
        const auto t0 = std::chrono::steady_clock::now();
        make_split();
        std::vector<std::thread> th;
        std::vector<std::exception_ptr> err(gs.size());
        for (int g = 0; g < ngpu(); g++)
            th.emplace_back([&, g] {
                try {
                    Scorer& s = *gs[g];
                    const int lo = split[g], hi = split[g + 1];
                    if (lo >= hi) { gpu_ms[g] = 0; return; }
                    if (tune && s.bands.size() > 1) {
                        for (size_t c = 0; c < s.bands.size(); c++) {
                            s.cur = (int)c;
                            s.run(hQX.data(), hQY.data(), hQZ.data(), hTH.data(), htix.data(), hP.data(), (int)(pred.size() / 6), lo, hi, part.data(), nb.data(), cost.data());
                            s.bands[c].ms = s.last_ms;
                        }
                        s.cur = 0;
                        for (size_t c = 1; c < s.bands.size(); c++) if (s.bands[c].ms < s.bands[s.cur].ms) s.cur = (int)c;
                    }
                    s.run(hQX.data(), hQY.data(), hQZ.data(), hTH.data(), htix.data(), hP.data(), (int)(pred.size() / 6), lo, hi, part.data(), nb.data(), cost.data());
                    gpu_ms[g] = s.last_ms;
                    reduce(lo, hi, out);
                } catch (...) { err[g] = std::current_exception(); }
            });
        for (auto& t : th) t.join();
        for (auto& e : err) if (e) std::rethrow_exception(e);
        last_ms = (float)(1e3 * std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        if (tune && !o.quiet)
            for (int g = 0; g < ngpu(); g++) {
                const Scorer& s = *gs[g];
                std::string all;
                for (auto& b : s.bands) all += " " + std::to_string(b.c.rb) + "r" + std::to_string(b.c.pb) + "p" + std::to_string(b.c.threads) + "t" + std::to_string(b.c.list / 1024) + "k:" + std::to_string((int)std::lround(b.ms)) + "ms";
                std::fprintf(stderr, "autotune GPU %d: bands of %d rows, batches of %d, %d threads, lists of %d |%s\n", s.dev, s.bands[s.cur].c.rb,
                             s.bands[s.cur].c.pb, s.bands[s.cur].c.threads, s.bands[s.cur].c.list, all.c_str());
            }
    }

    // X: [N][6] (x, y, z, theta, phi, v), world frame. Scores out [N][NM] (-1: every future of that move dies).
    void score(const std::vector<double>& X, std::vector<float>& out, bool tune = false) {
        prepare(X);
        out.assign((size_t)N * NM, -1.f);
        run_all(tune, out);
    }
    void tune(const std::vector<double>& X, std::vector<float>& out) { score(X, out, true); }
    float imbalance() const {  // slowest GPU / mean
        double mx = 0, mean = 0;
        for (float v : gpu_ms) { mx = std::max(mx, (double)v); mean += v; }
        mean /= gpu_ms.size();
        return mean > 0 ? (float)(mx / mean) : 1.f;
    }
};

}  // namespace fsm3d
