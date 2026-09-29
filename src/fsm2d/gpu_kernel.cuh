// The 2D score kernel (docs/OPTIMIZATIONS.md, "GPU 2D"). One thread block per agent (persistent blocks pulling
// agents from a queue); everything is in the agent's own frame, so the template's offsets are used directly and the
// other agents' futures are rotated once per level. Per tile of the template:
//   1. collision pruning of the edges into the tile (path counts per first move, as on the CPU);
//   2. per group of positions: paint one bitmask of B angular bins per reachable position (the other agents' discs),
//      list each mask's interval edges, and let every unit add +-(its path counts) at its edges, shifted by its
//      heading (a whole number of bins), into a difference array D of S_b.
// At the end S_b = the prefix sum of D, and W_b = (1 / (pi n_b)) sum_bins S_b (n_b - S_b) dtheta, in exact integer
// arithmetic until the final scaling, so scores do not depend on the launch configuration or the order of work.
#pragma once

#include <type_traits>

#include "../common/cuda_util.cuh"
#include "gpu_template.hpp"

namespace fsm2d {

struct DevTemplate {
    int tau, M;
    const int *lstart, *ltile, *node_unit, *unit_sh, *unit_pos;
    const float *node_x, *node_y, *pos_x, *pos_y;
    const int *tile_gstart, *tile_estart;
    const float *tile_bb, *lvl_bb;
    const int *grp_p0, *grp_u0, *grp_u1, *grp_n0, *grp_n1;
    const int *eparent, *echild, *emove;
};

struct KernelArgs {
    DevTemplate t;
    int N;
    float R, D2;
    int half_asin;
    const float *QX, *QY, *TH;  // futures of every agent, level-major [(tau+1) N]; headings
    int* mult;                  // [grid][M][NM] collision-free path counts per first move
    float2* lvl;                // [grid][N] other agents at the current level, agent frame (when not in shared memory)
    float4* cl;                 // [grid][N] collision candidates of the current level (positions at k-1 and k)
    float4* ct;                 // [grid][N] collision candidates of the current tile (overflow of shared memory)
    int* diff;                  // [grid][NM][B] difference arrays of S_b (when not in shared memory)
    int* edges;                 // [grid][G][B] edge lists that overflow shared memory
    int cand_cap, edge_cap;     // shared-memory capacities (float2 candidates, int edges)
    int d_shared;               // difference arrays in shared memory
    int* next_agent;
    int a0, a1;                 // this launch scores agents a0 .. a1 - 1 (a GPU's share)
    float* score;               // [N][NM]
};

// round to nearest without the (quarter-rate) float->int conversion; valid for |v| < 2^22
__device__ __forceinline__ int fast_rint(float v) { return __float_as_int(v + 12582912.f) - 0x4B400000; }

// Paint the angular interval subtended by a disc of radius R at offset (dx, dy) into a bitmask of B bins (bin q is set
// when its centre (q + 1/2) 2 pi / B lies inside the interval). ASIN: half-angle asin(R/d), else atan(R/d).
template <int B, bool ASIN>
__device__ __forceinline__ void paint(unsigned* bits, float dx, float dy, float R, float inv_bw) {
    const float d2 = dx * dx + dy * dy;
    if (d2 <= 0.f) return;
    const float x = R * rsqrtf(d2);
    float half;
    if (ASIN) {
        if (x >= 1.f) {  // inside another disc: everything occupied
            for (int w = 0; w < B / 32; w++) atomicOr(&bits[w], 0xffffffffu);
            return;
        }
        half = asinf(x);
    } else {
        half = fsm::half_atan(x);
    }
    const float ang = fsm::fast_atan2(dy, dx);
    // first/last bin whose centre is inside [ang - half, ang + half]: ceil(u - 1/2) = rint(u), floor(u - 1/2) =
    // rint(u - 1). Offset by B so both are non-negative (ang > -pi, half < pi/2).
    const float c = ang * inv_bw + (float)B, h = half * inv_bw;
    const int k0 = fast_rint(c - h), k1 = fast_rint(c + h - 1.f);
    if (k1 < k0) return;  // narrower than a bin and missing every centre
    // Branch-lean painting (lanes of a warp paint discs from different positions, so the word alignment varies
    // across lanes): the first and last words unconditionally, a loop only for the full words in between.
    constexpr int W = B / 32;
    const int w0 = k0 >> 5, w1 = k1 >> 5;
    const int wa = w0 >= W ? w0 - W : w0;  // k0 in (B/4, 3B/2): at most one wrap
    const unsigned mfirst = 0xffffffffu << (k0 & 31), mlast = 0xffffffffu >> (31 - (k1 & 31));
    atomicOr(&bits[wa], w1 == w0 ? (mfirst & mlast) : mfirst);
    if (w1 != w0) {
        int w = wa + 1 == W ? 0 : wa + 1;
        for (int q = w0 + 1; q < w1; q++) {
            atomicOr(&bits[w], 0xffffffffu);
            w = w + 1 == W ? 0 : w + 1;
        }
        atomicOr(&bits[w], mlast);
    }
}

// Dynamic shared memory: [ candidates: cand_cap float2 | D: NM B ints if d_shared | work area ], where the work area
// holds the tile's collision candidates (float4) or, during painting, G position masks followed by the group's edge
// lists (edge_cap ints).
template <int B>
struct SmemLayout {
    static constexpr int MW = B / 32 + 1;  // padded mask stride: lanes painting different masks hit different banks
    static constexpr size_t mask_bytes = (size_t)G * MW * 4;
};

// Units add +-(their path counts) at their mask's edges, rotated by their heading, into D; a unit whose mask covers
// its own bin 0 adds its counts to the base values instead (S_b = base_b + prefix sum of D).
template <int B, int THREADS>
__device__ __forceinline__ void scatter_units(int nu, const int* uslot, const int* slotmap, const int* ush, const int (*um)[NM],
                                              const unsigned* masks, const int* eoff, const int* elist, int d_shared,
                                              int* Ds, int* Dg, int* basep) {
    constexpr int MW = SmemLayout<B>::MW, NW = THREADS / 32;
    const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5;
    for (int c = warp; c < nu; c += NW) {  // one warp per unit
        const int slot = uslot[c];
        if (slot < 0) continue;
        const int g = slotmap[slot], sh = ush[c];
        int mb[NM];
#pragma unroll
        for (int b = 0; b < NM; b++) mb[b] = um[c][b];
        if (lane == 0 && ((masks[g * MW + (sh >> 5)] >> (sh & 31)) & 1u)) {
#pragma unroll
            for (int b = 0; b < NM; b++) basep[b] += mb[b];
        }
        auto scatter = [&](int* Dp) {
            for (int e = eoff[g] + lane; e < eoff[g + 1]; e += 32) {
                const int ev = elist[e];
                int bin = (ev >> 1) - sh;
                if (bin < 0) bin += B;
                if (bin == 0) continue;  // bin 0 is carried by the base count
                const bool end = ev & 1;
#pragma unroll
                for (int b = 0; b < NM; b++)
                    if (mb[b]) atomicAdd(&Dp[b * B + bin], end ? -mb[b] : mb[b]);
            }
        };
        if (d_shared) scatter(Ds); else scatter(Dg);
    }
}

// W_b from D, the base counts and n_b. Integers until the final scaling (deterministic).
template <int B, int THREADS>
__device__ __forceinline__ void finalize_scores(const int* D, const long long* nb, const int* base, long long (*redl)[THREADS / 32],
                                                int (*redi)[THREADS / 32], float* score) {
    constexpr int BPT = B / THREADS, NW = THREADS / 32;
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5;
#pragma unroll 1
    for (int b = 0; b < NM; b++) {
        // block exclusive scan of per-thread chunk sums (thread owns bins [tid BPT, (tid + 1) BPT))
        const int* Db = D + b * B + tid * BPT;
        int local = 0;
        for (int q = 0; q < BPT; q++) local += (tid == 0 && q == 0) ? 0 : Db[q];
        int incl = local;
        for (int o = 1; o < 32; o <<= 1) {
            const int v = __shfl_up_sync(0xffffffffu, incl, o);
            if (lane >= o) incl += v;
        }
        if (lane == 31) redi[0][warp] = incl;
        __syncthreads();
        int woff = 0;
        for (int w = 0; w < warp; w++) woff += redi[0][w];
        __syncthreads();
        long long S = base[b] + woff + incl - local, part = 0;
        for (int q = 0; q < BPT; q++) {
            if (!(tid == 0 && q == 0)) S += Db[q];
            part += S * (nb[b] - S);
        }
        for (int o = 16; o > 0; o >>= 1) part += __shfl_xor_sync(0xffffffffu, part, o);
        if (lane == 0) redl[b][warp] = part;
    }
    __syncthreads();
    if (tid < NM) {
        long long tot = 0;
        for (int w = 0; w < NW; w++) tot += redl[tid][w];
        score[tid] = nb[tid] > 0 ? (float)((double)tot * (2.0 / B) / (double)nb[tid]) : -1.f;
    }
}

template <int B, int THREADS>
__global__ void __launch_bounds__(THREADS) score_kernel(KernelArgs a) {
    constexpr int WORDS = B / 32, MW = SmemLayout<B>::MW, NW = THREADS / 32;
    static_assert(B % THREADS == 0 && B % 32 == 0, "bad bin/thread configuration");
    const DevTemplate& t = a.t;
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5, N = a.N, T = t.tau, nc = N - 1;
    int* mult = a.mult + (size_t)blockIdx.x * t.M * NM;
    int* edge_g = a.edges + (size_t)blockIdx.x * G * B;
    float2* lvlg = a.lvl + (size_t)blockIdx.x * N;
    float4* clv = a.cl + (size_t)blockIdx.x * N;
    float4* ctg = a.ct + (size_t)blockIdx.x * N;

    extern __shared__ __align__(16) unsigned char smem[];
    float2* s_cand = (float2*)smem;
    int* const Ds = (int*)(smem + (size_t)a.cand_cap * 8);
    int* const D = a.d_shared ? Ds : a.diff + (size_t)blockIdx.x * NM * B;
    unsigned char* work = smem + (size_t)a.cand_cap * 8 + (a.d_shared ? (size_t)NM * B * 4 : 0);
    unsigned* s_mask = (unsigned*)work;                      // [G][MW]
    int* s_edge = (int*)(work + SmemLayout<B>::mask_bytes);  // [edge_cap]
    float4* s_cc = (float4*)work;                            // [cc_cap] (collision phase only)
    const int cc_cap = (int)((SmemLayout<B>::mask_bytes + (size_t)a.edge_cap * 4) / 16);

    __shared__ float s_px[G], s_py[G];
    __shared__ int s_um[GU][NM];  // summed path counts of the group's units
    __shared__ int s_ush[GU], s_uslot[GU], s_greach[G], s_slotmap[G], s_ecnt[G], s_eoff[G + 1];
    __shared__ int s_n, s_ns, s_agent;
    __shared__ long long s_redl[NM][NW];
    __shared__ int s_redi[NM][NW];

    const float R = a.R, inv_bw = B / (2.f * (float)M_PI);
    for (;;) {
        __syncthreads();
        if (tid == 0) s_agent = atomicAdd(a.next_agent, 1);
        __syncthreads();
        const int i = a.a0 + s_agent;
        if (i >= a.a1) return;
        const float x0 = a.QX[i], y0 = a.QY[i], th0 = a.TH[i];
        float c0, s0;
        sincosf(th0, &s0, &c0);
        for (int q = tid; q < t.M * NM; q += THREADS) mult[q] = 0;
        for (int q = tid; q < NM * B; q += THREADS) D[q] = 0;
        int nbp[NM], basep[NM];  // per-thread partials: n_b, and S_b at bin 0
#pragma unroll
        for (int b = 0; b < NM; b++) nbp[b] = basep[b] = 0;
        float2* const cand = nc <= a.cand_cap ? s_cand : lvlg;

        for (int k = 1; k <= T; k++) {
            // ---- the other agents at level k in the agent frame (all of them are visible), and the collision
            //      candidates for the step k-1 -> k (segments meeting the bbox of this level's nodes) ----
            __syncthreads();
            if (tid == 0) s_n = 0;
            __syncthreads();
            {
                const float* lb = t.lvl_bb + 4 * k;
                const float m2 = 2.f * R;
                const float* ax = a.QX + (size_t)(k - 1) * N; const float* ay = a.QY + (size_t)(k - 1) * N;
                const float* bx = a.QX + (size_t)k * N;       const float* by = a.QY + (size_t)k * N;
                for (int j = tid; j < N; j += THREADS) {
                    if (j == i) continue;
                    const float dxa = ax[j] - x0, dya = ay[j] - y0, dxb = bx[j] - x0, dyb = by[j] - y0;
                    const float pax = c0 * dxa + s0 * dya, pay = -s0 * dxa + c0 * dya;
                    const float pbx = c0 * dxb + s0 * dyb, pby = -s0 * dxb + c0 * dyb;
                    cand[j < i ? j : j - 1] = make_float2(pbx, pby);
                    if (fmaxf(pax, pbx) < lb[0] - m2 || fminf(pax, pbx) > lb[1] + m2) continue;
                    if (fmaxf(pay, pby) < lb[2] - m2 || fminf(pay, pby) > lb[3] + m2) continue;
                    clv[atomicAdd(&s_n, 1)] = make_float4(pax, pay, pbx, pby);
                }
            }
            __syncthreads();
            const int nclv = s_n;

            for (int tt = t.ltile[k]; tt < t.ltile[k + 1]; tt++) {
                // ---- collision pruning of the edges into this tile ----
                __syncthreads();
                if (tid == 0) s_n = 0;
                __syncthreads();
                {
                    const float m2 = 2.f * R;
                    const float bx0 = t.tile_bb[4 * tt] - m2, bx1 = t.tile_bb[4 * tt + 1] + m2;
                    const float by0 = t.tile_bb[4 * tt + 2] - m2, by1 = t.tile_bb[4 * tt + 3] + m2;
                    for (int q = tid; q < nclv; q += THREADS) {
                        const float4 c = clv[q];
                        if (fmaxf(c.x, c.z) < bx0 || fminf(c.x, c.z) > bx1) continue;
                        if (fmaxf(c.y, c.w) < by0 || fminf(c.y, c.w) > by1) continue;
                        const int s = atomicAdd(&s_n, 1);
                        if (s < cc_cap) s_cc[s] = c; else ctg[s] = c;
                    }
                }
                __syncthreads();
                const int ncc = s_n;
                if (ncc > cc_cap) {  // rare: copy the shared part to global so one pointer serves all
                    for (int q = tid; q < cc_cap; q += THREADS) ctg[q] = s_cc[q];
                    __syncthreads();
                }
                const float4* cc = ncc > cc_cap ? ctg : s_cc;
                for (int e = t.tile_estart[tt] + tid; e < t.tile_estart[tt + 1]; e += THREADS) {
                    const int pa = t.eparent[e], ch = t.echild[e], mv = t.emove[e];
                    int ma[NM];
                    bool any = false;
                    if (k == 1) {
#pragma unroll
                        for (int b = 0; b < NM; b++) ma[b] = (b == mv);
                        any = true;
                    } else {
#pragma unroll
                        for (int b = 0; b < NM; b++) { ma[b] = mult[pa * NM + b]; any |= ma[b] != 0; }
                    }
                    if (!any) continue;
                    const float p0x = t.node_x[pa], p0y = t.node_y[pa], p1x = t.node_x[ch], p1y = t.node_y[ch];
                    bool hit = false;
                    for (int c = 0; c < ncc && !hit; c++) {
                        const float4 q = cc[c];
                        const float r1x = p1x - q.z, r1y = p1y - q.w;
                        if (r1x * r1x + r1y * r1y < a.D2) { hit = true; break; }
                        const float r0x = p0x - q.x, r0y = p0y - q.y;
                        if (r0x * r0x + r0y * r0y < a.D2) continue;  // overlapping at the start: ignore the pair
                        const float dx = r1x - r0x, dy = r1y - r0y, dd = dx * dx + dy * dy;
                        if (dd <= 0.f) continue;
                        const float s = -(r0x * dx + r0y * dy) / dd;
                        if (s > 0.f && s < 1.f) {
                            const float mx = r0x + s * dx, my = r0y + s * dy;
                            if (mx * mx + my * my < a.D2) hit = true;
                        }
                    }
                    if (hit) continue;
#pragma unroll
                    for (int b = 0; b < NM; b++) if (ma[b]) atomicAdd(&mult[ch * NM + b], ma[b]);
                }
                __syncthreads();

                // ---- visual states, one group of positions at a time ----
                for (int gi = t.tile_gstart[tt]; gi < t.tile_gstart[tt + 1]; gi++) {
                    const int p0 = t.grp_p0[gi], u0 = t.grp_u0[gi], nu = t.grp_u1[gi] - u0;
                    const int n0 = t.grp_n0[gi], n1 = t.grp_n1[gi];
                    for (int q = tid; q < GU * NM; q += THREADS) (&s_um[0][0])[q] = 0;
                    if (tid < G) s_greach[tid] = 0;
                    __syncthreads();
                    for (int n = n0 + tid; n < n1; n += THREADS) {  // sum node path counts into their units
                        const int u = t.node_unit[n] - u0;
#pragma unroll
                        for (int b = 0; b < NM; b++) {
                            const int m = mult[n * NM + b];
                            nbp[b] += m;
                            if (m) atomicAdd(&s_um[u][b], m);
                        }
                    }
                    __syncthreads();
                    if (tid < nu) {
                        bool any = false;
#pragma unroll
                        for (int b = 0; b < NM; b++) any |= s_um[tid][b] != 0;
                        const int slot = t.unit_pos[u0 + tid] - p0;
                        s_ush[tid] = t.unit_sh[u0 + tid];
                        s_uslot[tid] = any ? slot : -1;
                        if (any) s_greach[slot] = 1;
                    }
                    __syncthreads();
                    if (tid < 32) {  // compact the reached positions
                        const bool r = tid < G && s_greach[tid];
                        const unsigned bm = __ballot_sync(0xffffffffu, r);
                        if (tid < G) s_slotmap[tid] = r ? __popc(bm & ((1u << tid) - 1u)) : -1;
                        if (tid == 0) s_ns = __popc(bm);
                    }
                    __syncthreads();
                    const int ns = s_ns;
                    if (ns == 0) continue;
                    if (tid < G && s_slotmap[tid] >= 0) {
                        s_px[s_slotmap[tid]] = t.pos_x[p0 + tid];
                        s_py[s_slotmap[tid]] = t.pos_y[p0 + tid];
                    }
                    for (int w = tid; w < ns * MW; w += THREADS) s_mask[w] = 0u;
                    __syncthreads();
                    // ---- paint: the group's (candidate, position) pairs, p = j ns + g, are dealt to the threads in
                    //      turn, so every lane has work whatever ns is; neighbouring lanes share a few candidates
                    //      (broadcast loads). Two independent chains per thread: pairs tid + 2m THREADS and
                    //      tid + (2m + 1) THREADS ----
                    {
                        const int dq = 2 * THREADS / ns, dr = 2 * THREADS % ns;
                        auto step = [&](int& j, int& g) {
                            g += dr;
                            j += dq;
                            if (g >= ns) { g -= ns; j++; }
                        };
                        auto loop = [&](const float2* cv, auto asin_tag) {
                            constexpr bool AS = decltype(asin_tag)::value;
                            int ja = tid / ns, ga = tid - ja * ns;
                            int jb = (tid + THREADS) / ns, gb = tid + THREADS - jb * ns;  // jb >= ja
                            for (; jb < nc; step(ja, ga), step(jb, gb)) {
                                const float2 v = cv[ja], u = cv[jb];
                                paint<B, AS>(s_mask + ga * MW, v.x - s_px[ga], v.y - s_py[ga], R, inv_bw);
                                paint<B, AS>(s_mask + gb * MW, u.x - s_px[gb], u.y - s_py[gb], R, inv_bw);
                            }
                            if (ja < nc) {
                                const float2 v = cv[ja];
                                paint<B, AS>(s_mask + ga * MW, v.x - s_px[ga], v.y - s_py[ga], R, inv_bw);
                            }
                        };
                        auto dispatch = [&](const float2* cv) {
                            if (a.half_asin) loop(cv, std::true_type{}); else loop(cv, std::false_type{});
                        };
                        if (nc <= a.cand_cap) dispatch(s_cand); else dispatch(lvlg);
                    }
                    __syncthreads();
                    // ---- the edges of each position's mask (bin k starts an interval when k is set and k-1 is not,
                    //      ends one when k is clear and k-1 set): counted, then listed; one warp per position ----
                    for (int g = warp; g < ns; g += NW) {
                        const unsigned* bits = s_mask + g * MW;
                        int cnt = 0;
                        for (int w = lane; w < WORDS; w += 32) {
                            const unsigned m = bits[w], prev = bits[w == 0 ? WORDS - 1 : w - 1];
                            cnt += __popc(m ^ ((m << 1) | (prev >> 31)));
                        }
                        for (int o = 16; o > 0; o >>= 1) cnt += __shfl_xor_sync(0xffffffffu, cnt, o);
                        if (lane == 0) s_ecnt[g] = cnt;
                    }
                    __syncthreads();
                    if (tid == 0) {
                        s_eoff[0] = 0;
                        for (int g = 0; g < ns; g++) s_eoff[g + 1] = s_eoff[g] + s_ecnt[g];
                    }
                    __syncthreads();
                    int* elist = s_eoff[ns] <= a.edge_cap ? s_edge : edge_g;
                    for (int g = warp; g < ns; g += NW) {
                        const unsigned* bits = s_mask + g * MW;
                        int base = s_eoff[g];
                        for (int w0 = 0; w0 < WORDS; w0 += 32) {
                            const int w = w0 + lane;
                            unsigned ed = 0u, m = 0u;
                            if (w < WORDS) {
                                m = bits[w];
                                const unsigned prev = bits[w == 0 ? WORDS - 1 : w - 1];
                                ed = m ^ ((m << 1) | (prev >> 31));
                            }
                            const int c = __popc(ed);
                            int incl = c;  // warp inclusive scan of edge counts
                            for (int o = 1; o < 32; o <<= 1) {
                                const int v = __shfl_up_sync(0xffffffffu, incl, o);
                                if (lane >= o) incl += v;
                            }
                            int pos = base + incl - c;
                            while (ed) {  // entry: bin * 2 + (1 if the interval ends here)
                                const int bit = __ffs(ed) - 1;
                                ed &= ed - 1;
                                elist[pos++] = ((w * 32 + bit) << 1) | (int)(((m >> bit) & 1u) ^ 1u);
                            }
                            base += __shfl_sync(0xffffffffu, incl, 31);
                        }
                    }
                    __syncthreads();
                    // ---- units scatter +-path counts at their (heading-rotated) edges; one warp per unit ----
                    scatter_units<B, THREADS>(nu, s_uslot, s_slotmap, s_ush, s_um, s_mask, s_eoff, elist, a.d_shared, Ds, D, basep);
                    __syncthreads();
                }
            }
        }
        // ---- n_b and the base counts, then the scores ----
        __syncthreads();
#pragma unroll
        for (int b = 0; b < NM; b++) {
            int v = nbp[b], w = basep[b];
            for (int o = 16; o > 0; o >>= 1) { v += __shfl_xor_sync(0xffffffffu, v, o); w += __shfl_xor_sync(0xffffffffu, w, o); }
            if (lane == 0) { s_redl[b][warp] = v; s_redi[b][warp] = w; }
        }
        __syncthreads();
        long long nb[NM];
        int base[NM];
#pragma unroll
        for (int b = 0; b < NM; b++) {
            long long s = 0;
            int bs = 0;
            for (int w = 0; w < NW; w++) { s += s_redl[b][w]; bs += s_redi[b][w]; }
            nb[b] = s; base[b] = bs;
        }
        __syncthreads();
        finalize_scores<B, THREADS>(D, nb, base, s_redl, s_redi, a.score + (size_t)i * NM);
    }
}

}  // namespace fsm2d
