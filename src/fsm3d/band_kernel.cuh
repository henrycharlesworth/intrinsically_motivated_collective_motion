// Kernel 2 of the 3D scoring (docs/OPTIMIZATIONS.md, "GPU 3D"): painting and accumulation, one block per (agent,
// band of RB latitude rows); persistent blocks pulling tasks. For the positions reached by any first move, level by
// level and in batches of PB positions:
//   - the candidates whose discs can reach the band are culled (per level, then per batch);
//   - the band's rows of each position's mask are painted into shared bitmaps: a cheap pre-test per (candidate,
//     position) pair, then one longitude interval per disc and row;
//   - every (unit, first move b) of the batch adds +-(its path count) at the mask's row edges, shifted by the unit's
//     azimuth, into a shared difference array of S_b over the band (a mask is painted once for all moves).
// At the end S_b = the prefix sum along each row, and the band's exact integer row sums of S_b (n_b - S_b) are written
// out; the host weights them by the rows' pixel solid angles in a fixed order (deterministic).
// Keeping the difference arrays in shared memory is the point of the band split: per-agent arrays for all moves and
// rows (MBs) only fit in global memory, where the edge atomics are ~50x slower.
#pragma once

#include "device.cuh"

namespace fsm3d {

// ---- 2. band painting and accumulation ------------------------------------------------------------------------
// Grid geometry: longitude pixel q has its centre at (q + 1/2) 2pi / NLON; latitude row j at -pi/2 + (j + 1/2) pi / nlat.
// Shared memory (dynamic; offsets below): per band the difference arrays D[RB][NLON]; per batch of PB positions
// their band bitmaps, positions and units, the kept candidates, and a list of L entries holding first the (disc,
// row) painting tasks, then (once painted) the masks' edges. Lists that overflow fall back to slower direct paths.
template <int NLON, int RB, int PB, int THREADS, int L>
struct BandCfg {
    static constexpr int W = NLON / 32;          // words per row
    static constexpr int WS = W;                 // row stride of the bitmaps
    static constexpr int ENT = 16 * PB;          // (unit, first move) weights of a batch (more: slow path)
    static constexpr int CAND = 512;             // kept candidates per batch
    static constexpr int ECAP = L, TCAP = L;     // edges / painting tasks per batch (the same buffer, in turn)
    static constexpr int LCAP = PB + THREADS;    // pending reached positions
    static constexpr int NW = THREADS / 32;
    // byte offsets (all 16-byte aligned)
    static constexpr size_t oD = 0, oBits = oD + (size_t)NM * RB * NLON * 4, oBuf = oBits + (size_t)PB * RB * WS * 4 + 12;
    static constexpr size_t oBuf16 = (oBuf + 15) / 16 * 16, oTask = oBuf16 + (size_t)CAND * 16;
    static constexpr size_t oEnt = oTask + (size_t)TCAP * 4, oList = oEnt + (size_t)ENT * 8;
    static constexpr size_t oPos = oList + (size_t)LCAP * 4, oEoff = oPos + (size_t)PB * 16, oPp = oEoff + (size_t)(PB + 1) * 4;
    static constexpr size_t oCb = (oPp + (size_t)PB * 4 + 7) / 8 * 8, bytes = oCb + (size_t)CAND * 8;
};

// Reserve c slots of a shared list for this lane, one atomic per warp (the lanes' counts are prefix-summed with
// shuffles). All 32 lanes must call it (c = 0 for none); returns the lane's first slot.
__device__ __forceinline__ int warp_reserve(int* counter, int c) {
    const int lane = threadIdx.x & 31;
    int incl = c;
    for (int o = 1; o < 32; o <<= 1) { const int v = __shfl_up_sync(0xffffffffu, incl, o); if (lane >= o) incl += v; }
    const int total = __shfl_sync(0xffffffffu, incl, 31);
    int base = 0;
    if (lane == 31 && total) base = atomicAdd(counter, total);
    return __shfl_sync(0xffffffffu, base, 31) + incl - c;
}

// first/last integer q with centre (q + 1/2) inside [x0, x1] (in pixel units), i.e. ceil(x0 - 1/2), floor(x1 - 1/2)
__device__ __forceinline__ int first_px(float x) { return (int)ceilf(x - 0.5f); }
__device__ __forceinline__ int last_px(float x) { return (int)floorf(x - 0.5f); }

// Paint the longitude interval of a disc (centre longitude lonc, half-width dl) on one bitmap row (W words).
template <int NLON>
__device__ __forceinline__ void paint_row(unsigned* bits, float lonc, float dl, float inv_px_lon) {
    constexpr int W = NLON / 32;
    if (dl >= (float)M_PI) { for (int w = 0; w < W; w++) bits[w] = 0xffffffffu; return; }  // benign race: same value
    int k0 = first_px((lonc - dl) * inv_px_lon), k1 = last_px((lonc + dl) * inv_px_lon);
    if (k1 < k0) return;
    if (k1 - k0 >= NLON - 1) { for (int w = 0; w < W; w++) bits[w] = 0xffffffffu; return; }
    k0 += NLON; k1 += NLON;  // non-negative (lonc - dl > -2pi)
    const int w0 = k0 >> 5, w1 = k1 >> 5;
    const unsigned mf = 0xffffffffu << (k0 & 31), ml = 0xffffffffu >> (31 - (k1 & 31));
    int wa = w0 % W;
    if (w0 == w1) { atomicOr(&bits[wa], mf & ml); return; }
    atomicOr(&bits[wa], mf);
    for (int w = w0 + 1; w < w1; w++) { wa = wa + 1 == W ? 0 : wa + 1; atomicOr(&bits[wa], 0xffffffffu); }
    wa = wa + 1 == W ? 0 : wa + 1;
    atomicOr(&bits[wa], ml);
}

// Geometry of a disc seen from a position: offset (dx, dy, dz), squared distance d2. The half-angle's sine and
// cosine (no trig), and the sine / cosine of the centre's latitude.
struct Disc { float sa, ca, sl, cl, id; };
__device__ __forceinline__ Disc disc_of(float d2, float dz, float R, int half_asin) {
    Disc o;
    if (half_asin) {
        if (R * R >= d2) { o.sa = 1.f; o.ca = 0.f; }
        else { o.sa = R * rsqrtf(d2); o.ca = sqrtf(1.f - o.sa * o.sa); }
    } else {
        const float ih = rsqrtf(d2 + R * R);
        o.sa = R * ih; o.ca = sqrtf(d2) * ih;
    }
    o.id = rsqrtf(d2);
    o.sl = dz * o.id;
    o.cl = sqrtf(fmaxf(0.f, 1.f - o.sl * o.sl));
    return o;
}
// longitude half-width of the disc on a row of latitude (sin, cos) = (sr, cr); < 0: the row misses it
__device__ __forceinline__ float row_halfwidth(const Disc& o, float sr, float cr) {
    const float num = o.ca - sr * o.sl, den = cr * o.cl;
    if (den <= 1e-12f) return num > 0.f ? -1.f : (float)M_PI;
    const float rr = num / den;
    if (rr > 1.f) return -1.f;
    return rr <= -1.f ? (float)M_PI : fsm::fast_atan2(sqrtf(fmaxf(0.f, 1.f - rr * rr)), rr);  // acos
}

template <int NLON, int RB, int PB, int THREADS, int L>
__global__ void __launch_bounds__(THREADS) band_kernel(Args a, int nlat) {
    using C = BandCfg<NLON, RB, PB, THREADS, L>;
    constexpr int W = C::W, WS = C::WS, CAND = C::CAND;
    const int tid = threadIdx.x, N = a.N, T = a.tau;
    extern __shared__ __align__(16) unsigned char smem[];
    int (*const s_D)[RB][NLON] = (int (*)[RB][NLON])(smem + C::oD);             // S_b differences along each row, per b
    unsigned (*const s_bits)[RB][WS] = (unsigned (*)[RB][WS])(smem + C::oBits);  // band rows of the batch's masks
    float4* const s_cand = (float4*)(smem + C::oBuf16);  // candidates kept for the batch (agent frame), then
    int* const s_edge = (int*)(smem + C::oTask);         // edges (after painting): g << 24 | r << 16 | pixel << 1 | end
    unsigned* const s_ptask = (unsigned*)(smem + C::oTask);  // painting tasks: g << 26 | r << 22 | candidate
    int2* const s_ent = (int2*)(smem + C::oEnt);         // the batch's (unit shift | b << 12 | g << 16, weight), by g
    int* const s_list = (int*)(smem + C::oList);         // reached positions pending a batch
    float4* const s_pos = (float4*)(smem + C::oPos);
    int* const s_eoff = (int*)(smem + C::oEoff);         // [PB + 1] entry range of each position
    int* const s_pp = (int*)(smem + C::oPp);             // the batch's positions (template index)
    float2* const s_cb = (float2*)(smem + C::oCb);       // per painted candidate: bounds on the sine of its latitude
    __shared__ int s_base[NM][RB];                       // S_b at longitude pixel 0 of each row
    __shared__ long long s_nbb[NM];
    __shared__ int s_nc1, s_nk, s_task, s_nl, s_ne, s_nt, s_np, s_nr2, s_npair;
    __shared__ int s_wc[C::NW];
    __shared__ float s_bc[4];                            // batch centre and radius
    __shared__ float s_sinr[RB], s_cosr[RB];
    float4* cand1 = a.cand + (size_t)blockIdx.x * 3 * N;  // stage-1 candidates of the level (global)
    float4 *candA = cand1 + N, *candB = cand1 + 2 * N;    // stage-2 candidates still to paint (global, ping-pong)
    unsigned* const g_ovf = a.ovf + (size_t)blockIdx.x * OVF;  // list entries beyond shared memory (global)
    unsigned short* const g_pairs = a.pairs + (size_t)blockIdx.x * PAIRS;  // pairs past the pre-test: j << 6 | g
    const float inv_px_lon = NLON / (2.f * (float)M_PI), px_lat = (float)M_PI / nlat, inv_px_lat = nlat / (float)M_PI;
    for (;;) {
        __syncthreads();
        if (tid == 0) s_task = atomicAdd(a.next_task, 1);
        __syncthreads();
        const int task = s_task;
        const long long t_task = clock64();
        if (task >= a.na * a.nbands) return;
        const int band = task % a.nbands, slot = task / a.nbands;
        const int r0 = band * RB, nrows = min(RB, nlat - r0);
        long long* out = a.part + (size_t)slot * NM * nlat + r0;  // exact row sums of S (n - S), [b][nlat]
        bool anyb = false;
        for (int b = 0; b < NM; b++) anyb |= a.nb[(size_t)slot * NM + b] != 0;
        if (!anyb) { for (int q = tid; q < NM * nrows; q += THREADS) out[(q / nrows) * nlat + q % nrows] = 0; continue; }
        if (tid < NM) s_nbb[tid] = a.nb[(size_t)slot * NM + tid];
        const int i = a.a0 + slot;
        const DevTemplate& t = a.tpl[a.tix[i]];
        const int* um = a.um + (size_t)slot * a.Umax * NM;
        float x0, y0, z0, c0, s0;
        frame_of(a, i, x0, y0, z0, c0, s0);
        const float lat_lo = -0.5f * (float)M_PI + r0 * px_lat, lat_hi = lat_lo + nrows * px_lat;  // band edges
        // first and last row centres, for the trig-free test "can this disc cover a row centre of the band"
        float sinlo, coslo, sinhi, coshi;
        sincosf(lat_lo + 0.5f * px_lat, &sinlo, &coslo);
        sincosf(lat_hi - 0.5f * px_lat, &sinhi, &coshi);
        for (int q = tid; q < NM * RB * NLON; q += THREADS) (&s_D[0][0][0])[q] = 0;
        for (int q = tid; q < NM * RB; q += THREADS) (&s_base[0][0])[q] = 0;
        if (tid < RB) {
            const float lat = -0.5f * (float)M_PI + (r0 + tid + 0.5f) * px_lat;
            s_sinr[tid] = sinf(lat); s_cosr[tid] = cosf(lat);
        }
        // May a disc at distance d (> rad + 2R) from a centre reach the band from within rad of that centre?
        auto reaches = [&](float ez, float d, float rad) {
            if (d <= rad + 2.f * a.R) return true;
            const float delta = asinf(rad / d) + atanf(a.R / (d - rad)) + 1e-3f;
            const float latc = asinf(fminf(1.f, fmaxf(-1.f, ez / d)));
            return latc + delta >= lat_lo - 0.5f * px_lat && latc - delta <= lat_hi + 0.5f * px_lat;
        };
        for (int k = 1; k <= T; k++) {
            const int P0 = t.lpos[k], P1 = t.lpos[k + 1];
            // stage 1: candidates at level k whose disc can reach the band from anywhere in this level's bbox
            __syncthreads();
            if (tid == 0) { s_nc1 = 0; s_nl = 0; }
            __syncthreads();
            {
                const float* lb = t.lvl_bb + 6 * k;
                const float cx = 0.5f * (lb[0] + lb[1]), cy = 0.5f * (lb[2] + lb[3]), cz = 0.5f * (lb[4] + lb[5]);
                const float rad = 0.5f * sqrtf((lb[1] - lb[0]) * (lb[1] - lb[0]) + (lb[3] - lb[2]) * (lb[3] - lb[2]) + (lb[5] - lb[4]) * (lb[5] - lb[4]));
                const float *qx = a.QX + (size_t)k * N, *qy = a.QY + (size_t)k * N, *qz = a.QZ + (size_t)k * N;
                for (int j0 = tid & ~31; j0 < N; j0 += THREADS) {  // warp-uniform trips (warp_reserve)
                    const int j = j0 + (tid & 31);
                    bool keep = false;
                    float px = 0, py = 0, pz = 0;
                    if (j < N && j != i) {
                        const float dx = qx[j] - x0, dy = qy[j] - y0;
                        px = c0 * dx + s0 * dy; py = -s0 * dx + c0 * dy; pz = qz[j] - z0;
                        const float ex = px - cx, ey = py - cy, ez = pz - cz;
                        keep = reaches(ez, sqrtf(ex * ex + ey * ey + ez * ez), rad);
                    }
                    const int sl = warp_reserve(&s_nc1, keep);
                    if (keep) cand1[sl] = make_float4(px, py, pz, 0.f);
                }
            }
            __syncthreads();
            const int nc1 = s_nc1;
            // one batch: the first np pending positions
            auto run_batch = [&](int np) {
                // +-w of every (unit, first move) of position g at an edge of its mask (root pixel col), in the
                // unit's own frame (pixel col - sh; its pixel 0 is carried by the base count instead)
                auto scatter_edge = [&](int g, int r, int col, int end) {
                    if (s_eoff[np] <= C::ENT) {
                        for (int e = s_eoff[g]; e < s_eoff[g + 1]; e++) {
                            const int2 en = s_ent[e];
                            int cc = col - (en.x & 0xfff);
                            if (cc < 0) cc += NLON;
                            if (cc) atomicAdd(&s_D[(en.x >> 12) & 7][r][cc], end ? -en.y : en.y);
                        }
                    } else {
                        const int pp = s_pp[g];
                        for (int u = t.pos_u0[pp]; u < t.pos_u0[pp + 1]; u++) {
                            int cc = col - t.unit_sh[u];
                            if (cc < 0) cc += NLON;
                            if (!cc) continue;
                            for (int b = 0; b < NM; b++) { const int w = um[u * NM + b]; if (w) atomicAdd(&s_D[b][r][cc], end ? -w : w); }
                        }
                    }
                };
                if (tid < np) {
                    const int pp = s_list[tid];
                    s_pp[tid] = pp;
                    s_pos[tid] = make_float4(t.pos_x[pp], t.pos_y[pp], t.pos_z[pp], 0.f);
                    int ne = 0;
                    for (int u = t.pos_u0[pp]; u < t.pos_u0[pp + 1]; u++)
                        for (int b = 0; b < NM; b++) ne += um[u * NM + b] != 0;
                    s_eoff[tid + 1] = ne;
                }
                for (int q = tid; q < PB * RB * WS; q += THREADS) (&s_bits[0][0][0])[q] = 0u;
                if (tid == 0) { s_nk = 0; s_nt = 0; s_ne = 0; }
                __syncthreads();
                if (tid == 0) {  // entry ranges (np <= PB: serial)
                    s_eoff[0] = 0;
                    for (int g = 0; g < np; g++) s_eoff[g + 1] += s_eoff[g];
                }
                __syncthreads();
                if (tid < np && s_eoff[np] <= C::ENT) {  // the (unit, first move) weights, by position
                    const int pp = s_list[tid];
                    int e = s_eoff[tid];
                    for (int u = t.pos_u0[pp]; u < t.pos_u0[pp + 1]; u++)
                        for (int b = 0; b < NM; b++) {
                            const int w = um[u * NM + b];
                            if (w) s_ent[e++] = make_int2(t.unit_sh[u] | b << 12 | tid << 16, w);
                        }
                }
                if (tid < 32) {  // batch centre and radius (bounding box), one warp
                    float mn[3] = {1e30f, 1e30f, 1e30f}, mx[3] = {-1e30f, -1e30f, -1e30f};
                    for (int g = tid; g < np; g += 32) {
                        const float4 q = s_pos[g];
                        mn[0] = fminf(mn[0], q.x); mx[0] = fmaxf(mx[0], q.x); mn[1] = fminf(mn[1], q.y); mx[1] = fmaxf(mx[1], q.y);
                        mn[2] = fminf(mn[2], q.z); mx[2] = fmaxf(mx[2], q.z);
                    }
                    for (int o = 16; o > 0; o >>= 1)
                        for (int d = 0; d < 3; d++) {
                            mn[d] = fminf(mn[d], __shfl_xor_sync(0xffffffffu, mn[d], o));
                            mx[d] = fmaxf(mx[d], __shfl_xor_sync(0xffffffffu, mx[d], o));
                        }
                    if (tid == 0) {
                        for (int d = 0; d < 3; d++) s_bc[d] = 0.5f * (mn[d] + mx[d]);
                        s_bc[3] = 0.5f * sqrtf((mx[0] - mn[0]) * (mx[0] - mn[0]) + (mx[1] - mn[1]) * (mx[1] - mn[1]) + (mx[2] - mn[2]) * (mx[2] - mn[2]));
                    }
                }
                __syncthreads();
                // stage 2: candidates that can reach the band from this batch (to the block's global list)
                const float cx = s_bc[0], cy = s_bc[1], cz = s_bc[2], rad = s_bc[3];
                for (int q0 = tid & ~31; q0 < nc1; q0 += THREADS) {
                    const int q = q0 + (tid & 31);
                    bool keep = false;
                    float4 c = make_float4(0.f, 0.f, 0.f, 0.f);
                    if (q < nc1) {
                        c = cand1[q];
                        const float ex = c.x - cx, ey = c.y - cy, ez = c.z - cz;
                        keep = reaches(ez, sqrtf(ex * ex + ey * ey + ez * ez), rad);
                    }
                    const int sl = warp_reserve(&s_nk, keep);
                    if (keep) candA[sl] = c;
                }
                __syncthreads();
                const int nk = s_nk;
                // The kept candidates are painted in rounds of at most CAND (the shared list's capacity).
                float4 *cR = candA, *cR2 = candB;
                int nr = nk;
                while (nr > 0) {
                    __syncthreads();
                    if (tid == 0) { s_np = 0; s_nr2 = 0; }
                    __syncthreads();
                    for (int q0 = tid & ~31; q0 < nr; q0 += THREADS) {
                        const int q = q0 + (tid & 31);
                        float4 c = make_float4(0.f, 0.f, 0.f, 0.f);
                        int cls = 0;  // 0 none, 1 paint in this round, 2 in a later one
                        float d = 0.f;
                        if (q < nr) {
                            c = cR[q];
                            const float ex = c.x - cx, ey = c.y - cy, ez = c.z - cz;
                            d = sqrtf(ex * ex + ey * ey + ez * ez);
                            cls = 1;
                        }
                        const int sp = warp_reserve(&s_np, cls == 1);
                        if (cls == 1 && sp >= CAND) cls = 2;
                        const int sr = warp_reserve(&s_nr2, cls == 2);
                        if (cls == 1) {
                            s_cand[sp] = c;
                            // sine-of-latitude window of the pair pre-test: the disc's half-angle from any position
                            // of the batch is at most alm (distance >= d - rad)
                            float2 cb = make_float2(-2.f, 2.f);
                            if (d > rad + 1e-3f) {
                                const float ar = a.R / (d - rad);
                                const float alm = a.half_asin ? asinf(fminf(1.f, ar)) : atanf(ar);
                                const float la = lat_lo + 0.5f * px_lat - alm, lb = lat_hi - 0.5f * px_lat + alm;
                                if (la > -0.5f * (float)M_PI) cb.x = sinf(la) - 1e-4f;
                                if (lb < 0.5f * (float)M_PI) cb.y = sinf(lb) + 1e-4f;
                            }
                            s_cb[sp] = cb;
                        }
                        else if (cls == 2) cR2[sr] = c;
                    }
                    __syncthreads();
                    const int npaint = min(s_np, CAND);
                    nr = s_nr2;
                    { float4* tq = cR; cR = cR2; cR2 = tq; }
                    if (npaint == 0) continue;
                    // paint, pass 1: every (candidate, position) pair: a trig-free test whether the disc can cover
                    // a row centre of the band (latitude of its centre within half-angle al of the first / last row,
                    // compared as sines; a side whose angle passes a pole holds), then its rows, one task per row.
                    const float4* cl = s_cand;
                    const int ncand = npaint;
                    // pass 1a: every (candidate, position) pair: a cheap test on the sine of the latitude (the
                    // candidate's window), compacting the pairs that pass. A thread keeps one position in registers
                    // (lanes: consecutive positions) and steps through the candidates (broadcast loads).
                    if (tid == 0) s_npair = 0;
                    __syncthreads();
                    {
                        const int G = (np + 31) & ~31, groups = THREADS / G, g = tid % G, grp = tid / G;
                        const float4 q = s_pos[min(g, np - 1)];
                        for (int j = grp; j < ncand; j += groups) {
                            const float4 c = cl[j];
                            const float2 cb = s_cb[j];
                            const float dx = c.x - q.x, dy = c.y - q.y, dz = c.z - q.z, d2 = dx * dx + dy * dy + dz * dz;
                            const float sl = dz * rsqrtf(d2);
                            const bool on = g < np && d2 > 0.f && sl >= cb.x && sl <= cb.y;
                            const unsigned bal = __ballot_sync(0xffffffffu, on);
                            if (!bal) continue;
                            int base = 0;
                            if ((tid & 31) == 0) base = atomicAdd(&s_npair, __popc(bal));
                            base = __shfl_sync(0xffffffffu, base, 0);
                            if (on) g_pairs[base + __popc(bal & ((1u << (tid & 31)) - 1u))] = (unsigned short)(j << 6 | g);
                        }
                    }
                    __syncthreads();
                    const int npair = s_npair;
                    // pass 1b: the pairs that passed: the exact test and their rows, one task per row
                    for (int q0 = tid & ~31; q0 < npair; q0 += THREADS) {
                        const int q = q0 + (tid & 31);
                        int j = 0, g = 0;
                        int j0 = 0, j1 = -1;
                        Disc o{};
                        float dx = 0, dy = 0;
                        if (q < npair) {
                            const int e = g_pairs[q];
                            j = e >> 6; g = e & 63;
                            const float4 c = cl[j], qp = s_pos[g];
                            dx = c.x - qp.x; dy = c.y - qp.y;
                            const float dz = c.z - qp.z, d2 = dx * dx + dy * dy + dz * dz;
                            if (d2 > 0.f) {
                                o = disc_of(d2, dz, a.R, a.half_asin);
                                const float sx = sinlo * o.ca - coslo * o.sa, cx = coslo * o.ca + sinlo * o.sa;  // lo - al
                                const float sy = sinhi * o.ca + coshi * o.sa, cy = coshi * o.ca - sinhi * o.sa;  // hi + al
                                if (((cx < 0.f && sx < 0.f) || o.sl >= sx - 1e-5f) && (cy < 0.f || o.sl <= sy + 1e-5f)) {
                                    const float al = a.half_asin ? (o.sa >= 1.f ? 0.5f * (float)M_PI : fsm::fast_atan2(o.sa, o.ca)) : fsm::half_atan(a.R * o.id);
                                    const float latc = fsm::fast_atan2(o.sl, o.cl);
                                    j0 = max(r0, first_px((latc - al + 0.5f * (float)M_PI) * inv_px_lat));
                                    j1 = min(r0 + nrows - 1, last_px((latc + al + 0.5f * (float)M_PI) * inv_px_lat));
                                }
                            }
                        }
                        const int cnt = max(0, j1 - j0 + 1);
                        if (!__any_sync(0xffffffffu, cnt)) continue;
                        int s = warp_reserve(&s_nt, cnt);
                        for (int jj = j0; jj <= j1; jj++, s++) {
                            const unsigned tv = (unsigned)g << 26 | (unsigned)(jj - r0) << 22 | (unsigned)j;
                            if (s < C::TCAP) { s_ptask[s] = tv; continue; }
                            if (s < C::TCAP + OVF) { g_ovf[s - C::TCAP] = tv; continue; }
                            const float dl = row_halfwidth(o, s_sinr[jj - r0], s_cosr[jj - r0]);  // beyond both: paint here
                            if (dl >= 0.f) paint_row<NLON>(&s_bits[g][jj - r0][0], fsm::fast_atan2(dy, dx), dl, inv_px_lon);
                        }
                    }
                    __syncthreads();
                    // pass 2: the (disc, row) tasks, dealt to all threads
                    const int nt = min(s_nt, C::TCAP + OVF);
                    for (int q = tid; q < nt; q += THREADS) {
                        const unsigned tv = q < C::TCAP ? s_ptask[q] : g_ovf[q - C::TCAP];
                        const int g = tv >> 26, r = (tv >> 22) & 15, j = tv & 0x3fffff;
                        const float4 c = cl[j], p = s_pos[g];
                        const float dx = c.x - p.x, dy = c.y - p.y, dz = c.z - p.z;
                        const Disc o = disc_of(dx * dx + dy * dy + dz * dz, dz, a.R, a.half_asin);
                        const float dl = row_halfwidth(o, s_sinr[r], s_cosr[r]);
                        if (dl >= 0.f) paint_row<NLON>(&s_bits[g][r][0], fsm::fast_atan2(dy, dx), dl, inv_px_lon);
                    }
                    __syncthreads();
                    if (tid == 0) s_nt = 0;
                }
                __syncthreads();
                // scatter. Every (position, row, word) finds its interval edges (bit q set and q-1 clear: an
                // interval starts at q; clear and set: one ends) and lists them (slots reserved with one shared
                // atomic per word; the order doesn't matter, the sums are integers); then the edges are dealt to
                // all threads, and each adds +-w for every unit of its position at the edge's pixel in the unit's
                // own frame (root pixel q is the unit's pixel q - sh; its pixel 0 is carried by the base count
                // instead). Edges beyond the list's capacity are scattered directly.
                for (int q0 = tid & ~31; q0 < np * RB * W; q0 += THREADS) {  // warp-uniform trips
                    const int q = q0 + (tid & 31);
                    const int w = q % W, r = (q / W) % RB, g = q / (W * RB);
                    unsigned m = 0, ed = 0;
                    if (q < np * RB * W && r < nrows) {
                        const unsigned* bits = &s_bits[g][r][0];
                        m = bits[w];
                        ed = m ^ ((m << 1) | (bits[w == 0 ? W - 1 : w - 1] >> 31));
                    }
                    int slot = warp_reserve(&s_ne, __popc(ed));
                    while (ed) {
                        const int bit = __ffs(ed) - 1;
                        ed &= ed - 1;
                        const int col = w * 32 + bit, end = ((m >> bit) & 1u) ^ 1u;
                        const int ev = g << 24 | r << 16 | col << 1 | end;
                        if (slot < C::ECAP) { s_edge[slot++] = ev; continue; }
                        if (slot < C::ECAP + OVF) { g_ovf[slot++ - C::ECAP] = (unsigned)ev; continue; }
                        scatter_edge(g, r, col, end);
                    }
                }
                __syncthreads();
                const int ne = min(s_ne, C::ECAP + OVF);
                for (int e = tid; e < ne; e += THREADS) {
                    const int ev = e < C::ECAP ? s_edge[e] : (int)g_ovf[e - C::ECAP], g = ev >> 24, r = (ev >> 16) & 0xff, col = (ev >> 1) & 0x7fff;
                    scatter_edge(g, r, col, ev & 1);
                }
                // base: the unit covers its own pixel 0 (root pixel sh)
                const int nent = s_eoff[np];
                if (nent <= C::ENT) {
                    for (int q = tid; q < nent * nrows; q += THREADS) {
                        const int e = q / nrows, r = q % nrows;
                        const int2 en = s_ent[e];
                        const int sh = en.x & 0xfff, g = en.x >> 16;
                        if ((s_bits[g][r][sh >> 5] >> (sh & 31)) & 1u) atomicAdd(&s_base[(en.x >> 12) & 7][r], en.y);
                    }
                } else {
                    for (int q = tid; q < np * nrows; q += THREADS) {
                        const int g = q / nrows, r = q % nrows, pp = s_pp[g];
                        for (int u = t.pos_u0[pp]; u < t.pos_u0[pp + 1]; u++) {
                            const int sh = t.unit_sh[u];
                            if (!((s_bits[g][r][sh >> 5] >> (sh & 31)) & 1u)) continue;
                            for (int b = 0; b < NM; b++) { const int w = um[u * NM + b]; if (w) atomicAdd(&s_base[b][r], w); }
                        }
                    }
                }
                __syncthreads();
            };
            // positions of the level reached by move b, gathered a chunk at a time into the pending list
            for (int c0p = P0;; c0p += THREADS) {  // (one call site of run_batch, so that it is inlined)
                const bool last = c0p + THREADS >= P1;
                const int pp = c0p + tid;
                bool on = false;
                if (pp < P1)
                    for (int q = t.pos_u0[pp] * NM; q < t.pos_u0[pp + 1] * NM && !on; q++) on = um[q] != 0;
                const unsigned bal = __ballot_sync(0xffffffffu, on);  // order-preserving compaction
                if ((tid & 31) == 0) s_wc[tid >> 5] = __popc(bal);
                __syncthreads();
                int off = s_nl;
                for (int w = 0; w < (tid >> 5); w++) off += s_wc[w];
                if (on) s_list[off + __popc(bal & ((1u << (tid & 31)) - 1u))] = pp;
                __syncthreads();
                if (tid == THREADS - 1) s_nl = off + __popc(bal);
                __syncthreads();
                while (s_nl >= PB || (last && s_nl > 0)) {  // full batches (then the level's last, partial one);
                    const int np = min(s_nl, PB);             // the remainder moves to the front of the list
                    run_batch(np);
                    const int rest = s_nl - np;
                    int v[(C::LCAP + THREADS - 1) / THREADS];
                    for (int q = tid, z = 0; q < rest; q += THREADS, z++) v[z] = s_list[np + q];
                    __syncthreads();
                    for (int q = tid, z = 0; q < rest; q += THREADS, z++) s_list[q] = v[z];
                    if (tid == 0) s_nl = rest;
                    __syncthreads();
                }
                if (last) break;
            }
        }
        __syncthreads();
        // S_b along each row (warp per row: lanes own NLON / 32 consecutive pixels), and sum S (n - S) exactly
        for (int br = tid >> 5; br < NM * nrows; br += C::NW) {
            const int b = br / nrows, r = br % nrows, lane = tid & 31, per = NLON / 32;
            const long long nbb = s_nbb[b];
            if (nbb == 0) { if (lane == 0) out[(size_t)b * nlat + r] = 0; continue; }
            int local = 0;
            for (int q = 0; q < per; q++) { const int cc = lane * per + q; local += cc ? s_D[b][r][cc] : 0; }
            int incl = local;
            for (int o = 1; o < 32; o <<= 1) { const int v = __shfl_up_sync(0xffffffffu, incl, o); if (lane >= o) incl += v; }
            long long S = s_base[b][r] + incl - local, acc = 0;
            for (int q = 0; q < per; q++) {
                const int cc = lane * per + q;
                if (cc) S += s_D[b][r][cc];
                acc += S * (nbb - S);
            }
            for (int o = 16; o > 0; o >>= 1) acc += __shfl_xor_sync(0xffffffffu, acc, o);
            if (lane == 0) out[(size_t)b * nlat + r] = acc;
        }
        if (tid == 0) atomicAdd(&a.cost[slot], (unsigned long long)(clock64() - t_task));
    }
}

}  // namespace fsm3d
