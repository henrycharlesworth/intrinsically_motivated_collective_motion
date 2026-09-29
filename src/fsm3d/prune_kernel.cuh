// Kernel 1 of the 3D scoring (docs/OPTIMIZATIONS.md, "GPU 3D"): collision pruning. One block per agent (persistent
// blocks pulling agents from a queue): the 2D pruning with swept spheres, per-level and per-tile bounding boxes, and
// the predators (as each agent knows them) as extra candidates with their larger, growing contact distance. Output:
// the collision-free path counts of every unit per first move, um[u][b], and n_b.
#pragma once

#include "device.cuh"

namespace fsm3d {

// ---- 1. collision pruning: unit path counts per first move -------------------------------------------------
template <int THREADS>
__global__ void __launch_bounds__(THREADS) prune_kernel(Args a) {
    const int tid = threadIdx.x, N = a.N, T = a.tau;
    int* mult = a.mult + (size_t)blockIdx.x * a.Mmax * NM;
    float4* clv = a.cl + (size_t)blockIdx.x * 2 * N;  // two float4 per candidate: its positions at k - 1 and k
    float4* ctv = a.ct + (size_t)blockIdx.x * 2 * N;
    constexpr int CC = 2048;
    __shared__ float4 s_cc[CC];
    __shared__ int s_n, s_task;
    __shared__ long long s_nb[NM];
    for (;;) {
        __syncthreads();
        if (tid == 0) s_task = atomicAdd(a.next_task, 1);
        __syncthreads();
        const int slot = s_task;
        if (slot >= a.na) return;
        const int i = a.a0 + slot;
        const DevTemplate& t = a.tpl[a.tix[i]];
        float x0, y0, z0, c0, s0;
        frame_of(a, i, x0, y0, z0, c0, s0);
        for (int q = tid; q < t.M * NM; q += THREADS) mult[q] = 0;
        int* um = a.um + (size_t)slot * a.Umax * NM;
        for (int q = tid; q < t.U * NM; q += THREADS) um[q] = 0;
        if (tid < NM) s_nb[tid] = 0;
        for (int k = 1; k <= T; k++) {
            __syncthreads();
            if (tid == 0) s_n = 0;
            __syncthreads();
            {  // other agents whose step k-1 -> k meets the (inflated) bbox of this level's nodes
                const float* lb = t.lvl_bb + 6 * k;
                const float m2 = 2.f * a.R;
                const float *ax = a.QX + (size_t)(k - 1) * N, *ay = a.QY + (size_t)(k - 1) * N, *az = a.QZ + (size_t)(k - 1) * N;
                const float *bx = a.QX + (size_t)k * N, *by = a.QY + (size_t)k * N, *bz = a.QZ + (size_t)k * N;
                for (int j = tid; j < N; j += THREADS) {
                    if (j == i) continue;
                    const float dxa = ax[j] - x0, dya = ay[j] - y0, dxb = bx[j] - x0, dyb = by[j] - y0;
                    const float pax = c0 * dxa + s0 * dya, pay = -s0 * dxa + c0 * dya, paz = az[j] - z0;
                    const float pbx = c0 * dxb + s0 * dyb, pby = -s0 * dxb + c0 * dyb, pbz = bz[j] - z0;
                    if (fmaxf(pax, pbx) < lb[0] - m2 || fminf(pax, pbx) > lb[1] + m2) continue;
                    if (fmaxf(pay, pby) < lb[2] - m2 || fminf(pay, pby) > lb[3] + m2) continue;
                    if (fmaxf(paz, pbz) < lb[4] - m2 || fminf(paz, pbz) > lb[5] + m2) continue;
                    const int s = atomicAdd(&s_n, 1);
                    clv[2 * s] = make_float4(pax, pay, paz, m2);  // w: bbox margin, then the squared contact distance
                    clv[2 * s + 1] = make_float4(pbx, pby, pbz, a.D2);
                }
                // the predators, as this agent knows them (a thread each): the same test with their (larger,
                // growing) contact distance
                if (tid < a.npred && a.PA[(size_t)tid * N + i] >= 0.f) {
                    const size_t base = (size_t)tid * (T + 1) * N, qa = base + (size_t)(k - 1) * N + i, qb = base + (size_t)k * N + i;
                    const float dxa = a.PX[qa] - x0, dya = a.PY[qa] - y0, dxb = a.PX[qb] - x0, dyb = a.PY[qb] - y0;
                    const float pax = c0 * dxa + s0 * dya, pay = -s0 * dxa + c0 * dya, paz = a.PZ[qa] - z0;
                    const float pbx = c0 * dxb + s0 * dyb, pby = -s0 * dxb + c0 * dyb, pbz = a.PZ[qb] - z0;
                    const float mp = a.pm + a.pgrow * (k + a.PA[(size_t)tid * N + i]);
                    if (!(fmaxf(pax, pbx) < lb[0] - mp || fminf(pax, pbx) > lb[1] + mp || fmaxf(pay, pby) < lb[2] - mp ||
                          fminf(pay, pby) > lb[3] + mp || fmaxf(paz, pbz) < lb[4] - mp || fminf(paz, pbz) > lb[5] + mp)) {
                        const int s = atomicAdd(&s_n, 1);
                        clv[2 * s] = make_float4(pax, pay, paz, mp);
                        clv[2 * s + 1] = make_float4(pbx, pby, pbz, mp * mp);
                    }
                }
            }
            __syncthreads();
            const int nclv = s_n;
            for (int tt = t.ltile[k]; tt < t.ltile[k + 1]; tt++) {
                __syncthreads();
                if (tid == 0) s_n = 0;
                __syncthreads();
                const float* bb = t.tile_bb + 6 * tt;
                for (int q = tid; q < nclv; q += THREADS) {
                    const float4 A0 = clv[2 * q], A1 = clv[2 * q + 1];
                    const float m2 = A0.w;
                    if (fmaxf(A0.x, A1.x) < bb[0] - m2 || fminf(A0.x, A1.x) > bb[1] + m2) continue;
                    if (fmaxf(A0.y, A1.y) < bb[2] - m2 || fminf(A0.y, A1.y) > bb[3] + m2) continue;
                    if (fmaxf(A0.z, A1.z) < bb[4] - m2 || fminf(A0.z, A1.z) > bb[5] + m2) continue;
                    const int s = atomicAdd(&s_n, 1);
                    if (2 * s + 1 < CC) { s_cc[2 * s] = A0; s_cc[2 * s + 1] = A1; }
                    else { ctv[2 * s] = A0; ctv[2 * s + 1] = A1; }
                }
                __syncthreads();
                const int ncc = s_n;
                if (2 * ncc > CC) {  // rare: move the shared part to global so one pointer serves all
                    for (int q = tid; q < CC; q += THREADS) ctv[q] = s_cc[q];
                    __syncthreads();
                }
                const float4* cc = 2 * ncc > CC ? ctv : s_cc;
                for (int e = t.tile_estart[tt] + tid; e < t.tile_estart[tt + 1]; e += THREADS) {
                    const int pa = t.eparent[e], ch = t.echild[e], mv = t.emove[e];
                    int ma[NM];
                    bool any = false;
#pragma unroll
                    for (int b = 0; b < NM; b++) { ma[b] = k == 1 ? (b == mv) : mult[pa * NM + b]; any |= ma[b] != 0; }
                    if (!any) continue;
                    const float p0x = t.node_x[pa], p0y = t.node_y[pa], p0z = t.node_z[pa];
                    const float p1x = t.node_x[ch], p1y = t.node_y[ch], p1z = t.node_z[ch];
                    bool hit = false;
                    for (int c = 0; c < ncc && !hit; c++) {
                        const float4 q0 = cc[2 * c], q1 = cc[2 * c + 1];
                        const float r1x = p1x - q1.x, r1y = p1y - q1.y, r1z = p1z - q1.z;
                        const float D2 = q1.w;
                        if (r1x * r1x + r1y * r1y + r1z * r1z < D2) { hit = true; break; }
                        const float r0x = p0x - q0.x, r0y = p0y - q0.y, r0z = p0z - q0.z;
                        if (r0x * r0x + r0y * r0y + r0z * r0z < D2) continue;  // overlapping at the start: ignore
                        const float dx = r1x - r0x, dy = r1y - r0y, dz = r1z - r0z, dd = dx * dx + dy * dy + dz * dz;
                        if (dd <= 0.f) continue;
                        const float s = -(r0x * dx + r0y * dy + r0z * dz) / dd;
                        if (s > 0.f && s < 1.f) {
                            const float mx = r0x + s * dx, my = r0y + s * dy, mz = r0z + s * dz;
                            if (mx * mx + my * my + mz * mz < D2) hit = true;
                        }
                    }
                    if (hit) continue;
#pragma unroll
                    for (int b = 0; b < NM; b++) if (ma[b]) atomicAdd(&mult[ch * NM + b], ma[b]);
                }
                __syncthreads();
            }
        }
        __syncthreads();
        // units' path counts and n_b
        long long nbp[NM];
#pragma unroll
        for (int b = 0; b < NM; b++) nbp[b] = 0;
        for (int n = 1 + tid; n < t.M; n += THREADS) {
            const int u = t.node_unit[n];
#pragma unroll
            for (int b = 0; b < NM; b++) {
                const int m = mult[n * NM + b];
                if (m) { atomicAdd(&um[u * NM + b], m); nbp[b] += m; }
            }
        }
#pragma unroll
        for (int b = 0; b < NM; b++) {
            long long v = nbp[b];
            for (int o = 16; o > 0; o >>= 1) v += __shfl_xor_sync(0xffffffffu, v, o);
            if ((tid & 31) == 0 && v) atomicAdd((unsigned long long*)&s_nb[b], (unsigned long long)v);
        }
        __syncthreads();
        if (tid < NM) a.nb[(size_t)slot * NM + tid] = s_nb[tid];
    }
}

}  // namespace fsm3d
