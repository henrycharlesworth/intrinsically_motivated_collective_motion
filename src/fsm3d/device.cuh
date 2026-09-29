// Device-side data shared by the 3D kernels: the templates, the kernel arguments and per-block list capacities.
#pragma once

#include "../common/cuda_util.cuh"
#include "model.hpp"

namespace fsm3d {

constexpr int OVF = 1 << 16;  // per band block: entries of the task / edge lists beyond shared memory
constexpr int PAIRS = 512 * 64;  // per band block: (candidate, position) pairs of a round (CAND x PB at most)

struct DevTemplate {
    int M, P, U;
    const int *ltile, *lpos, *node_unit, *unit_sh, *pos_u0, *pos_n0, *tile_estart, *eparent, *echild, *emove;
    const float *node_x, *node_y, *node_z, *pos_x, *pos_y, *pos_z, *tile_bb, *lvl_bb;
};

struct Args {
    const DevTemplate* tpl;     // one per elevation index
    int tau, N, a0, na;         // agents a0 .. a0 + na - 1 in this pass
    float R, D2;
    int half_asin;
    const float *QX, *QY, *QZ;  // futures of every agent, level-major [(tau+1) N], shared world frame
    const float* TH;            // azimuths
    int npred;                  // predators on (some agents may still not know of them)
    const float *PX, *PY, *PZ;  // each agent's prediction of each, [npred][(tau+1)][N], same frame as QX..
    const float* PA;            // [npred][N] age of that information in steps (-1: the agent does not know of it)
    float pm, pgrow;            // lethal distance of a predator for futures: pm + pgrow * level
    const int* tix;             // template index of every agent
    int Umax, Mmax;
    int* mult;                  // prune: [grid][Mmax][NM]
    float4 *cl, *ct;            // prune: [grid][2 N] level / tile collision candidates
    int* um;                    // [na][Umax][NM] unit path counts per first move
    long long* nb;              // [na][NM]
    int nbands;
    float4* cand;               // band: [grid][3 N] candidates reaching the band from the level (stage 1), then the batch
    unsigned* ovf;              // band: [grid][OVF] overflow of the batch's task / edge lists (global, L2)
    unsigned short* pairs;      // band: [grid][PAIRS] pairs passing the pre-test
    long long* part;            // [na][NM][nlat] exact row sums of S (n - S)
    unsigned long long* cost;   // [na] band-kernel block cycles of each agent
    int* next_task;
};

__device__ __forceinline__ void frame_of(const Args& a, int i, float& x0, float& y0, float& z0, float& c0, float& s0) {
    x0 = a.QX[i]; y0 = a.QY[i]; z0 = a.QZ[i];
    sincosf(a.TH[i], &s0, &c0);
}

}  // namespace fsm3d
