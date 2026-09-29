// CUDA helpers shared by the GPU simulators: error checking, fast inverse trigonometry and GPU lists.
#pragma once

#include <cuda_runtime.h>

#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

#define CK(x)                                                                                              \
    do {                                                                                                   \
        cudaError_t e_ = (x);                                                                              \
        if (e_ != cudaSuccess) {                                                                           \
            std::fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); \
            std::exit(1);                                                                                  \
        }                                                                                                  \
    } while (0)

namespace fsm {

// --gpus 0,1,.. or all
inline std::vector<int> parse_gpus(const std::string& v) {
    std::vector<int> gpus;
    if (v == "all") {
        int n = 0;
        CK(cudaGetDeviceCount(&n));
        for (int g = 0; g < n; g++) gpus.push_back(g);
    } else {
        std::stringstream ss(v);
        std::string tok;
        while (std::getline(ss, tok, ',')) gpus.push_back(std::stoi(tok));
    }
    return gpus;
}

// atan2 by octant reduction and a fitted odd polynomial: max error 1.7e-6 rad, far below the angular resolution of
// the visual field (2e-3 rad for 3072 bins).
__device__ __forceinline__ float fast_atan2(float y, float x) {
    const float ax = fabsf(x), ay = fabsf(y), mx = fmaxf(ax, ay), mn = fminf(ax, ay);
    const float t = __fdividef(mn, fmaxf(mx, 1e-30f)), s = t * t;
    float r = fmaf(s, -1.171559314e-02f, 5.263565089e-02f);
    r = fmaf(r, s, -1.164120617e-01f);
    r = fmaf(r, s, 1.935322990e-01f);
    r = fmaf(r, s, -3.326208634e-01f);
    r = fmaf(r, s, 9.999770681e-01f);
    r *= t;
    if (ay > ax) r = 1.57079632679f - r;
    if (x < 0.f) r = 3.14159265359f - r;
    return y < 0.f ? -r : r;
}

// The half-angle atan(R / d) of a disc: a polynomial for x = R / d <= 1/2 (any neighbour that does not overlap us;
// max error 5e-7 rad), else the library atanf.
__device__ __forceinline__ float half_atan(float x) {
    if (x > 0.5f) return atanf(x);
    const float s = x * x;
    float r = fmaf(s, -9.578907725e-02f, 1.920534081e-01f);
    r = fmaf(r, s, -3.328121439e-01f);
    r = fmaf(r, s, 9.999906105e-01f);
    return r * x;
}

}  // namespace fsm
