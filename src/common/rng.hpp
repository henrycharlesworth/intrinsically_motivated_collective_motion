// Random streams. Every random draw belongs to a stream indexed by (seed, a, b) — for example (seed, step, agent) for
// an agent's tie-break at one step — so results do not depend on the number of threads or GPUs, on the order in
// which agents are processed, or on whether a run was interrupted and resumed.
#pragma once

#include <cstdint>
#include <random>

namespace fsm {

inline uint64_t splitmix64(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

inline std::mt19937_64 make_rng(uint64_t seed, uint64_t a, uint64_t b) {
    return std::mt19937_64(splitmix64(seed ^ splitmix64(a * 0x100000001B3ULL + b + 1)));
}

// Stream indices that are not (step, agent) pairs.
constexpr uint64_t kInitStream = ~0ULL;          // initial conditions: make_rng(seed, kInitStream, 0)
constexpr uint64_t kNoiseSeedXor = 0xABCDEFULL;  // motor noise: make_rng(seed ^ kNoiseSeedXor, step, agent)

}  // namespace fsm
