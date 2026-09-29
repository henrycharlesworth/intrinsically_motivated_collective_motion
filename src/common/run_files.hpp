// The output files of a run, written step by step (so a run can be inspected, or resumed, while it is going):
//   PREFIX.json    parameters and progress (steps done, wall time)
//   PREFIX.traj    float64 frames [steps+1][N][frame values]
//   PREFIX.moves   int8 [steps][N]: the move each agent took
//   PREFIX.scores  float32 [steps][N][moves]: each agent's score for every first move
// plus optional extra streams with one record per frame (e.g. the 3D predators) or per step.
//
// Resuming: if PREFIX.json exists with the same parameters, every file is cut back to the last step that is complete
// in all of them and the simulation continues from that step's frame. Because every random draw is indexed by step
// (common/rng.hpp), the result is bit-identical to an uninterrupted run.
#pragma once

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "json.hpp"

namespace fsm {

inline void write_text(const std::string& path, const std::string& s) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    f << s;
}

inline std::string read_text(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

class RunFiles {
public:
    struct Stream {
        std::string suffix;  // e.g. ".pred"
        size_t record;       // bytes per record
        bool per_frame;      // one record per frame (steps + 1 records), else one per step
    };

    // frame_bytes: one frame of .traj; step_moves / step_scores: bytes of .moves / .scores per step.
    RunFiles(std::string prefix, size_t frame_bytes, size_t step_moves, size_t step_scores, std::vector<Stream> extra = {})
        : prefix_(std::move(prefix)) {
        streams_ = {{".traj", frame_bytes, true}, {".moves", step_moves, false}, {".scores", step_scores, false}};
        streams_.insert(streams_.end(), extra.begin(), extra.end());
    }

    // Start the run: returns the first step to simulate (0 for a new run). params_json is the JSON the run will write
    // (its "steps" and "wall_seconds" lines are ignored in the comparison with an existing run). On resume, `frame`
    // receives the last complete frame and wall_seconds the time spent so far.
    int open(const std::string& params_json, bool resume, std::vector<char>& frame, double& wall_seconds) {
        namespace fs = std::filesystem;
        int start = 0;
        wall_seconds = 0;
        const std::string js = prefix_ + ".json";
        const bool have = std::all_of(streams_.begin(), streams_.end(), [&](const Stream& s) { return fs::exists(prefix_ + s.suffix); });
        if (resume && fs::exists(js) && have) {
            const std::string old = read_text(js);
            if (strip_progress(old) != strip_progress(params_json))
                throw std::runtime_error("--resume: " + js + " has different parameters");
            long long steps = -1;  // complete steps in every file
            for (const auto& s : streams_) {
                const long long n = (long long)(fs::file_size(prefix_ + s.suffix) / s.record) - (s.per_frame ? 1 : 0);
                steps = steps < 0 ? n : std::min(steps, n);
            }
            if (steps > 0) {
                start = (int)steps;
                for (const auto& s : streams_) fs::resize_file(prefix_ + s.suffix, (uintmax_t)s.record * (start + (s.per_frame ? 1 : 0)));
                std::ifstream f(prefix_ + ".traj", std::ios::binary);
                f.seekg((std::streamoff)(streams_[0].record * start));
                frame.resize(streams_[0].record);
                f.read(frame.data(), (std::streamsize)frame.size());
                wall_seconds = json_number(old, "wall_seconds", 0);
                std::fprintf(stderr, "resume: %s from step %d\n", prefix_.c_str(), start);
            }
        }
        const auto mode = start > 0 ? std::ios::binary | std::ios::app : std::ios::binary;
        for (const auto& s : streams_) {
            files_.emplace_back(prefix_ + s.suffix, mode);
            if (!files_.back()) throw std::runtime_error("cannot write " + prefix_ + s.suffix);
        }
        if (start == 0) write_json(params_json);
        return start;
    }

    // Stream index k: 0 traj, 1 moves, 2 scores, 3.. the extra streams in the order given.
    void write(int k, const void* data, size_t bytes) {
        files_[k].write((const char*)data, (std::streamsize)bytes);
        files_[k].flush();
    }
    void write_json(const std::string& text) const { write_text(prefix_ + ".json", text); }
    const std::string& prefix() const { return prefix_; }

private:
    static std::string strip_progress(const std::string& js) {
        std::istringstream in(js);
        std::string line, out;
        while (std::getline(in, line))
            if (line.find("\"steps\"") == std::string::npos && line.find("\"wall_seconds\"") == std::string::npos) out += line + "\n";
        return out;
    }

    std::string prefix_;
    std::vector<Stream> streams_;
    std::vector<std::ofstream> files_;
};

}  // namespace fsm
