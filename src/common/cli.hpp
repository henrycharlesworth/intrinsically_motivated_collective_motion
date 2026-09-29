// Minimal command-line parsing: `--flag value` pairs and bare `--switch`es, with typed conversion and an error for
// anything unknown or malformed.
#pragma once

#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fsm {

class Args {
public:
    Args(int argc, char** argv) : args_(argv + 1, argv + argc) {}

    bool done() const { return pos_ >= args_.size(); }
    // The next flag; every flag handler must then consume its value (if any) with one of the getters.
    const std::string& next() { return args_[pos_++]; }

    std::string str() {
        if (pos_ >= args_.size()) throw std::runtime_error("missing value for " + args_[pos_ - 1]);
        return args_[pos_++];
    }
    int integer() { return convert<int>([](const std::string& s) { return std::stoi(s); }); }
    uint64_t u64() { return convert<uint64_t>([](const std::string& s) { return std::stoull(s); }); }
    double real() { return convert<double>([](const std::string& s) { return std::stod(s); }); }
    std::vector<double> reals() {  // comma-separated list
        std::vector<double> out;
        for (const auto& item : split(str(), ',')) out.push_back(std::stod(item));
        return out;
    }
    // One of `choices` (for flags like --half-angle atan|asin).
    std::string choice(std::initializer_list<const char*> choices) {
        const std::string v = str();
        for (const char* c : choices)
            if (v == c) return v;
        std::string all;
        for (const char* c : choices) all += std::string(all.empty() ? "" : "|") + c;
        throw std::runtime_error(args_[pos_ - 2] + ": expected " + all + ", got " + v);
    }

    [[noreturn]] void unknown(const std::string& flag) const {
        throw std::runtime_error("unknown argument " + flag + " (see the usage at the top of the source, or README.md)");
    }

    static std::vector<std::string> split(const std::string& s, char sep) {
        std::vector<std::string> out;
        std::stringstream ss(s);
        std::string item;
        while (std::getline(ss, item, sep)) out.push_back(item);
        return out;
    }

private:
    template <class T, class F>
    T convert(F f) {
        const std::string flag = args_[pos_ - 1], v = str();
        try {
            return f(v);
        } catch (const std::exception&) {
            throw std::runtime_error(flag + ": cannot parse '" + v + "'");
        }
    }

    std::vector<std::string> args_;
    size_t pos_ = 0;
};

}  // namespace fsm
