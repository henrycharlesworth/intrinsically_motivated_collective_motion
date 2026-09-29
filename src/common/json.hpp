// A small JSON writer (ordered keys, one top-level key per line, doubles printed round-trip exact) and a reader for
// the flat numeric fields of our own run files.
#pragma once

#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

namespace fsm {

class Json {
public:
    Json& num(const std::string& key, double v) { return raw(key, number(v)); }
    Json& num(const std::string& key, long long v) { return raw(key, std::to_string(v)); }
    Json& num(const std::string& key, unsigned long long v) { return raw(key, std::to_string(v)); }
    Json& num(const std::string& key, int v) { return raw(key, std::to_string(v)); }
    Json& num(const std::string& key, size_t v) { return raw(key, std::to_string(v)); }
    Json& str(const std::string& key, const std::string& v) { return raw(key, quote(v)); }
    Json& boolean(const std::string& key, bool v) { return raw(key, v ? "true" : "false"); }
    Json& nums(const std::string& key, const std::vector<double>& v) {
        std::string s = "[";
        for (size_t i = 0; i < v.size(); i++) s += (i ? ", " : "") + number(v[i]);
        return raw(key, s + "]");
    }
    Json& obj(const std::string& key, const Json& v) { return raw(key, v.inline_text()); }

    // Multi-line text of a top-level object (one key per line).
    std::string text() const {
        std::string s = "{\n";
        for (size_t i = 0; i < items_.size(); i++) s += "  " + items_[i] + (i + 1 < items_.size() ? ",\n" : "\n");
        return s + "}\n";
    }
    std::string inline_text() const {
        std::string s = "{";
        for (size_t i = 0; i < items_.size(); i++) s += (i ? ", " : "") + items_[i];
        return s + "}";
    }

    static std::string number(double v) {
        std::ostringstream o;
        o.precision(17);
        o << v;
        return o.str();
    }
    static std::string quote(const std::string& v) {
        std::string s = "\"";
        for (char c : v) s += (c == '"' || c == '\\') ? std::string("\\") + c : std::string(1, c);
        return s + "\"";
    }

private:
    Json& raw(const std::string& key, const std::string& value) {
        items_.push_back(quote(key) + ": " + value);
        return *this;
    }
    std::vector<std::string> items_;
};

// Value of the first numeric field `key` in a JSON text (def if absent). Enough for our flat run files.
inline double json_number(const std::string& js, const std::string& key, double def) {
    auto q = js.find("\"" + key + "\"");
    if (q == std::string::npos) return def;
    q = js.find(':', q);
    return std::atof(js.c_str() + q + 1);
}

}  // namespace fsm
