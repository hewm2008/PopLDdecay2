#pragma once

#include <string>
#include <vector>

namespace pld2 {

// Splits a line on spaces and tabs, dropping empty fields -- the original's
// split(line, inf, " \t"). With `cr`, a carriage return separates fields too,
// for the readers that take CRLF files.
//
// `toks` is reused: the vector keeps its size and the token strings keep
// their capacity, so a loop over a file does not allocate per line. It is
// resized to the token count.
inline void split_ws_into(const char* data, size_t len,
                          std::vector<std::string>& toks, bool cr = false) {
    size_t n = 0;
    size_t start = 0;
    for (size_t i = 0; i <= len; ++i) {
        const char c = (i < len) ? data[i] : '\0';
        if (i == len || c == ' ' || c == '\t' || (cr && c == '\r')) {
            if (i > start) {
                if (n < toks.size()) {
                    toks[n].assign(data + start, i - start);
                } else {
                    toks.emplace_back(data + start, i - start);
                }
                ++n;
            }
            start = i + 1;
        }
    }
    toks.resize(n);
}

inline void split_ws_into(const std::string& line,
                          std::vector<std::string>& toks, bool cr = false) {
    split_ws_into(line.data(), line.size(), toks, cr);
}

// The same split, into a fresh vector.
inline std::vector<std::string> split_ws(const std::string& line,
                                         bool cr = false) {
    std::vector<std::string> toks;
    split_ws_into(line, toks, cr);
    return toks;
}

// Single-quotes a string for use in a system() shell command (embedded single
// quotes become the POSIX '\'' splice). Handles spaces and shell metacharacters
// in user-supplied output prefixes.
inline std::string shell_quote(const std::string& s) {
    std::string out = "'";
    for (char c : s) {
        if (c == '\'') {
            out += "'\\''";
        } else {
            out.push_back(c);
        }
    }
    out += "'";
    return out;
}

}  // namespace pld2
