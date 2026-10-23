#pragma once

#include <string>

#include "types.h"

namespace pld2 {

// Parses command line; returns true on success (may set options).
// On failure / help request, prints message and returns false.
bool parse_cli(int argc, char** argv, Options& opt);

// Strips a trailing ".gz" and, if the basename != "stat", a trailing ".stat"
// (byte-for-byte replication of original LD_Decay.cpp prefix handling).
std::string strip_stat_prefix(const std::string& out_stat);

// Strict numeric parsing: the whole token must be a number, because
// atoi()/atof() cannot tell "abc" from 0 (-MAF abc -> 0, the filter off).
bool parse_int_strict(const char* s, long long& out);

}  // namespace pld2
