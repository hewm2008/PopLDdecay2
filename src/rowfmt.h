#pragma once

#include <cstdarg>
#include <cstdio>
#include <string>
#include <vector>

#include <zlib.h>

namespace pld2 {

// One formatted text row, written without a fixed upper bound on its length.
//
// snprintf() returns the length the row WOULD have had with an unbounded
// buffer. When the row does not fit, that return value is LARGER than the
// buffer it was written into, so handing it on as a write length reads past
// the end of the buffer. These helpers keep the caller's stack buffer for the
// common case and re-format the row into a right-sized heap buffer when it
// does not fit, so the bytes written are always exactly the formatted row.
//
// The original PopLDdecay wrote these rows through a C++ stream, which has no
// length limit, so writing the full row is also what byte-identity with it
// requires: a chromosome name longer than the stack buffer used to make this
// tool read out of bounds and lose the tail of every row.

// Writes one formatted row to a gzip stream. Returns false if the row could
// not be formatted (encoding error) or not written.
inline bool write_row(gzFile gz, char* stackbuf, size_t cap, const char* fmt,
                      ...) {
    va_list ap;
    va_start(ap, fmt);
    va_list ap2;
    va_copy(ap2, ap);
    const int len = std::vsnprintf(stackbuf, cap, fmt, ap);
    va_end(ap);
    if (len < 0) {
        va_end(ap2);
        return false;
    }
    bool ok;
    if (static_cast<size_t>(len) < cap) {
        ok = gzwrite(gz, stackbuf, static_cast<unsigned>(len)) ==
             static_cast<int>(len);
    } else {
        std::vector<char> heap(static_cast<size_t>(len) + 1);
        std::vsnprintf(heap.data(), heap.size(), fmt, ap2);
        ok = gzwrite(gz, heap.data(), static_cast<unsigned>(len)) ==
             static_cast<int>(len);
    }
    va_end(ap2);
    return ok;
}

// Same, appending the row to a string instead of writing it to a stream.
inline bool append_row(std::string& out, char* stackbuf, size_t cap,
                       const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    va_list ap2;
    va_copy(ap2, ap);
    const int len = std::vsnprintf(stackbuf, cap, fmt, ap);
    va_end(ap);
    if (len < 0) {
        va_end(ap2);
        return false;
    }
    if (static_cast<size_t>(len) < cap) {
        out.append(stackbuf, static_cast<size_t>(len));
    } else {
        std::vector<char> heap(static_cast<size_t>(len) + 1);
        std::vsnprintf(heap.data(), heap.size(), fmt, ap2);
        out.append(heap.data(), static_cast<size_t>(len));
    }
    va_end(ap2);
    return true;
}

}  // namespace pld2
