# Minimal Makefile fallback for systems without CMake.
# Mirrors CMakeLists.txt's Release build: -O3 -std=c++17 -Wall -Wextra,
# automatic htslib detection, optional NATIVE=1 (-march=native) and LTO=1 (-flto).
# Used by make.sh only when `cmake` is not installed.
CXX      ?= g++
SRCS     := $(wildcard src/*.cpp)
OBJS     := $(patsubst src/%.cpp,build/%.o,$(SRCS))
HTS      ?=
STATIC   ?= 1

CXXFLAGS  = -O3 -DNDEBUG -std=c++17 -Wall -Wextra

ifeq ($(NATIVE),1)
CXXFLAGS += -march=native -mtune=native
endif
ifeq ($(LTO),1)
CXXFLAGS += -flto
endif

# ---- htslib detection (order: HTSLIB_ROOT -> common prefixes) ----
ifeq ($(strip $(HTS)),)
HTS := $(shell sh -c 'if [ -n "$${HTSLIB_ROOT:-}" ] && [ -d "$${HTSLIB_ROOT}/include/htslib" ]; then echo "$${HTSLIB_ROOT}"; \
  elif [ -d "$$HOME/01.Software/samtools-1.23/include/htslib" ]; then echo "$$HOME/01.Software/samtools-1.23"; \
  elif [ -d "$$HOME/01.Software/htslib-1.23/include/htslib" ]; then echo "$$HOME/01.Software/htslib-1.23"; \
  elif [ -d /usr/local/include/htslib ]; then echo /usr/local; \
  elif [ -d /usr/include/htslib ]; then echo /usr; \
  else echo ""; fi')
endif

ifeq ($(strip $(HTS)),)
$(error htslib not found. Set HTSLIB_ROOT=<prefix> (must contain include/htslib and lib/libhts) or edit HTS=)
endif

CXXFLAGS += -I$(HTS)/include -Isrc
# Header dependencies: without these, editing a header (types.h, cli.h, ...)
# leaves the stale objects in build/ in place and the "make" build silently
# disagrees with the source. CMake tracks them by default; this Makefile is
# the fallback used when cmake is missing, so it has to do it itself.
CXXFLAGS += -MMD -MP

# ---------------------------------------------------------------------------
# The kernel counts haplotypes with __builtin_popcountll: a strict x86-64
# baseline (and clang's default) emit software code for it, and the kernel
# measured 3.4x slower that way. popcount is exact integer arithmetic, so the
# hardware instruction changes the speed and not one output byte.
# ENABLE_MPOPCNT=0 turns it off; AArch64 needs nothing (cnt is a base
# instruction).
accepts = printf 'int main(){return 0;}' | $(CXX) -x c++ $(1) - -o /dev/null >/dev/null 2>&1 && echo YES
MPOPCNT_ARCH ?= $(shell uname -m)
ifeq ($(filter 0 OFF off,$(ENABLE_MPOPCNT)),)
ifneq ($(filter x86_64 amd64 i%86,$(MPOPCNT_ARCH)),)
ifeq ($(shell $(call accepts,-mpopcnt)),YES)
CXXFLAGS += -mpopcnt
endif
endif
endif

# FP contraction (a*b+c -> fma) changes the rounding and depends on the
# compiler's default target; turning it off makes the output independent of it.
ifeq ($(filter 0 OFF off,$(ENABLE_FP_CONTRACT_OFF)),)
ifeq ($(shell $(call accepts,-ffp-contract=off)),YES)
CXXFLAGS += -ffp-contract=off
endif
endif


# Prefer static libhts.a; fall back to shared libhts.so.
# Set STATIC=0 to force dynamic linking (used by make.sh fallback).
ifeq ($(STATIC),1)
ifneq ($(wildcard $(HTS)/lib/libhts.a),)
HTSLIB   := $(HTS)/lib/libhts.a
HTSDEPS  := -lm -lbz2 -llzma
# htslib may be built against libdeflate; probe HTS prefix then system paths.
HTS_DEF := $(or $(wildcard $(HTS)/lib/libdeflate.a), \
               $(wildcard $(HTS)/lib64/libdeflate.a), \
               $(wildcard $(shell echo /opt/homebrew/lib/libdeflate.a /usr/local/lib/libdeflate.a /usr/lib/libdeflate.a /usr/lib64/libdeflate.a 2>/dev/null)))
ifneq ($(HTS_DEF),)
HTSDEPS += $(HTS_DEF)
else
  HTS_DEF := $(or $(wildcard $(HTS)/lib/libdeflate.so), \
                 $(wildcard $(HTS)/lib64/libdeflate.so), \
                 $(wildcard $(shell echo /opt/homebrew/lib/libdeflate.so /usr/local/lib/libdeflate.so /usr/lib/libdeflate.so /usr/lib64/libdeflate.so /opt/homebrew/lib/libdeflate.dylib /usr/local/lib/libdeflate.dylib /usr/lib/libdeflate.dylib 2>/dev/null)))
  ifneq ($(HTS_DEF),)
    HTSDEPS += $(HTS_DEF)
  endif
endif
# A static htslib built with S3 support references libcurl, and its TLS
# backend may reference libcrypto. Probe for them the way the cmake build
# does and add only what is present: naming a library that is not installed
# breaks the link, and omitting one that htslib needs breaks it too. Without
# this the Makefile fallback could not link against the very htslib the cmake
# build links against.
HTS_EXTRA := $(shell sh -c 'for dep in curl crypto; do \
    for d in "$(HTS)/lib" "$(HTS)/lib64" /usr/local/lib /usr/lib64 /usr/lib; do \
      if [ -f "$$d/lib$$dep.a" ]; then echo "$$d/lib$$dep.a"; continue 2; fi; \
    done; \
    for d in /usr/local/lib /usr/lib64 /usr/lib; do \
      if ls "$$d"/lib$$dep.so* >/dev/null 2>&1; then echo "-L$$d -l$$dep"; continue 2; fi; \
    done; \
  done')
ifneq ($(HTS_EXTRA),)
  HTSDEPS += $(HTS_EXTRA)
endif
endif
endif
ifeq ($(HTSLIB),)
HTSLIB   := -L$(HTS)/lib -lhts
HTSDEPS  :=
endif

# Apple Clang does not ship libgcc; skip -static-libgcc on macOS.
# Only add -static-libstdc++ / -static-libgcc if supported by the toolchain.
_STATIC_FLAGS :=
_STATIC_TEST=$(shell echo 'int main(){return 0;}' | $(CXX) -x c++ -static-libstdc++ - -o /dev/null 2>/dev/null && echo YES)
ifeq ($(_STATIC_TEST),YES)
  _STATIC_FLAGS += -static-libstdc++
endif
ifeq ($(shell uname -s 2>/dev/null),Darwin)
  # Apple Clang has no libgcc
else
  _STATIC_TEST_GCC=$(shell echo 'int main(){return 0;}' | $(CXX) -x c++ -static-libgcc - -o /dev/null 2>/dev/null && echo YES)
  ifeq ($(_STATIC_TEST_GCC),YES)
    _STATIC_FLAGS += -static-libgcc
  endif
endif
LIBS     := $(HTSLIB) $(HTSDEPS) -lz -lpthread $(_STATIC_FLAGS)

all: build/PopLDdecay2

build/PopLDdecay2: $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

build/%.o: src/%.cpp | build
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Generated by -MMD -MP next to each object.
DEPS := $(OBJS:.o=.d)
-include $(DEPS)

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: all clean