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

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: all clean