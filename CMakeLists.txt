# 3.5 is the lowest version CMake 4 still accepts (CMake 4.0 removed the
# compatibility with anything older, and 3.27 warns about it). The range
# syntax (3.12+) sets the policy version to 3.28, so a current cmake does
# not warn about the policies introduced in between; older cmake ignores the
# "...3.28" part and behaves as 3.5.
cmake_minimum_required(VERSION 3.5...3.28)
project(PopLDdecay2 CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# CMAKE_CXX_STANDARD 17 is only honoured from CMake 3.8 on, and the code does
# need C++17 (if constexpr in ldkernel.cpp, <string_view> in region.h). On an
# older CMake the standard property above is silently ignored and the compiler
# falls back to its default (gnu++14 on GCC 7/8), which fails to compile. Pass
# the flag directly in that case.
if(CMAKE_VERSION VERSION_LESS 3.8)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

include(CheckCXXSourceCompiles)

# ---------------------------------------------------------------------------
# Build flags that decide how fast the LD kernel runs. Neither of them changes
# a single output byte; both are checked so that an old or unusual compiler is
# simply left alone.
# ---------------------------------------------------------------------------

# The LD kernel counts haplotypes with __builtin_popcountll. Whether that is
# the hardware POPCNT instruction or a software sequence depends on the
# compiler's default target, not on ours: a strict x86-64 baseline (and
# clang's default) emit software code, and the kernel measured 3.4x slower
# that way. popcount is exact integer arithmetic, so the instruction changes
# the speed and not one output byte. AArch64 needs nothing: cnt is a base
# instruction.
option(ENABLE_MPOPCNT "Use the hardware POPCNT instruction on x86" ON)
if(ENABLE_MPOPCNT AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86|amd64|AMD64|i[3-6]86")
  set(CMAKE_REQUIRED_FLAGS "-mpopcnt")
  check_cxx_source_compiles("int main(){return __builtin_popcountll(1)==1?0:1;}" HAVE_MPOPCNT)
  unset(CMAKE_REQUIRED_FLAGS)
  if(HAVE_MPOPCNT)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mpopcnt")
  else()
    message(STATUS "no -mpopcnt; the kernel will use the software popcount")
  endif()
endif()

# FP contraction (a*b+c -> fma) changes the rounding, and whether the compiler
# does it depends on its default target too (x86-64-v3, or -march=native).
# Turning it off makes the output independent of that.
option(ENABLE_FP_CONTRACT_OFF "Disable FP contraction" ON)
if(ENABLE_FP_CONTRACT_OFF)
  set(CMAKE_REQUIRED_FLAGS "-ffp-contract=off")
  check_cxx_source_compiles("int main(){return 0;}" HAVE_FFP_CONTRACT_OFF)
  unset(CMAKE_REQUIRED_FLAGS)
  if(HAVE_FFP_CONTRACT_OFF)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffp-contract=off")
  endif()
endif()

# Optional native-tuned build (-march=native): enables ISA-specific
# optimizations for the build machine. OFF by default for portability;
# enable via `NATIVE=1 ./make.sh` (or -DENABLE_NATIVE=ON).
option(ENABLE_NATIVE "Build with -march=native (tuned to the build machine)" OFF)
if(ENABLE_NATIVE)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native")
  message(STATUS "ENABLE_NATIVE: -march=native enabled (binary is machine-specific)")
endif()

# Optional link-time optimization (-flto): lets cross-TU inlining hit the LD
# hot loops in ldkernel.cpp. OFF by default for toolchain portability; enable
# via `LTO=1 ./make.sh` (or -DENABLE_LTO=ON). Semantics-neutral (strict IEEE
# FP is preserved), so output stays byte-identical.
option(ENABLE_LTO "Build with -flto (link-time optimization)" OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ---------------- htslib detection ----------------
# Order: HTSLIB_ROOT -> pkg-config -> common prefixes/system default.
# Prefer static libhts.a; fall back to shared libhts.so (with rpath).
set(HTSLIB_ROOT "" CACHE PATH "htslib installation prefix (contains include/htslib and lib/libhts*)")
set(HTSLIB_STATIC_DEPS "" CACHE STRING "extra link libs needed by a static libhts.a (e.g. 'm bz2 lzma deflate curl crypto')")
option(HTSLIB_PREFER_STATIC "Prefer static libhts.a when both static and shared are available" ON)

set(HTSLIB_FOUND FALSE)
set(HTSLIB_INCLUDE_DIRS "")
set(HTSLIB_LINK_LIBRARIES "")
set(HTSLIB_LINK_DIRS "")

if(HTSLIB_ROOT)
  set(HTSLIB_INCLUDE_DIRS "${HTSLIB_ROOT}/include")
  if(EXISTS "${HTSLIB_ROOT}/lib/libhts.a")
    set(HTSLIB_LINK_LIBRARIES "${HTSLIB_ROOT}/lib/libhts.a")
    set(HTSLIB_FOUND TRUE)
    set(HTSLIB_IS_STATIC TRUE)
  elseif(EXISTS "${HTSLIB_ROOT}/lib/libhts.so")
    set(HTSLIB_LINK_LIBRARIES "${HTSLIB_ROOT}/lib/libhts.so")
    set(HTSLIB_LINK_DIRS "${HTSLIB_ROOT}/lib")
    set(HTSLIB_FOUND TRUE)
  endif()
endif()

if(NOT HTSLIB_FOUND)
  find_package(PkgConfig QUIET)
  if(PKG_CONFIG_FOUND)
    pkg_check_modules(PC_HTSLIB QUIET htslib)
  endif()
  if(PC_HTSLIB_FOUND)
    set(HTSLIB_FOUND TRUE)
    set(HTSLIB_INCLUDE_DIRS ${PC_HTSLIB_INCLUDE_DIRS})
    set(HTSLIB_LINK_LIBRARIES ${PC_HTSLIB_LINK_LIBRARIES})
  endif()
endif()

if(NOT HTSLIB_FOUND)
  foreach(p "$ENV{HOME}/01.Software/samtools-1.23" "$ENV{HOME}/01.Software/htslib-1.23" /usr/local /usr)
    if(EXISTS "${p}/include/htslib/hts.h")
      if(EXISTS "${p}/lib/libhts.a")
        set(HTSLIB_INCLUDE_DIRS "${p}/include")
        set(HTSLIB_LINK_LIBRARIES "${p}/lib/libhts.a")
        set(HTSLIB_FOUND TRUE)
        set(HTSLIB_IS_STATIC TRUE)
        break()
      elseif(EXISTS "${p}/lib/libhts.so")
        set(HTSLIB_INCLUDE_DIRS "${p}/include")
        set(HTSLIB_LINK_LIBRARIES "${p}/lib/libhts.so")
        set(HTSLIB_LINK_DIRS "${p}/lib")
        set(HTSLIB_FOUND TRUE)
        break()
      endif()
    endif()
  endforeach()
endif()

if(NOT HTSLIB_FOUND)
  message(FATAL_ERROR
    "htslib not found.\n"
    "  Point to it, e.g.:\n"
    "    cmake -S . -B build -DHTSLIB_ROOT=<your-htslib-prefix>\n"
    "  or set env HTSLIB_ROOT to your htslib prefix (must contain include/htslib and lib/libhts).\n"
    "  htslib is typically installed via: samtools/htslib source build, Homebrew, apt, or conda.")
endif()

# htslib version. The check reads the HTS_VERSION macro from the header, so it
# neither links nor runs anything.
#   < 1.16  the fast path needs tbx_seqnames(), which does not exist yet
#   < 1.18  it works, but 1.18 and newer are noticeably faster
set(CMAKE_REQUIRED_INCLUDES ${HTSLIB_INCLUDE_DIRS})
check_cxx_source_compiles("#include <htslib/hts.h>
#if !defined(HTS_VERSION) || HTS_VERSION < 101600
#error htslib older than 1.16
#endif
int main(){return 0;}" HTSLIB_GE_116)
check_cxx_source_compiles("#include <htslib/hts.h>
#if !defined(HTS_VERSION) || HTS_VERSION < 101800
#error htslib older than 1.18
#endif
int main(){return 0;}" HTSLIB_GE_118)
unset(CMAKE_REQUIRED_INCLUDES)
if(NOT HTSLIB_GE_116)
  message(FATAL_ERROR
    "htslib 1.16 or newer is required (the fast path uses tbx_seqnames()).\n"
    "  Found an older one; install a newer htslib, or point the build at one with\n"
    "    cmake -S . -B build -DHTSLIB_ROOT=<your-htslib-prefix>")
elseif(NOT HTSLIB_GE_118)
  message(STATUS "htslib is older than 1.18: it works, but 1.18 and newer are faster")
endif()

message(STATUS "htslib include: ${HTSLIB_INCLUDE_DIRS}")
message(STATUS "htslib lib    : ${HTSLIB_LINK_LIBRARIES}")

find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

# Bundle the C++ runtime so the binary also runs on hosts whose system
# libstdc++ is older than the build compiler's (devtoolset-style trees).
# Only add -static-libstdc++ if the static library is available.
set(CMAKE_REQUIRED_FLAGS "-static-libstdc++")
check_cxx_source_compiles("int main(){return 0;}" HAVE_STATIC_LIBSTDCXX)
unset(CMAKE_REQUIRED_FLAGS)
if(HAVE_STATIC_LIBSTDCXX)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++")
endif()
# Apple Clang does not ship libgcc; skip -static-libgcc on macOS.
if(NOT APPLE)
  set(CMAKE_REQUIRED_FLAGS "-static-libgcc")
  check_cxx_source_compiles("int main(){return 0;}" HAVE_STATIC_LIBGCC)
  unset(CMAKE_REQUIRED_FLAGS)
  if(HAVE_STATIC_LIBGCC)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc")
  endif()
endif()

add_executable(PopLDdecay2
  src/main.cpp
  src/cli.cpp
  src/region.cpp
  src/subpop.cpp
  src/vcf_reader.cpp
  src/genotype_reader.cpp
  src/plink_reader.cpp
  src/ldkernel.cpp
  src/ldwriter.cpp
  src/output.cpp
  src/ehh.cpp
)
target_include_directories(PopLDdecay2 PRIVATE src ${HTSLIB_INCLUDE_DIRS})
target_compile_options(PopLDdecay2 PRIVATE -Wall -Wextra)

if(ENABLE_LTO)
  if(NOT CMAKE_VERSION VERSION_LESS 3.10)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT _ipo_supported)
    if(_ipo_supported)
      set_target_properties(PopLDdecay2 PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
      message(STATUS "ENABLE_LTO: -flto enabled")
    else()
      message(WARNING "ENABLE_LTO requested but the toolchain does not support IPO; ignored")
    endif()
  else()
    message(WARNING "ENABLE_LTO needs cmake >= 3.10; ignored on ${CMAKE_VERSION}")
  endif()
endif()

if(HTSLIB_IS_STATIC AND HTSLIB_PREFER_STATIC)
  if(NOT HTSLIB_STATIC_DEPS)
    # Probe the libraries htslib builds commonly reference and link only the
    # ones actually present (by full path, so no -L juggling and missing
    # optional deps are skipped instead of breaking the link step).
    # Also search common system paths (Homebrew /opt/homebrew, /usr/local)
    # for deps that live outside the HTS prefix (e.g. libdeflate on macOS).
    set(_HTS_SYS_LIB_DIRS "")
    if(HTSLIB_ROOT)
      list(APPEND _HTS_SYS_LIB_DIRS "${HTSLIB_ROOT}/lib" "${HTSLIB_ROOT}/lib64")
    endif()
    if(APPLE)
      list(APPEND _HTS_SYS_LIB_DIRS "/opt/homebrew/lib" "/usr/local/lib")
    endif()
    list(APPEND _HTS_SYS_LIB_DIRS "/usr/lib" "/usr/lib64")
    foreach(_dep m bz2 lzma deflate curl crypto)
      find_library(_HTS_DEP_${_dep} NAMES "lib${_dep}.a" ${_dep}
                   HINTS ${_HTS_SYS_LIB_DIRS})
      if(_HTS_DEP_${_dep})
        list(APPEND HTSLIB_STATIC_DEPS "${_HTS_DEP_${_dep}}")
      endif()
      unset(_HTS_DEP_${_dep} CACHE)
    endforeach()
  endif()
  message(STATUS "htslib static deps: ${HTSLIB_STATIC_DEPS}")
  target_link_libraries(PopLDdecay2 PRIVATE
    ${HTSLIB_LINK_LIBRARIES}
    ${HTSLIB_STATIC_DEPS}
    ZLIB::ZLIB
    Threads::Threads)
else()
  # Dynamic linking: ensure the linker can find libhts.so/dylib.
  if(HTSLIB_ROOT AND NOT HTSLIB_LINK_DIRS)
    set(HTSLIB_LINK_DIRS "${HTSLIB_ROOT}/lib")
  endif()
  target_link_libraries(PopLDdecay2 PRIVATE
    ${HTSLIB_LINK_LIBRARIES}
    ZLIB::ZLIB
    Threads::Threads)
  if(HTSLIB_LINK_DIRS)
    set_target_properties(PopLDdecay2 PROPERTIES
      INSTALL_RPATH "${HTSLIB_LINK_DIRS}"
      BUILD_RPATH "${HTSLIB_LINK_DIRS}")
  endif()
endif()

# ---------------- install ----------------
# The R plotting scripts are consumed standalone (see the "Used [perl
# ../bin/Plot_XX.pl ]" message): install them alongside the binary so the
# message matches a real location. `bin/mis/plink2genotype.pl` (PED/MAP ->
# genotype) emits the native -InGenotype format (one IUPAC code per sample).
install(TARGETS PopLDdecay2 RUNTIME DESTINATION bin)
install(PROGRAMS
  bin/Plot_OnePop.pl
  bin/Plot_MultiPop.pl
  DESTINATION bin)
install(FILES bin/mis/plink2genotype.pl
  DESTINATION bin/mis)

# ---------------- tests ----------------
# A standard entry point for the regression gate: `ctest` in the build
# directory runs the same script CI runs. The golden comparison against
# PopLDdecay 3.45 is NOT wired in here: 3.45 is not redistributable, so it is
# only present on an author's machine (RunTest/06.Verify/golden_diff.sh).
enable_testing()
add_test(verify bash ${CMAKE_CURRENT_SOURCE_DIR}/RunTest/06.Verify/step1_checkOutput.sh)
# The script defaults to bin/PopLDdecay2 (where make.sh puts the binary); in a
# cmake build tree the binary under test is the one just built, so point the
# test at it explicitly.
set_tests_properties(verify PROPERTIES
                     ENVIRONMENT "NEW=${CMAKE_CURRENT_BINARY_DIR}/PopLDdecay2")