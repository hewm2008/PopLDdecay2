cmake_minimum_required(VERSION 2.8.12)
project(PopLDdecay2 CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
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

message(STATUS "htslib include: ${HTSLIB_INCLUDE_DIRS}")
message(STATUS "htslib lib    : ${HTSLIB_LINK_LIBRARIES}")

find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

# Bundle the C++ runtime so the binary also runs on hosts whose system
# libstdc++ is older than the build compiler's (devtoolset-style trees).
# Only add -static-libstdc++ if the static library is available.
include(CheckCXXSourceCompiles)
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
install(FILES scripts/plink2genotype.pl
  DESTINATION bin/mis)