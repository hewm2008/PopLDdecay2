#!/bin/sh
# Build PopLDdecay2. POSIX-sh compatible: works as `sh make.sh` or `bash make.sh`.
# Uses CMake when available; otherwise falls back to the bundled Makefile
# (plain g++/make). After a successful build the binary is copied into bin/
# and the intermediate build/ tree is removed.
set -eu
cd "$(dirname "$0")"

# Support KEY=VALUE arguments (e.g. bash make.sh HTSLIB_ROOT=/path).
for _arg in "$@"; do
  case "$_arg" in
    HTSLIB_ROOT=*|CXX=*|NATIVE=*|LTO=*|KEEP_BUILD=*)
      eval "export $_arg" ;;
    *) echo "WARNING: ignored unknown argument: $_arg" >&2 ;;
  esac
done

# ---- htslib discovery (order: env HTSLIB_ROOT -> common prefixes) ----
HTS_ROOT=""
if [ -n "${HTSLIB_ROOT:-}" ] && [ -d "${HTSLIB_ROOT}/include/htslib" ]; then
  HTS_ROOT="$HTSLIB_ROOT"
fi
if [ -z "$HTS_ROOT" ]; then
  for p in "$HOME/01.Software/samtools-1.23" "$HOME/01.Software/htslib-1.23" /usr/local /usr; do
    if [ -d "$p/include/htslib" ]; then
      HTS_ROOT="$p"
      break
    fi
  done
fi

CLEAN_BUILD=1
[ "${KEEP_BUILD:-0}" = "1" ] && CLEAN_BUILD=0

# ---- pick a C++17-capable compiler (g++ >= 8; newest wins) ----
# Order tried: env CXX -> $HOME/01.Software/gcc-*/bin/g++ -> PATH c++/g++.

# POSIX-compatible realpath: resolves symlinks without GNU readlink -f.
# Handles bare command names (resolves via PATH first), absolute paths,
# and relative symlink chains.
_realpath() {
  _rp_target="$1"
  # Resolve bare names via PATH (e.g. "g++" -> "/usr/bin/g++").
  case "$_rp_target" in
    /*) ;;
    *) _rp_target=$(command -v "$_rp_target" 2>/dev/null) || _rp_target="$(pwd)/$1" ;;
  esac
  _rp_count=0
  while [ -L "$_rp_target" ]; do
    _rp_dir=$(dirname "$_rp_target")
    _rp_link=$(readlink "$_rp_target")
    case "$_rp_link" in
      /*) _rp_target="$_rp_link" ;;
      *)  _rp_target="$_rp_dir/$_rp_link" ;;
    esac
    _rp_count=$((_rp_count + 1))
    [ "$_rp_count" -gt 20 ] && break
  done
  _rp_dir=$(dirname "$_rp_target")
  printf '%s/%s' "$(cd "$_rp_dir" 2>/dev/null && pwd -P)" "$(basename "$_rp_target")"
}

pick_cxx() {
  _cands=""
  [ -n "${CXX:-}" ] && _cands="$_cands $CXX"
  for g in "$HOME"/01.Software/gcc-*/bin/g++; do
    [ -x "$g" ] && _cands="$_cands $g"
  done
  for g in c++ g++; do
    command -v "$g" >/dev/null 2>&1 && _cands="$_cands $g"
  done
  _best=""; _bestmaj=0; _bestver=""
  for c in $_cands; do
    v=$("$c" -dumpfullversion 2>/dev/null || "$c" -dumpversion 2>/dev/null) || continue
    case "$v" in ''|*[!0-9.]*) continue ;; esac
    maj=${v%%.*}
    case "$maj" in ''|*[!0-9]*) continue ;; esac
    if [ "$maj" -ge "${CXX_MIN_MAJOR:-8}" ] && { [ "$_bestmaj" -eq 0 ] || [ "$maj" -gt "$_bestmaj" ]; }; then
      _best="$c"; _bestmaj="$maj"; _bestver="$v"
    fi
  done
  if [ -n "$_best" ]; then
    # Resolve symlinks so the compiler's own runtime libs (libisl/mpc/mpfr...)
    # in prefix-style gcc trees are located next to the real cc1plus.
    _best=$(_realpath "$_best")
    CXX=$_best
    export CXX
    _cc=$(dirname "$_best")/gcc
    [ -x "$_cc" ] && CC=$_cc && export CC
    _gxx_root=$(dirname "$(dirname "$_best")")
    _ldlp=""
    for d in "$_gxx_root/lib64" "$_gxx_root/lib"; do
      [ -d "$d" ] || continue
      _ldlp="${_ldlp:+$_ldlp:}$d"
    done
    if [ -n "$_ldlp" ]; then
      LD_LIBRARY_PATH="${_ldlp}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
      export LD_LIBRARY_PATH
      # macOS uses DYLD_LIBRARY_PATH for runtime lib resolution.
      case "$(uname -s 2>/dev/null)" in
        Darwin)
          DYLD_LIBRARY_PATH="${_ldlp}${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}"
          export DYLD_LIBRARY_PATH
          ;;
      esac
    fi
    echo "compiler: $CXX ($_bestver)"
  else
    echo "WARNING: no C++17-capable g++ (>= ${CXX_MIN_MAJOR:-8}) found." >&2
    echo "  this project needs a newer compiler; point CXX at one and retry, e.g.:" >&2
    echo "    export CXX=\$HOME/01.Software/gcc-12.1.0/bin/g++ && bash make.sh" >&2
  fi
}
pick_cxx

# ---- detect parallel job count ----
_get_jobs() {
  _n=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
  [ "$_n" -gt 8 ] && _n=8
  printf '%s' "$_n"
}

# ---- probe whether static linking is feasible ----
# Checks if libhts.a exists in the HTS prefix or common system paths.
# Returns 0 = static possible, 1 = not.
probe_static_libs() {
  _hts="${HTS_ROOT:-}"
  _search=""
  [ -n "$_hts" ] && _search="$_search $_hts/lib $_hts/lib64"
  case "$(uname -s 2>/dev/null)" in
    Darwin) _search="$_search /opt/homebrew/lib /usr/local/lib" ;;
  esac
  _search="$_search /usr/lib /usr/lib64"
  for _p in $_search; do
    [ -f "$_p/libhts.a" ] && return 0
  done
  return 1
}

build_with_cmake() {
  # Concatenated flags (no bash arrays -> POSIX sh compatible).
  BASE_CMAKE_ARGS="-DCMAKE_BUILD_TYPE=Release"
  if [ "${NATIVE:-0}" = "1" ]; then
    echo "native build enabled (-march=native; machine-specific, may shift 1-ULP"
    echo "  boundary rounding -> golden byte-comparison requires NATIVE=0)"
    BASE_CMAKE_ARGS="$BASE_CMAKE_ARGS -DENABLE_NATIVE=ON"
  else
    # Explicitly OFF so a previous NATIVE=1 cache cannot silently persist.
    BASE_CMAKE_ARGS="$BASE_CMAKE_ARGS -DENABLE_NATIVE=OFF"
  fi
  if [ "${LTO:-0}" = "1" ]; then
    echo "LTO build enabled (-flto; semantics-neutral, byte-identical output)"
    BASE_CMAKE_ARGS="$BASE_CMAKE_ARGS -DENABLE_LTO=ON"
  else
    # Explicitly OFF so a previous LTO=1 cache cannot silently persist.
    BASE_CMAKE_ARGS="$BASE_CMAKE_ARGS -DENABLE_LTO=OFF"
  fi
  if [ -n "$HTS_ROOT" ]; then
    echo "htslib detected at: $HTS_ROOT"
    BASE_CMAKE_ARGS="$BASE_CMAKE_ARGS -DHTSLIB_ROOT=$HTS_ROOT"
    # static libhts.a may need sibling libs shipped in the same prefix (e.g. libdeflate)
    LIBRARY_PATH="${HTS_ROOT}/lib${LIBRARY_PATH:+:$LIBRARY_PATH}"
    export LIBRARY_PATH
  else
    echo "htslib not found in common prefixes; relying on pkg-config/system paths."
  fi

  mkdir -p build
  # Reconfigure from scratch if a previous configure used another compiler,
  # otherwise the old CMakeCache would silently keep building with it.
  if [ -f build/CMakeCache.txt ]; then
    cached=$(sed -n 's/^CMAKE_CXX_COMPILER:FILEPATH=//p' build/CMakeCache.txt | head -n 1)
    cur=$(command -v "${CXX:-g++}" 2>/dev/null || printf '%s' "${CXX:-g++}")
    if [ -n "$cached" ] && [ -n "$cur" ] && [ "$cached" != "$cur" ]; then
      echo "compiler changed ($cached -> $cur); wiping stale CMake cache"
      rm -rf build
      mkdir -p build
    fi
  fi

  # Try static linking first if feasible.
  CMAKE_ARGS="$BASE_CMAKE_ARGS"
  _static_log=""
  if probe_static_libs; then
    echo "attempting static build..."
    _static_log=$(mktemp)
    if cmake -S . -B build $CMAKE_ARGS >"$_static_log" 2>&1 && \
       cmake --build build -j"$(_get_jobs)" >"$_static_log" 2>&1; then
      rm -f "$_static_log"
      return 0
    fi
    rm -rf build && mkdir -p build
  fi

  # Fallback: dynamic linking.
  CMAKE_ARGS="$BASE_CMAKE_ARGS -DHTSLIB_PREFER_STATIC=OFF"
  echo "attempting dynamic build..."
  # shellcheck disable=SC2086 # intentional word splitting of flag string
  if cmake -S . -B build $CMAKE_ARGS && cmake --build build -j"$(_get_jobs)"; then
    rm -f "$_static_log"
    return 0
  fi

  # Both failed — show static error log if available.
  if [ -n "$_static_log" ] && [ -s "$_static_log" ]; then
    echo "--- static build log ---" >&2
    cat "$_static_log" >&2
  fi
  rm -f "$_static_log"
  return 1
}

build_with_make() {
  echo "cmake not found; falling back to the bundled Makefile (g++/make)."
  if [ "${NATIVE:-0}" = "1" ]; then
    echo "native build enabled (-march=native; machine-specific, may shift 1-ULP"
    echo "  boundary rounding -> golden byte-comparison requires NATIVE=0)"
    NATIVE=1
  else
    NATIVE=0
  fi
  if [ "${LTO:-0}" = "1" ]; then
    echo "LTO build enabled (-flto; semantics-neutral, byte-identical output)"
    LTO=1
  else
    LTO=0
  fi
  make clean >/dev/null 2>&1 || true

  # Try static linking first if feasible.
  _static_log=""
  if probe_static_libs; then
    echo "attempting static build..."
    _static_log=$(mktemp)
    if make -j"$(_get_jobs)" NATIVE="$NATIVE" LTO="$LTO" HTS="$HTS_ROOT" >"$_static_log" 2>&1; then
      rm -f "$_static_log"
      return 0
    fi
    make clean >/dev/null 2>&1 || true
  fi

  # Fallback: dynamic linking — force STATIC=0 so Makefile skips libhts.a.
  echo "attempting dynamic build..."
  if make -j"$(_get_jobs)" NATIVE="$NATIVE" LTO="$LTO" HTS="$HTS_ROOT" STATIC=0; then
    rm -f "$_static_log"
    return 0
  fi

  # Both failed — show static error log if available.
  if [ -n "$_static_log" ] && [ -s "$_static_log" ]; then
    echo "--- static build log ---" >&2
    cat "$_static_log" >&2
  fi
  rm -f "$_static_log"
  return 1
}

# cmake -S/-B needs >= 3.13; older cmake (e.g. 3.3 on CentOS 7) misparses it.
have_modern_cmake() {
  command -v cmake >/dev/null 2>&1 || return 1
  v=$(cmake --version | sed -n 's/^cmake version \([0-9][0-9]*\)\.\([0-9][0-9]*\).*/\1 \2/p')
  [ -n "$v" ] || return 1
  set -- $v
  [ "$1" -gt 3 ] && return 0
  { [ "$1" -eq 3 ] && [ "${2:-0}" -ge 13 ]; } && return 0
  return 1
}

if have_modern_cmake; then
  if ! build_with_cmake; then
    echo "WARNING: cmake configure/build failed; falling back to the bundled Makefile." >&2
    build_with_make
  fi
else
  build_with_make
fi

mkdir -p bin bin/mis
[ -f build/PopLDdecay2 ] || { echo "ERROR: build did not produce build/PopLDdecay2" >&2; exit 1; }
mv -f build/PopLDdecay2 bin/PopLDdecay2
echo "built: $(pwd)/bin/PopLDdecay2"

# Smoke-test BEFORE removing intermediates: only a binary that actually runs
# counts as a successful install. 126/127 = exec or loader failure.
./bin/PopLDdecay2 --help </dev/null >/dev/null 2>&1 && rc=0 || rc=$?
if [ "$rc" -eq 126 ] || [ "$rc" -eq 127 ]; then
  echo "ERROR: bin/PopLDdecay2 fails to execute (rc=$rc, missing/new lib?); keeping build/ for inspection." >&2
  exit 1
fi

if [ "$CLEAN_BUILD" = "1" ]; then
  rm -rf build
  echo "removed intermediate build/ (set KEEP_BUILD=1 to keep it)"
fi
