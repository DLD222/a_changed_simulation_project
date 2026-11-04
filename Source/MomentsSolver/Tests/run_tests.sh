#!/usr/bin/env sh
# M1 unit tests: compile against AMReX headers and run.
set -e
cd "$(dirname "$0")/../../.."   # -> ~/WarpX_fluid

INC="-ISource -ISource/MomentsSolver"
CXX="g++ -std=c++17 -O2"

for t in euler_sod ten_moment maxwell langmuir em_dispersion ten_moment_k0 ten_moment_wave friction eos fv_nd; do
  echo "=== $t ==="
  $CXX $INC "Source/MomentsSolver/Tests/test_${t}.cpp" -o "/tmp/test_${t}"
  "/tmp/test_${t}"
done
