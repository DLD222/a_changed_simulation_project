#!/usr/bin/env bash
# M2a-1b L2 criteria for the 2D XZ TE mode.
#
#   L2a  regular domain (EB outside of the box): ECT == Yee to round-off
#   L2b  grid-aligned PEC plane: ECT puts the wall at the true position, the
#        staircase Yee solver one cell inside
#
# Both criteria are described in ECT_2D_TE_REPORT.md (section "L2 criteria").
set -e
cd "$(dirname "$0")"
BIN=${WARPX_BIN:-build_2d/bin/warpx.2d.MPI.OMP.DP.PDP.EB}
RUNS=${RUNS:-/tmp/l2}
rm -rf "$RUNS"
mkdir -p "$RUNS"

run_case () {  # $1 = case name, $2 = input file
  mkdir -p "$RUNS/$1"
  cp "$2" "$RUNS/$1/inputs"
  ( cd "$RUNS/$1" && OMP_NUM_THREADS=${OMP_NUM_THREADS:-8} "$BIN" inputs > run.log 2>&1 )
}

echo "=== running the four cases (this takes a few minutes) ==="
run_case ect_noeb inputs_ect_2d_noeb
run_case yee_noeb inputs_yee_2d_noeb
run_case ect_pec  inputs_ect_2d_degen
run_case yee_pec  inputs_yee_2d_degen

echo
echo "=== L2a: ECT vs Yee on a regular domain (expect bitwise agreement) ==="
python3 verify_l2a.py "$RUNS/ect_noeb" "$RUNS/yee_noeb"

echo
echo "=== L2b: grid-aligned PEC plane (expect ECT wall at z=0, Yee one cell inside) ==="
python3 verify_l2b.py "$RUNS/ect_pec" "$RUNS/yee_pec"
