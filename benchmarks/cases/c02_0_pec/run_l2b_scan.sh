#!/usr/bin/env bash
# L2b resolution scan: grid-aligned PEC plane, ECT vs Yee, n = 64/128/256.
#
# For every resolution the analytic round-trip time 2*(z_wall - z_probe)/c is
# compared with the measured t_ref - t_inc. The round trip (not the absolute
# arrival time) is the robust quantity: the probe reads a carrier peak whose
# phase drifts with dz, and both t_inc and t_ref drift together.
#
# Usage:  ./run_l2b_scan.sh [n1 n2 ...]        (default: 64 128 256)
# Env:    WARPX_BIN, RUNS (default /tmp/l2scan), OMP_NUM_THREADS
set -e
cd "$(dirname "$0")"
BIN=${WARPX_BIN:-build_2d/bin/warpx.2d.MPI.OMP.DP.PDP.EB}
RUNS=${RUNS:-/tmp/l2scan}
NS=${*:-64 128 256}
mkdir -p "$RUNS"

for n in $NS; do
  # keep the physical time fixed: dt scales like 1/n, so does the step count
  ms=$((200 * n / 64))
  for v in ect yee; do
    d="$RUNS/${v}${n}"
    mkdir -p "$d"
    sed -e "s/^amr.n_cell = 64 64/amr.n_cell = $n $n/" \
        -e "s/^max_step = 200/max_step = $ms/" \
        "inputs_${v}_2d_degen" > "$d/inputs"
    ( cd "$d" && OMP_NUM_THREADS=${OMP_NUM_THREADS:-8} "$BIN" inputs > run.log 2>&1 )
  done
  echo "=== n = $n (max_step = $ms) ==="
  python3 verify_l2b.py "$RUNS/ect${n}" "$RUNS/yee${n}"
  echo
done
