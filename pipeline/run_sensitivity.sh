#!/usr/bin/env bash
# One-at-a-time sensitivity runs around the base case. Each run is a full set of 16 slices (ITERS=2: about 10-12 min each on a laptop).
# usage: bash pipeline/run_sensitivity.sh            (all runs, ~1.5 h)   |   RUNS="base car55" bash pipeline/run_sensitivity.sh
# Results: data/work/sens/<run>/slice_stats.csv, then: python3 pipeline/summarize_sensitivity.py data/work/sens
set -euo pipefail
W=${W:-data/work}; S=$W/sens; mkdir -p "$S"; export ITERS=${ITERS:-2}
RUNS=${RUNS:-"base car35 car55 early30 mix cap80 cap120 beta_lo beta_hi"}
for r in $RUNS; do
  case $r in
    base)    env_="" ;;
    car35)   env_="CAR_WORK=0.35" ;;
    car55)   env_="CAR_WORK=0.55" ;;
    early30) env_="WORK_PROFILE=early30" ;;
    mix)     env_="WORK_PROFILE=mix" ;;
    cap80)   env_="CAP_SCALE=0.8" ;;
    cap120)  env_="CAP_SCALE=1.2" ;;
    beta_lo) env_="BETA_WORK=0.05" ;;
    beta_hi) env_="BETA_WORK=0.12" ;;
    *) echo "unknown run $r"; exit 1 ;;
  esac
  [ -f "$S/$r/slice_stats.csv" ] && { echo "skip $r (done)"; continue; }
  echo "=== $r ($env_)"; env $env_ python3 pipeline/congestion_slices.py "$W/road_edges.csv.gz" "$W/mahalla_zones.csv" "$W/students/pa_vectors.npz" "$S/$r" > "$S/$r.log"
done
python3 pipeline/summarize_sensitivity.py "$S"
