#!/usr/bin/env bash

set -euo pipefail

# Two-phase script for 3D HP benchmark sequences:
# Phase 1: Optimize settings for each sequence (find best L, seed, params)
# Phase 2: Run solver with optimized settings and longer time limit

# ============ Configuration ============
# Phase 1: Optimization settings
OPT_TIME_BUDGET=${OPT_TIME_BUDGET:-300}     # 5 min optimization per sequence
OPT_SEEDS=${OPT_SEEDS:-5}                    # Seeds to try during optimization

# Phase 2: Final run settings
RUN_TIME_LIMIT=${RUN_TIME_LIMIT:-300}        # 5 min solve time for final run
WORKERS=${WORKERS:-8}

# Directories
OPT_DIR="out/optimized"
RESULTS_DIR="out/results"
mkdir -p "$OPT_DIR" "$RESULTS_DIR"

# ============ Benchmark Sequences ============
# Format: "id:length:sequence:best_known"
cases=(
  "3d1:20:(HP)2PH(HP)2(PH)2HP(PH)2:-11"
  "3d2:24:H2P2(HP2)6H2:-13"
  "3d3:25:P2HP2(H2P4)3H2:-9"
  "3d4:36:P(P2H2)2P5H5(H2P2)2P2H(HP2)2:-18"
  "3d5:46:P2H3PH3P3HPH2PH2P2HPH4PHP2H5PHPH2P2H2P:-35"
  "3d6:48:P2H(P2H2)2P5H10P6(H2P2)2HP2H5:-31"
  "3d7:50:H2(PH)3PH4PH(P3H)2P4(HP3)2HPH4(PH)3PH2:-34"
  "3d8:58:PH(PH3)2P(PH2PH)2H(HP)3(H2P2H)2PHP4(H(P2H)2)2:-44"
  "3d9:60:P(PH3)3H5P3H10PHP3H12P4H6PH2PH:-55"
)

# ============ Helper Functions ============
get_optimized_L() {
  local json_file="$1"
  python3 -c "import json; print(json.load(open('${json_file}'))['best_result']['L'])" 2>/dev/null || echo ""
}

get_params_file() {
  local json_file="$1"
  local params_file="${json_file%.json}_params.json"
  if [ -f "$params_file" ]; then
    echo "$params_file"
  else
    echo ""
  fi
}

compute_default_L() {
  local len=$1
  python3 -c "import math; print(int(math.ceil(2 * ${len}**(1/3))))"
}

# ============ Main Script ============
echo "========================================================"
echo "  3D HP Benchmark - Optimized Run"
echo "========================================================"
echo "Phase 1: Optimization (${OPT_TIME_BUDGET}s budget, ${OPT_SEEDS} seeds)"
echo "Phase 2: Final solve (${RUN_TIME_LIMIT}s time limit)"
echo "Workers: ${WORKERS}"
echo "========================================================"
echo

# Summary file
SUMMARY="${RESULTS_DIR}/summary.txt"
echo "3D HP Benchmark Results (Optimized)" > "$SUMMARY"
echo "====================================" >> "$SUMMARY"
echo "Optimization: ${OPT_TIME_BUDGET}s budget, ${OPT_SEEDS} seeds" >> "$SUMMARY"
echo "Final run: ${RUN_TIME_LIMIT}s time limit, ${WORKERS} workers" >> "$SUMMARY"
echo "" >> "$SUMMARY"
printf "%-6s %-4s %-4s %-8s %-10s %-6s %-10s %s\n" \
  "ID" "Len" "L" "Energy" "BestKnown" "Gap" "Status" "Time" >> "$SUMMARY"
echo "---------------------------------------------------------------" >> "$SUMMARY"

total_start=$(date +%s)

for entry in "${cases[@]}"; do
  # Parse entry
  id="${entry%%:*}"
  rest="${entry#*:}"
  len="${rest%%:*}"
  rest2="${rest#*:}"
  seq="${rest2%%:*}"
  best_known="${rest2##*:}"

  opt_file="${OPT_DIR}/${id}.json"
  params_file="${OPT_DIR}/${id}_params.json"
  result_json="${RESULTS_DIR}/${id}.json"
  result_png="${RESULTS_DIR}/${id}.png"

  echo "========================================================"
  echo "  Processing ${id} (len=${len}, best_known=${best_known})"
  echo "========================================================"

  # -------- Phase 1: Optimize (if not already done) --------
  if [ ! -f "$opt_file" ]; then
    echo "[Phase 1] Optimizing settings..."
    uv run python main.py optimize-settings "${seq}" \
      --time-budget "${OPT_TIME_BUDGET}" \
      --workers "${WORKERS}" \
      --seeds "${OPT_SEEDS}" \
      --output "${opt_file}" \
      2>&1
  else
    echo "[Phase 1] Using cached optimization from ${opt_file}"
  fi

  # -------- Extract optimized settings --------
  opt_L=$(get_optimized_L "$opt_file")
  if [ -z "$opt_L" ]; then
    opt_L=$(compute_default_L "$len")
    echo "[Warning] Could not read optimized L, using default: ${opt_L}"
  fi

  params_arg=""
  if [ -f "$params_file" ]; then
    params_arg="--params ${params_file}"
    echo "[Phase 2] Using tuned params from ${params_file}"
  fi

  # -------- Phase 2: Run with optimized settings --------
  echo "[Phase 2] Running solver with L=${opt_L}, time=${RUN_TIME_LIMIT}s..."

  run_start=$(date +%s)

  uv run python main.py run-hp "${seq}" \
    -L "${opt_L}" \
    --dim 3 \
    -t "${RUN_TIME_LIMIT}" \
    -w "${WORKERS}" \
    --viz none \
    --snap "${result_png}" \
    --save-solution "${result_json}" \
    ${params_arg} \
    2>&1

  run_end=$(date +%s)
  run_elapsed=$((run_end - run_start))

  # -------- Extract results --------
  if [ -f "$result_json" ]; then
    energy=$(python3 -c "import json; print(json.load(open('${result_json}'))['energy_epsHH'])")
    contacts=$(python3 -c "import json; print(json.load(open('${result_json}'))['contacts'])")
    gap=$((energy - best_known))

    # Determine status
    if [ "$energy" -eq "$best_known" ]; then
      status="OPTIMAL"
    elif [ "$energy" -lt "$best_known" ]; then
      status="BETTER!"
    else
      status="gap=${gap}"
    fi
  else
    energy="N/A"
    contacts="N/A"
    gap="N/A"
    status="FAILED"
  fi

  printf "%-6s %-4s %-4s %-8s %-10s %-6s %-10s %ss\n" \
    "$id" "$len" "$opt_L" "$energy" "$best_known" "$gap" "$status" "$run_elapsed" >> "$SUMMARY"

  echo ""
  echo "Result: energy=${energy}, best_known=${best_known}, gap=${gap}, status=${status}"
  echo "Saved: ${result_json}, ${result_png}"
  echo ""
done

total_end=$(date +%s)
total_elapsed=$((total_end - total_start))

echo "" >> "$SUMMARY"
echo "Total time: ${total_elapsed}s ($(( total_elapsed / 60 ))m $(( total_elapsed % 60 ))s)" >> "$SUMMARY"

echo "========================================================"
echo "  Complete!"
echo "========================================================"
echo "Total time: ${total_elapsed}s ($(( total_elapsed / 60 ))m $(( total_elapsed % 60 ))s)"
echo ""
echo "Results: ${RESULTS_DIR}/"
echo "Optimized settings: ${OPT_DIR}/"
echo ""
echo "Summary:"
echo "--------"
cat "$SUMMARY"
