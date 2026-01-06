#!/usr/bin/env bash

set -euo pipefail

# 3D HP Benchmark - Full Pipeline with Optimization
# Total budget: 7 hours (25200 seconds) for 9 sequences
# Per sequence: ~46 minutes (2800 seconds)
#   - Optimization: ~10 minutes (600 seconds)
#   - Final solve: ~36 minutes (2200 seconds)

# ============ Configuration ============
TOTAL_BUDGET_S=25200  # 7 hours in seconds
NUM_SEQUENCES=9

# Calculate per-sequence budgets
PER_SEQ_BUDGET=$((TOTAL_BUDGET_S / NUM_SEQUENCES))  # ~2800s per sequence
OPT_BUDGET=$((PER_SEQ_BUDGET * 20 / 100))           # 20% for optimization (~560s)
SOLVE_BUDGET=$((PER_SEQ_BUDGET * 80 / 100))         # 80% for final solve (~2240s)

# Override with environment variables if set
OPT_TIME_BUDGET=${OPT_TIME_BUDGET:-$OPT_BUDGET}
RUN_TIME_LIMIT=${RUN_TIME_LIMIT:-$SOLVE_BUDGET}
WORKERS=${WORKERS:-8}
OPT_SEEDS=${OPT_SEEDS:-7}

# Directories
OUT_DIR="out"
OPT_DIR="${OUT_DIR}/optimized"
RESULTS_DIR="${OUT_DIR}/results"
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

compute_default_L() {
  local len=$1
  python3 -c "import math; print(int(math.ceil(2 * ${len}**(1/3))))"
}

format_time() {
  local seconds=$1
  printf "%dh %dm %ds" $((seconds/3600)) $((seconds%3600/60)) $((seconds%60))
}

# ============ Main Script ============
echo "========================================================"
echo "  3D HP Benchmark - Full Pipeline"
echo "========================================================"
echo "Total budget: $(format_time $TOTAL_BUDGET_S)"
echo "Per sequence: $(format_time $PER_SEQ_BUDGET)"
echo "  - Optimization: $(format_time $OPT_TIME_BUDGET) (seeds: $OPT_SEEDS)"
echo "  - Final solve:  $(format_time $RUN_TIME_LIMIT)"
echo "Workers: ${WORKERS}"
echo "========================================================"
echo

# Summary file
SUMMARY="${RESULTS_DIR}/summary.txt"
cat > "$SUMMARY" << EOF
3D HP Benchmark Results
=======================
Total budget: $(format_time $TOTAL_BUDGET_S)
Per sequence: $(format_time $PER_SEQ_BUDGET)
  - Optimization: $(format_time $OPT_TIME_BUDGET)
  - Final solve:  $(format_time $RUN_TIME_LIMIT)
Workers: ${WORKERS}

EOF
printf "%-6s %-4s %-4s %-8s %-10s %-6s %-10s %-12s %s\n" \
  "ID" "Len" "L" "Energy" "BestKnown" "Gap" "Status" "SolveTime" "TotalTime" >> "$SUMMARY"
echo "-------------------------------------------------------------------------" >> "$SUMMARY"

total_start=$(date +%s)
sequences_done=0
total_contacts=0
optimal_count=0

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

  seq_start=$(date +%s)

  echo ""
  echo "========================================================"
  echo "  [$((sequences_done + 1))/${NUM_SEQUENCES}] ${id} (len=${len}, best_known=${best_known})"
  echo "========================================================"

  # -------- Phase 1: Optimize --------
  if [ ! -f "$opt_file" ] || [ "${FORCE_REOPTIMIZE:-}" = "1" ]; then
    echo "[Phase 1] Optimizing settings (budget: $(format_time $OPT_TIME_BUDGET))..."

    uv run python main.py optimize-settings "${seq}" \
      --time-budget "${OPT_TIME_BUDGET}" \
      --workers "${WORKERS}" \
      --seeds "${OPT_SEEDS}" \
      --output "${opt_file}" \
      2>&1 | grep -v "^Building\|^Declaring\|^Adding\|^Fixing\|^Edge-based\|^Setting objective\|^Objective set"
  else
    echo "[Phase 1] Using cached optimization: ${opt_file}"
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
    echo "[Phase 2] Using tuned params: ${params_file}"
  fi

  # -------- Phase 2: Final solve --------
  echo "[Phase 2] Running solver (L=${opt_L}, time: $(format_time $RUN_TIME_LIMIT))..."

  solve_start=$(date +%s)

  uv run python main.py run-hp "${seq}" \
    -L "${opt_L}" \
    --dim 3 \
    -t "${RUN_TIME_LIMIT}" \
    -w "${WORKERS}" \
    --viz none \
    --snap "${result_png}" \
    --save-solution "${result_json}" \
    ${params_arg} \
    2>&1 | grep -v "^Building\|^Declaring\|^Adding\|^Fixing\|^Edge-based\|^Setting objective\|^Objective set"

  solve_end=$(date +%s)
  solve_elapsed=$((solve_end - solve_start))
  seq_end=$(date +%s)
  seq_elapsed=$((seq_end - seq_start))

  # -------- Extract and report results --------
  if [ -f "$result_json" ]; then
    energy=$(python3 -c "import json; print(json.load(open('${result_json}'))['energy_epsHH'])")
    contacts=$(python3 -c "import json; print(json.load(open('${result_json}'))['contacts'])")
    gap=$((energy - best_known))
    total_contacts=$((total_contacts + contacts))

    if [ "$energy" -eq "$best_known" ]; then
      status="OPTIMAL"
      optimal_count=$((optimal_count + 1))
    elif [ "$energy" -lt "$best_known" ]; then
      status="BETTER!"
      optimal_count=$((optimal_count + 1))
    else
      status="gap=${gap}"
    fi
  else
    energy="N/A"
    contacts="0"
    gap="N/A"
    status="FAILED"
  fi

  printf "%-6s %-4s %-4s %-8s %-10s %-6s %-10s %-12s %s\n" \
    "$id" "$len" "$opt_L" "$energy" "$best_known" "$gap" "$status" \
    "$(format_time $solve_elapsed)" "$(format_time $seq_elapsed)" >> "$SUMMARY"

  sequences_done=$((sequences_done + 1))
  elapsed_total=$((seq_end - total_start))
  remaining_est=$(( (TOTAL_BUDGET_S - elapsed_total) > 0 ? (TOTAL_BUDGET_S - elapsed_total) : 0 ))

  echo ""
  echo "Result: energy=${energy}, best_known=${best_known}, status=${status}"
  echo "Time: solve=$(format_time $solve_elapsed), total=$(format_time $seq_elapsed)"
  echo "Progress: ${sequences_done}/${NUM_SEQUENCES}, elapsed=$(format_time $elapsed_total)"
  echo "Saved: ${result_json}"
done

total_end=$(date +%s)
total_elapsed=$((total_end - total_start))

# -------- Final Summary --------
cat >> "$SUMMARY" << EOF
-------------------------------------------------------------------------

Total time: $(format_time $total_elapsed)
Sequences: ${sequences_done}/${NUM_SEQUENCES}
Optimal/Better: ${optimal_count}/${NUM_SEQUENCES}
Total contacts: ${total_contacts}
EOF

echo ""
echo "========================================================"
echo "  COMPLETE"
echo "========================================================"
echo "Total time: $(format_time $total_elapsed)"
echo "Optimal/Better: ${optimal_count}/${NUM_SEQUENCES}"
echo "Total contacts: ${total_contacts}"
echo ""
echo "Results: ${RESULTS_DIR}/"
echo "Optimized settings: ${OPT_DIR}/"
echo ""
echo "Summary:"
echo "--------"
cat "$SUMMARY"
