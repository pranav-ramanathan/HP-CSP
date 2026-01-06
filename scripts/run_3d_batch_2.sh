#!/usr/bin/env bash

set -euo pipefail

# Template: 3D HP Benchmark - Full Pipeline with Optimization
# Fill the `cases` array with your sequences and adjust budgets/paths below.
# Each case entry uses the format: "id:length:sequence:best_known"

# ============ User Configuration ============

# Example cases (replace with your own)
cases=(
  "case1:60:P(PH3)2H5P3H10PHP3H12P4H6PH2PHP:-55"
)

# Total time budget (seconds) for all sequences
TOTAL_BUDGET_S=${TOTAL_BUDGET_S:-3600}
# Percent of per-sequence budget to spend on optimization (rest for final solve)
OPT_RATIO=${OPT_RATIO:-50}
# Solver workers and seeds for optimization
WORKERS=${WORKERS:-8}
OPT_SEEDS=${OPT_SEEDS:-7}

# Commands (switch PY_CMD to "python" if you are not using uv)
UV_BIN=${UV_BIN:-uv}
PY_CMD=${PY_CMD:-"$UV_BIN run python"}

# Directories (override via environment variables if desired)
OUT_DIR=${OUT_DIR:-"out/correction"}
OPT_DIR=${OPT_DIR:-"${OUT_DIR}/optimized"}
RESULTS_DIR=${RESULTS_DIR:-"${OUT_DIR}/results"}
mkdir -p "$OPT_DIR" "$RESULTS_DIR"

# Allow overrides for phase budgets
OPT_TIME_BUDGET=${OPT_TIME_BUDGET:-0}  # 0 means compute from TOTAL_BUDGET_S
RUN_TIME_LIMIT=${RUN_TIME_LIMIT:-0}    # 0 means compute from TOTAL_BUDGET_S

# ============ Helper Functions ============
get_optimized_L() {
  local json_file="$1"
  python3 - <<'PY' "$json_file" 2>/dev/null || true
import json, sys
path = sys.argv[1]
print(json.load(open(path))['best_result']['L'])
PY
}

compute_default_L() {
  local len=$1
  python3 - <<'PY' "$len"
import math, sys
length = int(sys.argv[1])
print(int(math.ceil(2 * length ** (1/3))))
PY
}

format_time() {
  local seconds=$1
  printf "%dh %dm %ds" $((seconds/3600)) $((seconds%3600/60)) $((seconds%60))
}

# ============ Budget Calculation ============
NUM_SEQUENCES=${#cases[@]}
if [ "$NUM_SEQUENCES" -eq 0 ]; then
  echo "Error: No cases defined. Populate the 'cases' array." >&2
  exit 1
fi

PER_SEQ_BUDGET=$((TOTAL_BUDGET_S / NUM_SEQUENCES))
OPT_BUDGET=$((PER_SEQ_BUDGET * OPT_RATIO / 100))
SOLVE_BUDGET=$((PER_SEQ_BUDGET * (100 - OPT_RATIO) / 100))

if [ "$OPT_TIME_BUDGET" -eq 0 ]; then OPT_TIME_BUDGET=$OPT_BUDGET; fi
if [ "$RUN_TIME_LIMIT" -eq 0 ]; then RUN_TIME_LIMIT=$SOLVE_BUDGET; fi

# ============ Intro ============
echo "========================================================"
echo "  3D HP Benchmark (Template)"
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

# ============ Main Loop ============
for entry in "${cases[@]}"; do
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
    ${PY_CMD} main.py optimize-settings "${seq}" \
      --time-budget "${OPT_TIME_BUDGET}" \
      --workers "${WORKERS}" \
      --seeds "${OPT_SEEDS}" \
      --output "${opt_file}" \
      2>&1 | grep -v "^Building\\|^Declaring\\|^Adding\\|^Fixing\\|^Edge-based\\|^Setting objective\\|^Objective set"
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

  ${PY_CMD} main.py run-hp "${seq}" \
    -L "${opt_L}" \
    --dim 3 \
    -t "${RUN_TIME_LIMIT}" \
    -w "${WORKERS}" \
    --viz none \
    --snap "${result_png}" \
    --save-solution "${result_json}" \
    ${params_arg} \
    2>&1 | grep -v "^Building\\|^Declaring\\|^Adding\\|^Fixing\\|^Edge-based\\|^Setting objective\\|^Objective set"

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


