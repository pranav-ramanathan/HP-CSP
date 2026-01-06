#!/usr/bin/env bash

set -euo pipefail

# Optimize settings for all 3D benchmark sequences.
# Each sequence gets a configurable time budget for optimization.
# Results are saved to out/optimized/

# Default time budget per sequence (5 minutes)
TIME_BUDGET=${TIME_BUDGET:-300}
WORKERS=${WORKERS:-8}
N_SEEDS=${N_SEEDS:-5}

# Output directory
OUT_DIR="out/optimized"
mkdir -p "$OUT_DIR"

# Format: "id:length:sequence:best_known"
# Sequences from standard 3D HP benchmark set
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

echo "========================================"
echo "3D HP Benchmark Settings Optimization"
echo "========================================"
echo "Time budget per sequence: ${TIME_BUDGET}s"
echo "Workers: ${WORKERS}"
echo "Seeds to try: ${N_SEEDS}"
echo "Output directory: ${OUT_DIR}"
echo "========================================"
echo

# Summary file
SUMMARY="${OUT_DIR}/summary.txt"
echo "3D HP Optimization Summary" > "$SUMMARY"
echo "=========================" >> "$SUMMARY"
echo "Time budget: ${TIME_BUDGET}s, Workers: ${WORKERS}, Seeds: ${N_SEEDS}" >> "$SUMMARY"
echo "" >> "$SUMMARY"
printf "%-6s %-6s %-10s %-10s %-6s %-8s %s\n" "ID" "Len" "Best" "BestKnown" "Gap" "L" "Time" >> "$SUMMARY"
echo "--------------------------------------------------------------" >> "$SUMMARY"

total_start=$(date +%s)

for entry in "${cases[@]}"; do
  # Parse entry
  id="${entry%%:*}"
  rest="${entry#*:}"
  len="${rest%%:*}"
  rest2="${rest#*:}"
  seq="${rest2%%:*}"
  best_known="${rest2##*:}"

  output_file="${OUT_DIR}/${id}.json"

  echo "=== Optimizing ${id} (len=${len}, best_known=${best_known}) ==="

  start=$(date +%s)

  uv run python main.py optimize-settings "${seq}" \
    --time-budget "${TIME_BUDGET}" \
    --workers "${WORKERS}" \
    --seeds "${N_SEEDS}" \
    --output "${output_file}" \
    2>&1

  end=$(date +%s)
  elapsed=$((end - start))

  # Extract best energy from JSON
  if [ -f "$output_file" ]; then
    best_energy=$(python3 -c "import json; print(json.load(open('${output_file}'))['best_result']['energy'])")
    best_L=$(python3 -c "import json; print(json.load(open('${output_file}'))['best_result']['L'])")
    gap=$((best_energy - best_known))
  else
    best_energy="N/A"
    best_L="N/A"
    gap="N/A"
  fi

  printf "%-6s %-6s %-10s %-10s %-6s %-8s %ss\n" \
    "$id" "$len" "$best_energy" "$best_known" "$gap" "$best_L" "$elapsed" >> "$SUMMARY"

  echo "=== ${id} done: best=${best_energy}, best_known=${best_known}, gap=${gap} ==="
  echo
done

total_end=$(date +%s)
total_elapsed=$((total_end - total_start))

echo "" >> "$SUMMARY"
echo "Total time: ${total_elapsed}s" >> "$SUMMARY"

echo "========================================"
echo "Optimization Complete!"
echo "========================================"
echo "Total time: ${total_elapsed}s"
echo "Results saved to: ${OUT_DIR}/"
echo "Summary: ${SUMMARY}"
echo
cat "$SUMMARY"
