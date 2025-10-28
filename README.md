## HP Protein Folding with OR-Tools CP-SAT (2D/3D)

Constraint Programming solver for the HP (Hydrophobic–Polar) model on 2D and 3D lattices, with a modern CLI, compressed-sequence input, edge-based contacts, native 3D visualization, batch evaluation, and solve-time reporting.

### Highlights
- 2D (4-neighborhood) and 3D (6-neighborhood) lattices
- One-hot placement per residue; exact-one and at-most-one constraints
- Chain adjacency via neighbor implications; symmetry fixing (start at center, next “up”)
- Edge-based H–H contacts (memory-lean): one contact var per lattice edge (not per (i,j,u,v))
- Compressed HP notation (ASCII and Unicode subscripts), with strict mode and length assertion
- Native 3D viewer (PyVista preferred; vedo fallback) and Plotly HTML export
- Batch evaluator with CSV/JSON/LaTeX outputs; energy reporting E/ε_HH = −(contacts)
- Typer + Rich CLI; solve wall time printed after each solve

## Install

```bash
uv sync
```

## Quickstart

### 2D (default)
```bash
uv run main.py run-hp "HPPH" -L 4
```

### 3D
```bash
uv run main.py run-hp "HPPH" -L 4 --dim 3
```

### With visualization
```bash
uv run main.py run-hp "HPPH" -L 6 --dim 3 --viz 3d
uv run main.py run-hp "HPPH" -L 6 --dim 3 --viz native --snap out/fold.png
```

### Use tuned CP-SAT parameters

```bash
# Tune on a sequence set (20 min budget; 60s per solve)
uv run main.py tune-cpsat -s HPHPPH -L 6 --dim 3 --time 60 --workers 8 --budget 1200

# Run with tuned params (auto-resolve by dim/L/n)
uv run main.py run-hp "HPHPPH" -L 6 --dim 3 --time 60 --use-tuned

# Or provide an explicit params file
uv run main.py run-hp "HPHPPH" -L 6 --dim 3 --time 60 --params out/tuned/hp_3d_L6_n6.json

# Batch with tuned params
uv run main.py batch --file sequences.txt -L 6 --dim 3 --time 60 --use-tuned
```

## Compressed Input (ASCII and Unicode subscripts)

Items: `H` or `P` optionally followed by a count (ASCII `0-9` or Unicode subscripts). Whitespace/newlines/commas are ignored by default.

Examples:
- `HPH2PH9P3H10` → expands to `H P H H P H…`
- `HPH₂PH₉P₃H₁₀` → identical expansion
- Spaced: `HPH2 PH9, P3\nH10` → identical expansion

Flags:
- `--print-expanded` prints the expanded string (clip to 200 chars)
- `--strict` errors on non-HP/non-digit/subscript chars
- `--expect-length N` asserts expanded length is `N`

## Energy & Objective

- Contacts are per undirected lattice edge `{u,v}` with both endpoints H and occupied, and not used by the chain.
- Objective: maximize `Σ C_edge`; energy in ε_HH units is `E = − Σ C_edge`.
- Console prints contacts and energy; also prints `Solve time: …s` after each solve.

## Batch Evaluation

```bash
uv run main.py batch --file sequences.txt -L 6 --dim 3 --time 120 --workers 8 --runs 3 \
  --csv out/hp_results.csv --json out/hp_results.json --latex out/hp_results.tex
```

Per sequence summary: best energy, contacts, optimal?, median contacts, best-run time; CSV/JSON/LaTeX exports supported.

## Choosing L, Workers, and Time

- 2D: `L ≈ ceil(sqrt(n)) + 1`
- 3D: `L ≈ ceil(n^(1/3)) + 2`
- Workers: typical laptop 6–8; reduce when multitasking.
- Time limit is a cap; solver can stop earlier if proven optimal/infeasible.

## Files

- `main.py` — CLI (2D/3D, viz, batch); prints wall time
- `hp2d.py` — 2D grid, neighbors, chain, edge-based contacts
- `hp3d.py` — 3D grid, neighbors, chain, edge-based contacts
- `hp_parse.py` — Compressed input parser and summary
- `hp_eval.py` — Batch evaluator (energy, CSV/JSON, LaTeX)
- `viz_native.py` — Native 3D viz (PyVista/vedo), energy overlays; PNG snaps

## Troubleshooting

- OOM / `basic_string`: grid too large or sequence very long. Reduce `-L`, lower `--workers`, or try 2D.
- Native viewer requires `pyvista`/`vtk` or uses `vedo` fallback.
- Plotly HTML uses CDN by default; ensure it’s allowed.

## License

See LICENSE.

