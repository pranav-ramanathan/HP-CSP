# CSP Codebase Explained

This guide explains every part of the CSP (Constraint Programming) HP protein folding repository. It covers what each file does, how the system works end-to-end, the constraint model, key data structures, algorithms, CLI usage, and extensibility points.

The project implements the classic HP (Hydrophobic–Polar) lattice protein folding model via OR-Tools CP-SAT in both 2D and 3D.


## High-level Overview

- Purpose: Given an HP sequence, place residues on a lattice to maximize the number of non-consecutive H–H contacts (equivalently, minimize energy E/ε_HH = −contacts).
- Approach: Constraint Programming with Boolean one-hot placement, adjacency constraints for the backbone, at-most-one occupancy, symmetry breaking, and an edge-based contact formulation for memory efficiency. The objective is to maximize the sum of H–H contact variables.
- Dimensions: 2D (4-neighborhood) and 3D (6-neighborhood) are supported with parallel APIs.
- CLI: A modern Typer-based interface in `main.py` provides single-run solving, batch evaluation, visualization, and saved solution loading.
- Visualization: Plotly HTML for web-like 3D viewing and a native desktop 3D viewer using PyVista (preferred) or vedo.


## Repository Map

- `main.py`: CLI entrypoint. Orchestrates parsing, building and solving models, visualization, batch mode, and viewing saved solutions.
- `hp_parse.py`: Compressed HP input parser (supports ASCII digits and Unicode subscripts). Also provides summary metrics.
- `hp2d.py`: 2D lattice helpers, placement/chain/contact constraints, solver configuration, and extraction utilities.
- `hp3d.py`: 3D lattice helpers, 3D analogs of constraints, Plotly HTML 3D renderer, and debug slice printer.
- `hp_eval.py`: Batch evaluation, energy computation, timing, results aggregation, and CSV/JSON/LaTeX export.
- `viz_native.py`: Native 3D visualization (PyVista or vedo) and offscreen PNG snapshots.
- `tune_params.py`: Utilities to apply/load/save tuned CP-SAT parameter dictionaries and key resolution.
- `tools/tune_cpsat.py`: Tuning harness integrating `cpsat-autotune` with an evaluation loop.
- `pyproject.toml`: Project metadata and dependencies.
- `README.md`: Quickstart, feature highlights, flags, and troubleshooting.


## End-to-End Execution Flow

### Single run (2D/3D)
1. CLI command:
   - `uv run main.py run-hp "HPPH" -L 6 --dim 3 --viz native --snap out/fold.png`
2. `main.run_hp(...)`:
   - Parses compressed/expanded input via `hp_parse.expand_hp_sequence` (prints hydrophobic fraction; optional strict mode and length assertion).
   - Branches on `dim` to 2D or 3D.
3. Model construction (2D via `hp2d`, 3D via `hp3d`):
   - Build lattice cells and neighbor map: `build_grid`/`build_grid_3d`, `neighbors_map`/`neighbors_map_3d`.
   - Symmetry anchors: `center_and_up`/`center_and_up_3d`.
   - Declare one-hot placement variables `X[i][c]`.
   - Add constraints: exactly-one per residue, at-most-one per cell, chain adjacency implications, symmetry fixing of first two residues.
   - Declare edge-based contact variables C_edge subject to occupancy, H-only, and excluding chain edges.
   - Objective: maximize `Σ C_edge` (equivalent to minimizing energy E/ε_HH = −Σ C_edge).
4. Solve with CP-SAT:
   - Configure solver parameters (time limit, workers).
   - Solve and report wall time and status.
5. Extract solution:
   - Positions `pos[i]` from one-hot variables.
   - Realized contacts from contact edge variables.
   - Energy computed as `-len(contacts)`.
6. Output/visualize:
   - Print summary to console; for 2D, also a simple grid.
   - Optional snapshot PNG via `viz_native.save_native_3d`.
   - Optional interactive native viewer (`viz_native.render_native_3d`) or Plotly HTML (`hp3d.render_3d_html`).
   - Optional solution JSON with positions and contacts.

### Batch evaluation
1. `main.batch(...)` loads sequences (inline or from file), validates them, and calls `hp_eval.evaluate_batch(...)`.
2. `hp_eval.evaluate_batch(...)` runs `evaluate_sequence(...)` multiple times per sequence (supports seeds), aggregates best/median results. Optional `params`, `use_tuned`, and `tuned_dir` allow applying tuned CP-SAT parameter sets.
3. Prints a Rich table summary and optionally exports CSV, JSON, and LaTeX.

### View saved solution
- `main.view_solution(path, viz=...)` loads the JSON produced by `--save-solution` and renders via Plotly HTML or the native viewer. 2D solutions are lifted to z=0 for 3D rendering.


### Control/Data Flow Diagram (Run)

```
User CLI → main.run_hp
  ├─ hp_parse.expand_hp_sequence → expanded seq, warnings, hydrophobic fraction
  └─ (dim==2? hp2d : hp3d)
      ├─ build_grid[_3d], neighbors_map[_3d], center_and_up[_3d]
      ├─ declare_vars[_3d] → X[i][c] one-hot
      ├─ add_exactly_one_cell_per_residue
      ├─ add_at_most_one_residue_per_cell
      ├─ add_chain_adjacency_allowed_pairs[_3d]
      ├─ fix_start_and_orientation[_3d]
      ├─ declare_hh_contact_vars[_3d] → C_edge
      ├─ set_objective_max_contacts(Σ C_edge)
      ├─ configure_solver → time_limit, workers
      ├─ CP-SAT solve → status, wall time
      ├─ extract positions, realized contacts, energy = -|contacts|
      └─ output: console, (optional) JSON, (optional) Plotly/native viz, (optional) PNG
```


## HP Model and Constraint Formulation

### Variables
- Placement: `X[i][c] ∈ {0,1}` — residue `i` is placed at cell `c`.
- Chain edge usage indicators: `B[(i,u,v)] ∈ {0,1}` — backbone uses directed neighbor edge `(u→v)` between residues `i` and `i+1`.
- Occupancy: `Occ[u] ∈ {0,1}` — any residue occupies cell `u`.
- Hydrophobic occupancy: `H_at[u] ∈ {0,1}` — an `H` residue occupies `u`.
- Edge-based contact variables: `C_edge[{u,v}] ∈ {0,1}` — `u` and `v` are adjacent lattice cells, both occupied by `H`, and not used by the chain.

### Constraints
- Exactly one cell per residue: `∀i, Σ_c X[i][c] == 1`.
- At most one residue per cell: `∀c, Σ_i X[i][c] ≤ 1`.
- Chain adjacency (one-hot implication form): `X[i][u] == 1 ⇒ Σ_{v∈N(u)} X[i+1][v] == 1`.
- Symmetry breaking: Fix residue 0 at center; residue 1 at an “up” neighbor — reduces rotational/reflective degeneracy.
- Edge-based contacts: `C_edge[{u,v}] ≤ Occ[u], Occ[v], H_at[u], H_at[v]` and `C_edge[{u,v}] ≤ 1 - B[(i,u,v)]` and `≤ 1 - B[(i,v,u)]` for all chain edges — excludes counting backbone edges as contacts.

### Objective and Energy
- Objective: Maximize `Σ C_edge`.
- Energy: `E/ε_HH = - Σ C_edge` (most negative is best).

### 2D vs 3D
- 2D lattice uses 4-neighborhood; 3D uses 6-neighborhood with identical patterns for variables and constraints.
- Symmetry anchors differ in axis: 2D fixes “up” along Y; 3D fixes “up” along Z.


## Algorithms & Heuristics

- Model class: Constraint Programming with CP-SAT; one-hot encoding of residue placements.
- Chain enforcement: Implication constraints ensure each consecutive pair occupies adjacent cells in the lattice.
- Contacts (edge-based): One Boolean per undirected lattice edge `{u,v}`; activated only if both endpoints are occupied by `H` and the edge is not used by the chain. This avoids O(n^2|cells|^2) pair-location variables and scales with the number of lattice edges.
- Objective: Maximize sum of contact edge variables; energy is defined as `E/ε_HH = −Σ C_edge`.
- Symmetry breaking: Fixes residue 0 at lattice center and residue 1 “up” to reduce symmetric duplicates and search space.
- Complexity: Still exponential in the worst case; performance depends on `n`, `L`, and density of H’s. Edge-based contacts are memory-lean relative to per-pair formulations. Parallel workers can accelerate search.
- Practical tips: Use modest `L` just above the minimal embedding requirement, constrain time, and leverage batch runs to sample seeds when looking for good solutions quickly.


## Module Deep Dives

### `main.py` — CLI Orchestrator
- Typer app exposing commands:
  - `run-hp`: Solve a single sequence with options for grid size `-L`, dimension `-d {2,3}`, time limit, workers, visualization (`none|3d|native`), snapshot path, strict parsing, expected length, and saving a solution JSON.
  - `batch`: Evaluate multiple sequences (inline or from file) with runs, seeds, and export options (CSV/JSON/LaTeX). Aggregates best energy and median contacts.
  - `view_solution`: Load a solution JSON and visualize via Plotly HTML or native viewer; lifts 2D to 3D by adding `z=0`.
- 2D path (`run_hp_2d`): Uses `hp2d` helpers and constraints; extracts positions and realized contact edges; computes energy via `hp_eval.compute_energy_contacts` and prints a 2D grid.
- 3D path (`run_hp_3d`): Uses `hp3d` helpers and constraints; extracts positions and realized contact edges; computes energy as `-contacts_count`. Also prints 3D Z-slices for debugging. Visualization via Plotly or native.
- Solution saving: JSON contains `sequence`, `L`, `dim`, `positions`, `contacts_edges`, `contacts`, `energy_epsHH`, and `status`.

### `hp_parse.py` — Compressed Sequence Parser
- Supports ASCII digits and Unicode subscripts (e.g., `H₂P₁₀`), allows separators/whitespace by default, and has a strict mode to error on unexpected characters.
- `expand_hp_sequence(raw, strict=False) -> (expanded, warnings)`: Expands to a pure `H`/`P` string; rejects mixed digit styles, zero repeats, and invalid letters.
- `hydrophobic_fraction(expanded)`: `#H / len(expanded)`.
- Raises `ParseError` on invalid inputs.

### `hp2d.py` — 2D Model Components
- `build_grid(L)`: All `(x,y)` cells in `[0..L-1]^2`.
- `neighbors_map(cells)`: 4-neighborhood map.
- `center_and_up(L)`: Returns the grid center and an “up” neighbor cell.
- `declare_vars(model, n, cells)`: Builds one-hot terms `X[i][c]` for all residues and cells.
- Constraints: `add_exactly_one_cell_per_residue`, `add_at_most_one_residue_per_cell`, `add_chain_adjacency_allowed_pairs` (implication pattern), `fix_start_and_orientation` (symmetry breaking).
- Contacts: `declare_hh_contact_vars` builds `B[(i,u,v)]`, `Occ[u]`, `H_at[u]`, and `C_edge[{u,v}]` with chain-edge exclusion.
- Objective: `set_objective_max_contacts(model, C_edge.values())`.
- Solver and extraction: `configure_solver(solver, time_limit_s, workers)`, `extract_positions(solver, X)`.

### `hp3d.py` — 3D Model Components
- `build_grid_3d(L)`: All `(x,y,z)` cells in `[0..L-1]^3`.
- `neighbors_map_3d(cells)`: 6-neighborhood map.
- `center_and_up_3d(L)`: Center `(cx,cy,cz)` and “up” `(cx,cy,cz-1)` (or +1).
- `declare_vars_3d(model, n, cells)`: 3D one-hot terms.
- 3D constraints mirror 2D: `add_exactly_one_cell_per_residue_3d`, `add_at_most_one_residue_per_cell_3d`, `add_chain_adjacency_allowed_pairs_3d`, `fix_start_and_orientation_3d`.
- 3D contacts: `declare_hh_contact_vars_3d` with the same edge-based scheme.
- Visualization: `render_3d_html` generates an interactive Plotly HTML; `print_3d_slices` prints Z-slices to the console.

### `hp_eval.py` — Energy, Timing, Batch, and Exports
- Energy: `compute_energy_contacts(contacts) = -len(contacts)`.
- Timing: `timeit(label)` decorator and `solve_with_timer` measure wall time for solving.
- Single evaluation: `evaluate_sequence(seq, L, dim, time_limit_s, workers, seed)` builds a model (2D or 3D), solves it, extracts positions and contacts, computes energy, and returns a detailed dictionary.
- Batch mode: `evaluate_batch(sequences, L, dim, time_limit_s, workers, seeds, runs_per_seq)` — multiple runs per sequence; selects best energy; computes medians and timings.
- Exports: `save_results_csv`, `save_results_json`, and `latex_table` (tabular for LaTeX). `print_summary_table` renders a Rich console table.

### `viz_native.py` — Native 3D Viewer and Snapshots
- Backend selection: PyVista preferred; vedo fallback; otherwise the functions become no-ops with user guidance.
- `render_native_3d(seq, positions, contacts_hh, ...)`: Interactive window with spheres (residues), tubes (backbone), optional contact tubes, overlays, and camera controls.
- `save_native_3d(seq, positions, path_png, ...)`: Offscreen snapshot to PNG with overlays (supports both backends).
- Internal implementations: `_render_pyvista`, `_save_pyvista`, `_render_vedo`, `_save_vedo`.


## Key Data Structures and Conventions

- Sequence string: `'H'`/`'P'` only after expansion. Use `hp_parse.expand_hp_sequence` to normalize inputs consistently.
- Indices: Residues are zero-indexed consistently across the codebase.
- Positions: `Dict[int, Tuple[int,int]]` in 2D or `Dict[int, Tuple[int,int,int]]` in 3D. 2D is lifted to 3D as `(x,y,0)` for visualization.
- Contacts:
  - Edge-based representation keyed by undirected lattice edges `{u,v}` after normalization `(min(u,v), max(u,v))`.
  - Realized contacts are extracted by checking the value of `C_edge` after solving.
- Energy: Reported as `E/ε_HH` (integer), where `ε_HH` is the unit for each H–H contact.
- Logs: Rich console messages describe each model build step for transparency.


## CLI Reference and Usage

### Installation
```
uv sync
```

### Run a single sequence
- 2D (default):
```
uv run main.py run-hp "HPPH" -L 4
```
- 3D:
```
uv run main.py run-hp "HPPH" -L 4 --dim 3
```
- With visualization:
```
uv run main.py run-hp "HPPH" -L 6 --dim 3 --viz 3d
uv run main.py run-hp "HPPH" -L 6 --dim 3 --viz native --snap out/fold.png
```

Key flags:
- `--grid-size, -L`: Lattice dimension (2D: `L×L`, 3D: `L×L×L`).
- `--time, -t`: CP-SAT time limit in seconds.
- `--workers, -w`: Number of parallel search workers.
- `--dim, -d {2,3}`: 2D or 3D lattice.
- `--viz {none,3d,native}`: No viz, Plotly HTML, or native viewer.
- `--snap PATH`: Save PNG snapshot (native viewer offscreen render).
- `--strict`: Strict parsing of compressed input.
- `--print-expanded`: Print expanded sequence (clipped to 200 chars).
- `--expect-length N`: Assert expansion length equals `N`.
- `--save-solution PATH`: Save solution JSON (positions, contacts, energy).

### Batch evaluation
```
uv run main.py batch --file sequences.txt -L 6 --dim 3 --time 120 --workers 8 --runs 3 \
  --csv out/hp_results.csv --json out/hp_results.json --latex out/hp_results.tex
```

Outputs per sequence: best energy/contacts/status/optimal?, median contacts, best-run time. Optional CSV/JSON/LaTeX exports.

### Load and view a saved solution
```
uv run main.py view-solution ./out/solution.json --viz native --snap out/view.png
```


## Choosing Parameters

- Suggested lattice sizes:
  - 2D: `L ≈ ceil(sqrt(n)) + 1`
  - 3D: `L ≈ ceil(n^(1/3)) + 2`
- `workers`: Use 6–8 on a typical laptop; reduce when multitasking.
- `time`: A cap; solver stops early if proven optimal/infeasible.


## Extensibility Guide

- Alternative objectives:
  - Replace `set_objective_max_contacts` with weighted contacts or additional structural terms.
- Additional constraints:
  - Add angle/bend penalties, compactness constraints, or exclude specific motifs.
- Heuristics/assumptions:
  - Tighter symmetry breaking; preclude inverted duplicates; pin additional residues.
- New visualizations:
  - Extend `viz_native.py` for new overlays; add more export formats (SVG, GLB) in `hp3d.py`.
- New metrics/exports:
  - Add CSV columns or additional aggregations in `hp_eval.py`.


## Troubleshooting and Limitations

- Memory blowups / very long sequences:
  - The one-hot and edge-based contact model scales with `n * |cells|` and number of lattice edges. Use conservative `L`, fewer `workers`, or switch to 2D for faster runs.
- Visualization dependencies:
  - Native viewer requires `pyvista`/`vtk` or falls back to `vedo`. Plotly HTML uses CDN by default.
- Parser pitfalls:
  - Mixed digit styles (`H₂3`) are invalid; zero repeats are not allowed; strict mode errors on unknown non-whitespace characters.
- 2D vs 3D energy reporting:
  - Both use `E = -contacts` (in 2D via `hp_eval.compute_energy_contacts`, in 3D energy is computed inline as `-contacts_count`). They are consistent in meaning.


## Appendix: Solution JSON Schema

A typical saved solution (2D example) contains:

```json
{
  "sequence": "HPPH...",
  "L": 6,
  "dim": 2,
  "positions": {"0": [3,3], "1": [3,2], ...},
  "contacts_edges": [ [[x1,y1],[x2,y2]], ... ],
  "contacts": 5,
  "energy_epsHH": -5,
  "status": "saved"
}
```

3D `positions` and `contacts_edges` include `z` coordinates.


---

If you are new to the codebase, start with `main.py` to see how a single run is orchestrated, then read `hp2d.py` or `hp3d.py` to understand the constraint model, followed by `hp_eval.py` for batch workflows, `hp_parse.py` for input handling, and `viz_native.py`/`hp3d.py` for rendering.

## Per-Module API Summary

### `main.py`
- `run_hp(...)`: Unified CLI for 2D/3D solving; handles parsing, solving, visualization, saving.
- `run_hp_2d(...)`: Builds 2D model with `hp2d`, solves, computes energy, prints grid, handles native snapshot/viewer.
- `run_hp_3d(...)`: Builds 3D model with `hp3d`, solves, prints Z-slices, renders Plotly/native, handles snapshots.
- `batch(...)`: Batch evaluation and optional CSV/JSON/LaTeX exports via `hp_eval`.
- `view_solution(...)`: Load solution JSON and display via Plotly HTML or native viewer.

### `hp2d.py`
- Grid & neighbors: `build_grid`, `neighbors_map`, `center_and_up`.
- Variables: `declare_vars` → `Dict[int, Dict[(x,y), BoolVar]]`.
- Constraints: `add_exactly_one_cell_per_residue`, `add_at_most_one_residue_per_cell`, `add_chain_adjacency_allowed_pairs`, `fix_start_and_orientation`.
- Contacts: `declare_hh_contact_vars` → contact edge vars with chain exclusion.
- Objective: `set_objective_max_contacts`.
- Solver utils: `configure_solver`, `extract_positions`.

### `hp3d.py`
- Grid & neighbors: `build_grid_3d`, `neighbors_map_3d`, `center_and_up_3d`.
- Variables: `declare_vars_3d`.
- Constraints: `add_exactly_one_cell_per_residue_3d`, `add_at_most_one_residue_per_cell_3d`, `add_chain_adjacency_allowed_pairs_3d`, `fix_start_and_orientation_3d`.
- Contacts: `declare_hh_contact_vars_3d`.
- Visualization: `render_3d_html`, `print_3d_slices`.

### `hp_eval.py`
- Energy: `compute_energy_contacts` (= -len(contacts)).
- Timing: `timeit`, `solve_with_timer`.
- Single run: `evaluate_sequence`.
- Batch: `evaluate_batch`.
- Export: `save_results_csv`, `save_results_json`, `latex_table`, `print_summary_table`.

### `hp_parse.py`
- Errors: `ParseError`.
- Parsing: `expand_hp_sequence(raw, strict=False)`.
- Metrics: `hydrophobic_fraction`.

### `viz_native.py`
- Rendering: `render_native_3d` (PyVista or vedo).
- Snapshot: `save_native_3d` (offscreen).
- Internals: `_render_pyvista/_save_pyvista`, `_render_vedo/_save_vedo`.

## Usage & Outputs

- Install deps: `uv sync`
- Single run (2D): `uv run main.py run-hp "HPPH" -L 4`
- Single run (3D): `uv run main.py run-hp "HPPH" -L 4 --dim 3`
- Visualization:
  - Plotly HTML: `--viz 3d` → `folding_3d_<prefix>.html`
  - Native viewer: `--viz native` (interactive window)
  - Snapshot: `--snap out/fold.png` (offscreen PNG via native backend)
- Batch: see `README.md` for CSV/JSON/LaTeX exports; outputs saved under `out/` when you specify paths.
- Saved solutions: `--save-solution out/solution.json` and later load via `view-solution`.
- Typical artifacts: PNGs and HTMLs in `out/` (e.g., `out/fold.png`, `folding_2d_<seq>.png`, `folding_3d_<seq>.html`).
