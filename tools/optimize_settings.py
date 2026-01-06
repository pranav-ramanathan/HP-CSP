"""
Optimization Settings Finder for 3D HP Protein Folding

Orchestrates multiple strategies to find the best solver settings for a given sequence:
1. Grid size optimization - Find minimal L that allows good solutions
2. Multi-seed search - Run with different random seeds, take best
3. CP-SAT parameter tuning - Use cpsat-autotune for solver parameters
4. Symmetry variant testing - Compare different symmetry breaking strategies
"""

from __future__ import annotations

import json
import math
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from ortools.sat.python import cp_model
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TimeElapsedColumn

import hp3d
import hp2d
import tune_params

console = Console()


@dataclass
class RunResult:
    """Result from a single solver run."""
    contacts: int
    energy: int
    L: int
    seed: Optional[int]
    time_s: float
    status: str
    is_optimal: bool
    positions: Optional[Dict[int, Tuple[int, int, int]]] = None


@dataclass
class OptimizationResult:
    """Aggregated result from optimization process."""
    sequence: str
    expanded_length: int
    dim: int
    best_energy: int
    best_contacts: int
    best_L: int
    best_seed: Optional[int]
    best_params: Dict[str, Any]
    optimization_time_s: float
    all_runs: List[Dict[str, Any]] = field(default_factory=list)
    recommendations: Dict[str, Any] = field(default_factory=dict)
    best_positions: Optional[Dict[int, Tuple[int, int, int]]] = None


def compute_min_grid_size(n: int) -> int:
    """Compute minimum grid size L = ceil(2 * n^(1/3))."""
    return int(math.ceil(2 * (n ** (1/3))))


def _build_and_solve(
    seq: str,
    L: int,
    dim: int,
    time_limit_s: float,
    workers: int,
    seed: Optional[int] = None,
    params: Optional[Dict[str, Any]] = None
) -> RunResult:
    """
    Build model and solve, returning a RunResult.

    This is a low-level function used by the optimization strategies.
    """
    model = cp_model.CpModel()

    if dim == 3:
        cells = hp3d.build_grid_3d(L)
        nbrs = hp3d.neighbors_map_3d(cells)
        start, up = hp3d.center_and_up_3d(L)
        X = hp3d.declare_vars_3d(model, n=len(seq), cells=cells)
        hp3d.add_exactly_one_cell_per_residue_3d(model, X, cells)
        hp3d.add_at_most_one_residue_per_cell_3d(model, X, cells)
        hp3d.add_chain_adjacency_allowed_pairs_3d(model, X, nbrs, cells)
        hp3d.fix_start_and_orientation_3d(model, X, start, up)
        C = hp3d.declare_hh_contact_vars_3d(model, seq, X, nbrs, cells)
    else:
        cells = hp2d.build_grid(L)
        nbrs = hp2d.neighbors_map(cells)
        start, up = hp2d.center_and_up(L)
        X = hp2d.declare_vars(model, n=len(seq), cells=cells)
        hp2d.add_exactly_one_cell_per_residue(model, X, cells)
        hp2d.add_at_most_one_residue_per_cell(model, X, cells)
        hp2d.add_chain_adjacency_allowed_pairs(model, X, nbrs, cells)
        hp2d.fix_start_and_orientation(model, X, start, up)
        C = hp2d.declare_hh_contact_vars(model, seq, X, nbrs, cells)

    hp2d.set_objective_max_contacts(model, C.values())

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit_s
    solver.parameters.num_search_workers = workers
    if seed is not None:
        solver.parameters.random_seed = seed
    tune_params.apply_solver_params(solver, params)

    t0 = time.time()
    status = solver.Solve(model)
    elapsed = time.time() - t0

    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        contacts = sum(int(solver.Value(cij)) for cij in C.values())
        positions = {i: next(c for c in X[i].keys() if solver.Value(X[i][c]) == 1)
                     for i in range(len(seq))}
        return RunResult(
            contacts=contacts,
            energy=-contacts,
            L=L,
            seed=seed,
            time_s=elapsed,
            status=solver.StatusName(status),
            is_optimal=(status == cp_model.OPTIMAL),
            positions=positions
        )
    else:
        return RunResult(
            contacts=0,
            energy=0,
            L=L,
            seed=seed,
            time_s=elapsed,
            status=solver.StatusName(status),
            is_optimal=False,
            positions=None
        )


def find_optimal_grid_size(
    seq: str,
    dim: int = 3,
    time_per_L: float = 30.0,
    workers: int = 8,
    max_L_delta: int = 3
) -> Tuple[int, List[Dict[str, Any]]]:
    """
    Find the optimal grid size for a sequence.

    Tries L, L+1, L+2 (and optionally L+3) where L is the minimum computed size.
    Returns the smallest L that achieves the best contact count.

    Args:
        seq: HP sequence string
        dim: Dimension (3 for 3D)
        time_per_L: Time limit per grid size test
        workers: Number of parallel workers
        max_L_delta: Maximum L values to try above minimum

    Returns:
        (best_L, results_log) tuple
    """
    n = len(seq)
    min_L = compute_min_grid_size(n)

    results = []
    best_contacts = -1
    best_L = min_L

    console.print(f"[dim]Grid size search: n={n}, min_L={min_L}, testing L={min_L} to L={min_L + max_L_delta - 1}[/dim]")

    for delta in range(max_L_delta):
        L = min_L + delta
        result = _build_and_solve(seq, L, dim, time_per_L, workers)

        log_entry = {
            "L": L,
            "contacts": result.contacts,
            "energy": result.energy,
            "time_s": result.time_s,
            "status": result.status,
            "is_optimal": result.is_optimal
        }
        results.append(log_entry)

        console.print(f"  L={L}: contacts={result.contacts}, energy={result.energy}, "
                      f"status={result.status}, time={result.time_s:.1f}s")

        # Track best (prefer smaller L on tie)
        if result.contacts > best_contacts:
            best_contacts = result.contacts
            best_L = L

        # Early stop if optimal found
        if result.is_optimal:
            console.print(f"[green]  Optimal found at L={L}, stopping grid search[/green]")
            break

    return best_L, results


def multi_seed_search(
    seq: str,
    L: int,
    dim: int = 3,
    n_seeds: int = 5,
    time_per_seed: float = 30.0,
    workers: int = 8,
    params: Optional[Dict[str, Any]] = None
) -> Tuple[RunResult, List[int], List[Dict[str, Any]]]:
    """
    Run solver with multiple random seeds and return the best result.

    Args:
        seq: HP sequence string
        L: Grid size
        dim: Dimension (3 for 3D)
        n_seeds: Number of seeds to try
        time_per_seed: Time limit per seed
        workers: Number of parallel workers
        params: Optional CP-SAT parameters

    Returns:
        (best_result, top_seeds, results_log) tuple
    """
    # Use diverse seeds
    seed_list = [42, 17, 99, 1337, 2024, 7, 123, 456, 789, 314][:n_seeds]

    results = []
    best_result = None

    console.print(f"[dim]Multi-seed search: L={L}, {n_seeds} seeds, {time_per_seed}s each[/dim]")

    for seed in seed_list:
        result = _build_and_solve(seq, L, dim, time_per_seed, workers, seed=seed, params=params)

        log_entry = {
            "seed": seed,
            "contacts": result.contacts,
            "energy": result.energy,
            "time_s": result.time_s,
            "status": result.status,
            "is_optimal": result.is_optimal
        }
        results.append(log_entry)

        console.print(f"  seed={seed}: contacts={result.contacts}, energy={result.energy}, "
                      f"status={result.status}, time={result.time_s:.1f}s")

        if best_result is None or result.contacts > best_result.contacts:
            best_result = result

        # Early stop if optimal
        if result.is_optimal:
            console.print(f"[green]  Optimal found with seed={seed}, stopping seed search[/green]")
            break

    # Sort seeds by contacts achieved (descending)
    sorted_results = sorted(results, key=lambda r: r["contacts"], reverse=True)
    top_seeds = [r["seed"] for r in sorted_results[:3]]

    return best_result, top_seeds, results


def quick_tune_params(
    seq: str,
    L: int,
    dim: int = 3,
    budget_s: float = 120.0,
    time_limit_s: float = 30.0,
    workers: int = 8,
    n_trials: int = 10
) -> Dict[str, Any]:
    """
    Quick parameter tuning using cpsat-autotune.

    Uses fewer trials than full tuning for faster results.

    Args:
        seq: HP sequence string
        L: Grid size
        dim: Dimension
        budget_s: Total tuning budget (informational, not directly used)
        time_limit_s: Time limit per solve during tuning
        workers: Number of workers
        n_trials: Number of parameter configurations to try

    Returns:
        Dictionary of tuned parameters (empty if tuning fails/unavailable)
    """
    try:
        from cpsat_autotune.tune import tune_for_quality_within_timelimit
    except ImportError:
        console.print("[yellow]cpsat-autotune not available, skipping parameter tuning[/yellow]")
        return {}

    console.print(f"[dim]Quick param tuning: {n_trials} trials, {time_limit_s}s per solve[/dim]")

    # Build model
    model = cp_model.CpModel()
    if dim == 3:
        cells = hp3d.build_grid_3d(L)
        nbrs = hp3d.neighbors_map_3d(cells)
        start, up = hp3d.center_and_up_3d(L)
        X = hp3d.declare_vars_3d(model, n=len(seq), cells=cells)
        hp3d.add_exactly_one_cell_per_residue_3d(model, X, cells)
        hp3d.add_at_most_one_residue_per_cell_3d(model, X, cells)
        hp3d.add_chain_adjacency_allowed_pairs_3d(model, X, nbrs, cells)
        hp3d.fix_start_and_orientation_3d(model, X, start, up)
        C = hp3d.declare_hh_contact_vars_3d(model, seq, X, nbrs, cells)
    else:
        cells = hp2d.build_grid(L)
        nbrs = hp2d.neighbors_map(cells)
        start, up = hp2d.center_and_up(L)
        X = hp2d.declare_vars(model, n=len(seq), cells=cells)
        hp2d.add_exactly_one_cell_per_residue(model, X, cells)
        hp2d.add_at_most_one_residue_per_cell(model, X, cells)
        hp2d.add_chain_adjacency_allowed_pairs(model, X, nbrs, cells)
        hp2d.fix_start_and_orientation(model, X, start, up)
        C = hp2d.declare_hh_contact_vars(model, seq, X, nbrs, cells)

    hp2d.set_objective_max_contacts(model, C.values())

    try:
        best_params = tune_for_quality_within_timelimit(
            model=model,
            max_time_in_seconds=time_limit_s,
            obj_for_timeout=0,
            direction='maximize',
            n_samples_for_trial=2,  # Reduced for speed
            n_samples_for_verification=3,  # Reduced for speed
            n_trials=n_trials
        )

        if best_params:
            console.print(f"[green]  Found {len(best_params)} tuned parameters[/green]")
            return best_params
        else:
            console.print("[yellow]  No improvements found over defaults[/yellow]")
            return {}

    except Exception as e:
        console.print(f"[red]  Tuning failed: {e}[/red]")
        return {}


def optimize_for_sequence(
    seq: str,
    dim: int = 3,
    time_budget_s: float = 300.0,
    workers: int = 8,
    skip_grid_search: bool = False,
    skip_tuning: bool = False,
    n_seeds: int = 5
) -> OptimizationResult:
    """
    Main optimization orchestrator.

    Runs multiple strategies to find the best settings for a sequence:
    1. Grid size optimization (20% of budget)
    2. Multi-seed search (40% of budget)
    3. Parameter tuning (40% of budget)

    Args:
        seq: HP sequence string
        dim: Dimension (3 for 3D)
        time_budget_s: Total time budget in seconds
        workers: Number of parallel workers
        skip_grid_search: Skip grid size optimization
        skip_tuning: Skip CP-SAT parameter tuning
        n_seeds: Number of random seeds to try

    Returns:
        OptimizationResult with best settings and recommendations
    """
    start_time = time.time()
    n = len(seq)

    console.print(f"\n[bold cyan]Optimizing settings for sequence (n={n}, dim={dim}D)[/bold cyan]")
    console.print(f"[dim]Total budget: {time_budget_s}s, workers: {workers}[/dim]")

    all_runs = []
    best_overall: Optional[RunResult] = None
    best_params: Dict[str, Any] = {}
    recommended_seeds: List[int] = []

    # Phase 1: Grid size search (20% of budget)
    if not skip_grid_search:
        grid_budget = time_budget_s * 0.2
        time_per_L = grid_budget / 3  # Test 3 grid sizes

        console.print(f"\n[bold]Phase 1: Grid Size Search[/bold] (budget: {grid_budget:.0f}s)")
        best_L, grid_results = find_optimal_grid_size(
            seq, dim, time_per_L=time_per_L, workers=workers, max_L_delta=3
        )

        for r in grid_results:
            r["phase"] = "grid_search"
            all_runs.append(r)

        # Get best result from grid search
        best_grid = max(grid_results, key=lambda r: r["contacts"])
        if best_grid["contacts"] > 0:
            best_overall = RunResult(
                contacts=best_grid["contacts"],
                energy=best_grid["energy"],
                L=best_grid["L"],
                seed=None,
                time_s=best_grid["time_s"],
                status=best_grid["status"],
                is_optimal=best_grid["is_optimal"]
            )
    else:
        best_L = compute_min_grid_size(n)
        console.print(f"[dim]Skipping grid search, using L={best_L}[/dim]")

    # Phase 2: Multi-seed search (40% of budget)
    seed_budget = time_budget_s * 0.4
    time_per_seed = seed_budget / n_seeds

    console.print(f"\n[bold]Phase 2: Multi-Seed Search[/bold] (budget: {seed_budget:.0f}s)")
    seed_result, top_seeds, seed_results = multi_seed_search(
        seq, best_L, dim, n_seeds=n_seeds, time_per_seed=time_per_seed, workers=workers
    )

    for r in seed_results:
        r["phase"] = "seed_search"
        r["L"] = best_L
        all_runs.append(r)

    recommended_seeds = top_seeds

    if seed_result and (best_overall is None or seed_result.contacts > best_overall.contacts):
        best_overall = seed_result

    # Phase 3: Parameter tuning (40% of budget)
    if not skip_tuning:
        tune_budget = time_budget_s * 0.4
        # Use shorter time per solve during tuning
        tune_time_per_solve = min(30.0, tune_budget / 20)  # At least 20 solves possible
        n_trials = max(5, int(tune_budget / (tune_time_per_solve * 5)))  # 5 samples per trial

        console.print(f"\n[bold]Phase 3: Parameter Tuning[/bold] (budget: {tune_budget:.0f}s)")
        best_params = quick_tune_params(
            seq, best_L, dim,
            budget_s=tune_budget,
            time_limit_s=tune_time_per_solve,
            workers=workers,
            n_trials=n_trials
        )

        # Verify tuned params improve results
        if best_params:
            console.print("[dim]Verifying tuned parameters...[/dim]")
            verify_result = _build_and_solve(
                seq, best_L, dim,
                time_limit_s=time_per_seed * 2,  # Give more time for verification
                workers=workers,
                seed=recommended_seeds[0] if recommended_seeds else 42,
                params=best_params
            )

            verify_log = {
                "phase": "param_verify",
                "L": best_L,
                "seed": recommended_seeds[0] if recommended_seeds else 42,
                "contacts": verify_result.contacts,
                "energy": verify_result.energy,
                "time_s": verify_result.time_s,
                "status": verify_result.status,
                "is_optimal": verify_result.is_optimal,
                "params_used": True
            }
            all_runs.append(verify_log)

            console.print(f"  With tuned params: contacts={verify_result.contacts}, "
                          f"energy={verify_result.energy}")

            if verify_result.contacts > best_overall.contacts:
                best_overall = verify_result
                console.print("[green]  Tuned params improved result![/green]")
            else:
                console.print("[yellow]  Tuned params did not improve result[/yellow]")

    elapsed = time.time() - start_time

    # Build recommendations
    recommendations = {
        "grid_size": best_L,
        "recommended_seeds": recommended_seeds,
        "use_tuned_params": len(best_params) > 0,
        "params": best_params if best_params else None,
    }

    result = OptimizationResult(
        sequence=seq,
        expanded_length=n,
        dim=dim,
        best_energy=best_overall.energy if best_overall else 0,
        best_contacts=best_overall.contacts if best_overall else 0,
        best_L=best_L,
        best_seed=best_overall.seed if best_overall else None,
        best_params=best_params,
        optimization_time_s=elapsed,
        all_runs=all_runs,
        recommendations=recommendations,
        best_positions=best_overall.positions if best_overall else None
    )

    # Summary
    console.print(f"\n[bold green]Optimization Complete[/bold green]")
    console.print(f"  Best energy: {result.best_energy} (contacts: {result.best_contacts})")
    console.print(f"  Best L: {result.best_L}")
    console.print(f"  Best seed: {result.best_seed}")
    console.print(f"  Tuned params: {len(best_params)} parameters")
    console.print(f"  Total time: {elapsed:.1f}s")

    return result


def save_optimization_result(result: OptimizationResult, path: str) -> str:
    """
    Save optimization result to JSON file.

    Args:
        result: OptimizationResult to save
        path: Output file path

    Returns:
        Path to saved file
    """
    Path(path).parent.mkdir(parents=True, exist_ok=True)

    # Convert positions to JSON-serializable format
    positions_json = None
    if result.best_positions:
        positions_json = {str(k): list(v) for k, v in result.best_positions.items()}

    payload = {
        "sequence": result.sequence,
        "expanded_length": result.expanded_length,
        "dim": result.dim,
        "optimization_time_s": result.optimization_time_s,
        "best_result": {
            "energy": result.best_energy,
            "contacts": result.best_contacts,
            "L": result.best_L,
            "seed": result.best_seed,
            "positions": positions_json
        },
        "recommendations": result.recommendations,
        "exploration_log": result.all_runs
    }

    with open(path, 'w') as f:
        json.dump(payload, f, indent=2)

    console.print(f"[green]Saved optimization result to: {path}[/green]")
    return path
