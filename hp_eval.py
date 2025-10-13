"""
HP Protein Folding Batch Evaluator

Provides batch evaluation, energy computation, and results export for HP folding problems.
"""

from typing import Dict, List, Tuple, Optional, Any, Callable
import time
import json
import csv
import logging
from pathlib import Path
from statistics import median

from ortools.sat.python import cp_model
from rich.console import Console
from rich.table import Table

import hp2d
import hp3d

logger = logging.getLogger(__name__)
console = Console()


def compute_energy_contacts(contacts_realized: List[Tuple[int, int]]) -> int:
    """
    Compute energy in ε_HH units from realized H-H contacts.
    
    Energy formula: E/ε_HH = -(# of non-consecutive H-H contacts)
    Each contact contributes -1 in ε_HH units.
    
    Args:
        contacts_realized: List of (i, j) pairs representing H-H contacts
        
    Returns:
        Energy in ε_HH units (negative integer)
    """
    return -len(contacts_realized)


def timeit(label: str) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
    """
    Decorator to measure wall-clock time of a function call.

    Prints a short log with the measured time and attaches an `elapsed_s`
    attribute (float seconds) to the wrapped function for introspection.

    Args:
        label: Human-friendly label for the timed section

    Returns:
        Decorated function that reports elapsed time after execution
    """
    def _decorator(func: Callable[..., Any]) -> Callable[..., Any]:
        def _wrapper(*args: Any, **kwargs: Any) -> Any:
            t0 = time.time()
            try:
                return func(*args, **kwargs)
            finally:
                elapsed = time.time() - t0
                # Attach dynamic attribute for callers that want to read it
                setattr(_wrapper, "elapsed_s", elapsed)
                console.print(f"[dim]{label} took {elapsed:.3f}s[/dim]")
        return _wrapper
    return _decorator


def solve_with_timer(solver: cp_model.CpSolver, model: cp_model.CpModel) -> Tuple[int, float]:
    """
    Solve the model with CP-SAT and measure the solve time using the timeit decorator.

    Args:
        solver: CP-SAT solver instance
        model: Model to solve

    Returns:
        (status, elapsed_seconds)
    """
    @timeit("CP-SAT Solve")
    def _solve(_solver: cp_model.CpSolver, _model: cp_model.CpModel) -> int:
        return _solver.Solve(_model)

    status = _solve(solver, model)
    elapsed_s = getattr(_solve, "elapsed_s", 0.0)
    return status, elapsed_s


def evaluate_sequence(
    seq: str,
    L: int,
    dim: int,
    time_limit_s: float,
    workers: int,
    seed: Optional[int] = None
) -> Dict[str, Any]:
    """
    Evaluate a single HP sequence and return detailed results.
    
    Args:
        seq: Protein sequence string (H/P letters)
        L: Grid dimension
        dim: Lattice dimension (2 or 3)
        time_limit_s: Solver time limit in seconds
        workers: Number of parallel workers
        seed: Random seed for solver (optional)
        
    Returns:
        Dictionary with keys:
        - sequence: Input sequence
        - n: Sequence length
        - dim: Lattice dimension
        - L: Grid size
        - status: Solver status string
        - contacts: Number of H-H contacts achieved
        - energy_epsHH: Energy in ε_HH units (= -contacts)
        - positions: Residue positions dict
        - objective_value: Objective value from solver
        - best_bound: Best bound if available
        - is_optimal: Whether solution is proven optimal
        - solve_time_s: Solve time in seconds
        - seed: Random seed used
    """
    logger.info(f"Evaluating {seq} (n={len(seq)}, dim={dim}, L={L})")
    console.print(f"[dim]Evaluating {seq[:20]}{'...' if len(seq) > 20 else ''} "
                  f"(n={len(seq)}, dim={dim}D, L={L})...[/dim]")
    
    model = cp_model.CpModel()
    
    # Build grid and helpers based on dimension
    if dim == 2:
        cells = hp2d.build_grid(L)
        nbrs = hp2d.neighbors_map(cells)
        start, up = hp2d.center_and_up(L)
        X = hp2d.declare_vars(model, n=len(seq), cells=cells)
        hp2d.add_exactly_one_cell_per_residue(model, X, cells)
        hp2d.add_at_most_one_residue_per_cell(model, X, cells)
        hp2d.add_chain_adjacency_allowed_pairs(model, X, nbrs, cells)
        hp2d.fix_start_and_orientation(model, X, start, up)
        C = hp2d.declare_hh_contact_vars(model, seq, X, nbrs, cells)
    elif dim == 3:
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
        raise ValueError(f"Invalid dimension: {dim}")
    
    # Set objective
    hp2d.set_objective_max_contacts(model, C.values())
    
    # Configure solver
    solver = cp_model.CpSolver()
    hp2d.configure_solver(solver, time_limit_s=time_limit_s, workers=workers)
    if seed is not None:
        solver.parameters.random_seed = seed
    
    # Solve
    status, solve_time = solve_with_timer(solver, model)
    
    # Extract results
    status_name = solver.StatusName(status)
    
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return {
            "sequence": seq,
            "n": len(seq),
            "dim": dim,
            "L": L,
            "status": status_name,
            "contacts": 0,
            "energy_epsHH": 0,
            "positions": {},
            "objective_value": 0,
            "best_bound": None,
            "is_optimal": False,
            "solve_time_s": solve_time,
            "seed": seed,
        }
    
    # Extract positions
    pos = {i: next(c for c in X[i].keys() if solver.Value(X[i][c]) == 1) 
           for i in range(len(seq))}
    
    # Extract realized contacts
    contacts_realized = [(i, j) for (i, j), var in C.items() if solver.Value(var) == 1]
    num_contacts = len(contacts_realized)
    energy = compute_energy_contacts(contacts_realized)
    
    # Get bounds
    obj_value = solver.ObjectiveValue()
    best_bound = solver.BestObjectiveBound() if status == cp_model.OPTIMAL else None
    
    result = {
        "sequence": seq,
        "n": len(seq),
        "dim": dim,
        "L": L,
        "status": status_name,
        "contacts": num_contacts,
        "energy_epsHH": energy,
        "positions": pos,
        "objective_value": int(obj_value),
        "best_bound": int(best_bound) if best_bound is not None else None,
        "is_optimal": status == cp_model.OPTIMAL,
        "solve_time_s": solve_time,
        "seed": seed,
        "contacts_realized": contacts_realized,  # For visualization
    }
    
    logger.info(f"Completed {seq}: E/ε_HH = {energy}, contacts = {num_contacts}, "
                f"status = {status_name}, time = {solve_time:.2f}s")
    
    return result


def evaluate_batch(
    sequences: List[str],
    L: int,
    dim: int,
    time_limit_s: float,
    workers: int,
    seeds: Optional[List[int]] = None,
    runs_per_seq: int = 1
) -> List[Dict[str, Any]]:
    """
    Evaluate multiple sequences in batch mode.
    
    Runs each sequence runs_per_seq times (with different seeds if provided),
    keeps the best (most negative energy) result per sequence.
    
    Args:
        sequences: List of protein sequences
        L: Grid dimension
        dim: Lattice dimension (2 or 3)
        time_limit_s: Solver time limit per run
        workers: Number of parallel workers
        seeds: Optional list of random seeds (one per run)
        runs_per_seq: Number of runs per sequence
        
    Returns:
        List of per-sequence summary dictionaries with keys:
        - sequence: Sequence string
        - n: Sequence length
        - dim: Lattice dimension
        - L: Grid size
        - best_energy_epsHH: Best (most negative) energy found
        - best_contacts: Contacts in best run
        - best_status: Status of best run
        - is_optimal: Whether best run was proven optimal
        - median_contacts: Median contacts across runs
        - best_run_solve_time_s: Solve time of best run
        - runs: Number of runs performed
    """
    console.print(f"\n[bold cyan]Batch Evaluation: {len(sequences)} sequences, "
                  f"{runs_per_seq} run(s) each[/bold cyan]")
    
    results = []
    
    for seq_idx, seq in enumerate(sequences):
        console.print(f"\n[bold]Sequence {seq_idx + 1}/{len(sequences)}:[/bold] {seq}")
        
        run_results = []
        for run in range(runs_per_seq):
            seed = seeds[run] if seeds and run < len(seeds) else None
            result = evaluate_sequence(seq, L, dim, time_limit_s, workers, seed)
            run_results.append(result)
            
            console.print(f"  Run {run + 1}: E/ε_HH = {result['energy_epsHH']}, "
                          f"contacts = {result['contacts']}, "
                          f"status = {result['status']}")
        
        # Find best run (most negative energy = most contacts)
        best_run = min(run_results, key=lambda r: r['energy_epsHH'])
        
        # Aggregate statistics
        all_contacts = [r['contacts'] for r in run_results]
        median_contacts_val = median(all_contacts) if all_contacts else 0
        
        summary = {
            "sequence": seq,
            "n": len(seq),
            "dim": dim,
            "L": L,
            "best_energy_epsHH": best_run['energy_epsHH'],
            "best_contacts": best_run['contacts'],
            "best_status": best_run['status'],
            "is_optimal": best_run['is_optimal'],
            "median_contacts": median_contacts_val,
            "best_run_solve_time_s": best_run['solve_time_s'],
            "runs": runs_per_seq,
            "best_positions": best_run.get('positions', {}),
            "best_contacts_realized": best_run.get('contacts_realized', []),
        }
        
        results.append(summary)
        
        console.print(f"[green]✓ Best: E/ε_HH = {summary['best_energy_epsHH']}, "
                      f"contacts = {summary['best_contacts']}, "
                      f"optimal = {summary['is_optimal']}[/green]")
    
    return results


def save_results_csv(results: List[Dict[str, Any]], path: str) -> str:
    """
    Save batch results to CSV file.
    
    Args:
        results: List of result dictionaries
        path: Output CSV file path
        
    Returns:
        Path to saved CSV file
    """
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    
    with open(path, 'w', newline='') as f:
        if not results:
            return path
        
        fieldnames = [
            'sequence', 'n', 'dim', 'L',
            'best_energy_epsHH', 'best_contacts',
            'best_status', 'is_optimal',
            'median_contacts', 'best_run_solve_time_s', 'runs'
        ]
        
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(results)
    
    console.print(f"[green]✓ CSV saved to: {path}[/green]")
    return path


def save_results_json(results: List[Dict[str, Any]], path: str) -> str:
    """
    Save batch results to JSON file.
    
    Args:
        results: List of result dictionaries
        path: Output JSON file path
        
    Returns:
        Path to saved JSON file
    """
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    
    # Convert tuples to lists for JSON serialization
    json_safe_results = []
    for r in results:
        r_copy = r.copy()
        # Convert position tuples to lists
        if 'best_positions' in r_copy:
            r_copy['best_positions'] = {
                str(k): list(v) for k, v in r_copy['best_positions'].items()
            }
        # Remove contacts_realized for cleaner JSON
        r_copy.pop('best_contacts_realized', None)
        json_safe_results.append(r_copy)
    
    with open(path, 'w') as f:
        json.dump(json_safe_results, f, indent=2)
    
    console.print(f"[green]✓ JSON saved to: {path}[/green]")
    return path


def latex_table(results: List[Dict[str, Any]]) -> str:
    """
    Generate a LaTeX tabular from batch results.
    
    Creates a table with columns:
    - Sequence (truncated if long)
    - n (length)
    - dim
    - L
    - Best E/ε_HH
    - Best contacts
    - Optimal?
    - Time (s)
    
    Args:
        results: List of result dictionaries
        
    Returns:
        LaTeX tabular code as string
    """
    lines = [
        r"\begin{tabular}{llcccccc}",
        r"\hline",
        r"Sequence & $n$ & Dim & $L$ & Best $E/\varepsilon_{HH}$ & Contacts & Optimal & Time (s) \\",
        r"\hline",
    ]
    
    for r in results:
        seq_short = r['sequence'][:12] + "..." if len(r['sequence']) > 15 else r['sequence']
        seq_short = seq_short.replace('_', r'\_')  # Escape underscores
        
        optimal_mark = r"\checkmark" if r['is_optimal'] else ""
        
        line = (f"{seq_short} & {r['n']} & {r['dim']} & {r['L']} & "
                f"{r['best_energy_epsHH']} & {r['best_contacts']} & "
                f"{optimal_mark} & {r['best_run_solve_time_s']:.2f} \\\\")
        lines.append(line)
    
    lines.extend([
        r"\hline",
        r"\end{tabular}",
    ])
    
    return "\n".join(lines)


def print_summary_table(results: List[Dict[str, Any]]) -> None:
    """
    Print a Rich table summarizing batch results to console.
    
    Args:
        results: List of result dictionaries
    """
    table = Table(title="Batch Evaluation Results", show_header=True, header_style="bold cyan")
    
    table.add_column("Sequence", style="dim", max_width=20)
    table.add_column("n", justify="right")
    table.add_column("Dim", justify="center")
    table.add_column("L", justify="right")
    table.add_column("Best E/ε_HH", justify="right", style="bold")
    table.add_column("Contacts", justify="right")
    table.add_column("Optimal", justify="center")
    table.add_column("Time (s)", justify="right")
    
    for r in results:
        seq_short = r['sequence'][:17] + "..." if len(r['sequence']) > 20 else r['sequence']
        optimal_str = "✓" if r['is_optimal'] else ""
        
        energy_style = "green" if r['best_energy_epsHH'] < 0 else "yellow"
        
        table.add_row(
            seq_short,
            str(r['n']),
            str(r['dim']),
            str(r['L']),
            f"[{energy_style}]{r['best_energy_epsHH']}[/{energy_style}]",
            str(r['best_contacts']),
            optimal_str,
            f"{r['best_run_solve_time_s']:.2f}"
        )
    
    console.print("\n")
    console.print(table)

