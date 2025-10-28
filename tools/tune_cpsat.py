from __future__ import annotations

from typing import Any, Dict, List, Optional, Mapping, Tuple
from dataclasses import dataclass
from pathlib import Path
import time

from ortools.sat.python import cp_model

import hp2d
import hp3d
import tune_params


@dataclass
class TuningInstance:
    sequences: List[str]
    L: int
    dim: int
    time_limit_s: float
    workers: int


def _build_model(seq: str, L: int, dim: int) -> Tuple[cp_model.CpModel, Dict[int, Dict[Any, cp_model.IntVar]], Dict[Any, cp_model.IntVar]]:
    model = cp_model.CpModel()
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
        raise ValueError("dim must be 2 or 3")
    hp2d.set_objective_max_contacts(model, C.values())
    return model, X, C


def evaluate_config_for_sequence(seq: str,
                                 L: int,
                                 dim: int,
                                 time_limit_s: float,
                                 workers: int,
                                 params: Optional[Mapping[str, Any]]) -> Dict[str, Any]:
    model, X, C = _build_model(seq, L, dim)
    solver = cp_model.CpSolver()
    hp2d.configure_solver(solver, time_limit_s=time_limit_s, workers=workers, params=params)
    t0 = time.time()
    status = solver.Solve(model)
    elapsed = time.time() - t0
    contacts = 0
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        contacts = sum(int(solver.Value(cij)) for cij in C.values())
    return {
        "status": status,
        "status_name": solver.StatusName(status),
        "contacts": int(contacts),
        "objective": int(solver.ObjectiveValue()) if status in (cp_model.OPTIMAL, cp_model.FEASIBLE) else 0,
        "time_s": elapsed,
        "wall_time_s": solver.WallTime(),
        "best_bound": solver.BestObjectiveBound() if status == cp_model.OPTIMAL else None,
    }


def aggregate_score(per_seq: List[Dict[str, Any]]) -> Tuple[int, float]:
    # Primary: maximize contacts (use negative for minimization-based tuners)
    total_contacts = sum(r.get("contacts", 0) for r in per_seq)
    # Secondary: minimize time for tie-breaks
    total_time = sum(r.get("time_s", 0.0) for r in per_seq)
    return total_contacts, total_time


def tune(instance: TuningInstance,
         budget_seconds: float,
         max_trials: Optional[int] = None,
         search_space: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    Run parameter tuning using cpsat-autotune, returning best params.
    
    Uses cpsat-autotune's tune_for_quality_within_timelimit to find optimal
    solver parameters for maximizing H-H contacts within the time limit.
    
    Args:
        instance: TuningInstance with sequences, grid size, dimension, etc.
        budget_seconds: Total budget (not directly used by cpsat-autotune)
        max_trials: Number of parameter configurations to try (default: 30)
        search_space: Currently unused (cpsat-autotune uses its own parameter space)
    
    Returns:
        Dictionary of optimized CP-SAT parameters, or empty dict if tuning fails
    """
    try:
        from cpsat_autotune.tune import tune_for_quality_within_timelimit
        print("Starting CP-SAT parameter tuning using cpsat-autotune...")
    except ImportError as e:
        print(f"cpsat-autotune not available: {e}")
        return {}

    # Use the first sequence as representative for tuning
    seq = instance.sequences[0]
    print(f"Tuning on sequence: {seq[:50]}{'...' if len(seq) > 50 else ''}")
    print(f"Grid: {instance.L}x{instance.L}{'x'+str(instance.L) if instance.dim == 3 else ''}, Time limit: {instance.time_limit_s}s")
    
    # Build model
    model, X, C = _build_model(seq, instance.L, instance.dim)
    
    # For maximization, worst case is 0 contacts
    obj_for_timeout = 0
    
    # Compute reasonable number of samples based on budget
    # Each trial runs n_samples_for_trial + n_samples_for_verification solves
    # Conservative estimate: each solve takes time_limit_s
    n_trials = max_trials if max_trials else 30
    
    try:
        best_params = tune_for_quality_within_timelimit(
            model=model,
            max_time_in_seconds=instance.time_limit_s,
            obj_for_timeout=obj_for_timeout,
            direction='maximize',
            n_samples_for_trial=3,
            n_samples_for_verification=5,
            n_trials=n_trials
        )
        
        if best_params:
            print(f"\n✓ Tuning complete! Found {len(best_params)} optimized parameters")
            return best_params
        else:
            print("\n⚠ Tuning found no improvements over defaults")
            return {}
            
    except Exception as e:
        print(f"✗ Tuning failed: {e}")
        import traceback
        traceback.print_exc()
        return {}


def save_best_params(params: Mapping[str, Any], dim: int, L: int, n: int, out_dir: str = "out/tuned") -> str:
    key = tune_params.build_problem_key(dim, L, n)
    path = Path(out_dir) / f"{key}.json"
    meta = {"dim": dim, "L": L, "n": n}
    return tune_params.save_params(params, str(path), meta)


