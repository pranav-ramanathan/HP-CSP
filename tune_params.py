from __future__ import annotations

from typing import Mapping, Any, Optional, Dict
from pathlib import Path
import json
from datetime import datetime

from ortools.sat.python import cp_model


def apply_solver_params(solver: cp_model.CpSolver, params: Optional[Mapping[str, Any]]) -> None:
    """
    Apply a dictionary of CP-SAT parameters to an existing solver instance.
    Safe no-op if params is None.
    """
    if not params:
        return
    for key, value in params.items():
        try:
            setattr(solver.parameters, key, value)
        except Exception:
            # Ignore unknown/invalid keys to be robust across OR-Tools versions
            continue


def build_problem_key(dim: int, L: int, n: int) -> str:
    return f"hp_{dim}d_L{L}_n{n}"


def save_params(params: Mapping[str, Any], out_path: str, meta: Optional[Dict[str, Any]] = None) -> str:
    """
    Save params (plus optional metadata) to JSON file, creating directories as needed.
    """
    payload = {
        "params": dict(params),
        "metadata": {
            "created_utc": datetime.utcnow().isoformat(timespec="seconds") + "Z",
            **(meta or {}),
        }
    }
    p = Path(out_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload, indent=2))
    return str(p)


def load_params(path: str) -> Optional[Dict[str, Any]]:
    """
    Load params dict from JSON file. Supports either {"params": {...}} or a flat dict.
    Returns None if file missing.
    """
    p = Path(path)
    if not p.exists():
        return None
    try:
        data = json.loads(p.read_text())
    except Exception:
        return None
    if isinstance(data, dict) and "params" in data and isinstance(data["params"], dict):
        return dict(data["params"])
    if isinstance(data, dict):
        return dict(data)
    return None


def resolve_params_for(dim: int, L: int, n: int, directory: str = "out/tuned") -> Optional[Dict[str, Any]]:
    """
    Attempt to find a tuned-params JSON for the given (dim, L, n) key.
    Currently exact match only: out/tuned/hp_{dim}d_L{L}_n{n}.json
    Returns the params dict or None.
    """
    key = build_problem_key(dim, L, n)
    candidate = Path(directory) / f"{key}.json"
    return load_params(str(candidate))


