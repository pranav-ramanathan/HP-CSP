from typing import Dict, List, Tuple, Iterable, Sequence, Optional
import logging
from ortools.sat.python import cp_model
from rich.console import Console
import typer

logger = logging.getLogger(__name__)
console = Console()
app = typer.Typer()

Cell = Tuple[int, int]
Residue = int

# Grid helpers
def build_grid(L: int) -> List[Cell]:
    """
    Return all cells of an LxL grid as (x, y), in [0..L-1]
    """
    

    logger.info(f"Building grid for {L}x{L}")
    console.print(f"[dim]Building grid for {L}*{L}...[/dim]")

    return [(x, y) for y in range(L) for x in range(L)]


def neighbors_map(cells: Sequence[Cell]) -> Dict[Cell, List[Cell]]:
    """
    Return a 4-neighborhood map: for each cell u, list its Manhattan neighbors v.
    """

    logger.info(f"Building neighborhood map for {len(cells)} cells")
    console.print(f"[dim]Building neighborhood map for {len(cells)} cells...[/dim]")

    cell_set = set(cells)
    offsets = [(0, 1), (0, -1), (1, 0), (-1, 0)]
    return {u: [(u[0] + dx, u[1] + dy) for dx, dy in offsets if (u[0] + dx, u[1] + dy) in cell_set] for u in cells}

def center_and_up(L: int) -> Tuple[Cell, Cell]:
    """
    Return (center_cell, up_cell) for symmetry fixing.
    'up_cell' should be the cell directly above center if it exists;
    otherwise choose a deterministic fallback neighbor (e.g., below).
    """

    logger.info(f"Finding center and up cell for {L}x{L} grid")
    console.print(f"[dim]Finding center and up cell for {L}x{L} grid...[/dim]")

    center_x, center_y = L // 2, L // 2
    center_cell = (center_x, center_y)
    # "up" means decreasing y (standard grid convention)
    up_cell = (center_x, center_y - 1) if center_y > 0 else (center_x, center_y + 1)

    return center_cell, up_cell

# Variable declaration

def declare_vars(model: cp_model.CpModel, n: int, cells: Sequence[Cell]) -> Dict[Residue, Dict[Cell, cp_model.IntVar]]:
    """
    Create one-hot Bool vars X[i][c] ∈ {0,1} meaning residue i is at cell c.
    Use: model.NewBoolVar(name)
    Return: Dict[Residue, Dict[Cell, cp_model.IntVar]]
    """

    logger.info(f"Declaring {n} one-hot Bool vars for {len(cells)} cells")
    console.print(f"[dim]Declaring {n} one-hot Bool vars for {len(cells)} cells...[/dim]")

    return {i: {c: model.NewBoolVar(f"X_{i}_{c}") for c in cells} for i in range(n)}

# Core constraints

def add_exactly_one_cell_per_residue(model: cp_model.CpModel, X: Dict[Residue, Dict[Cell, cp_model.IntVar]], cells: Sequence[Cell]) -> None:
    """
    For each residue i: sum_c X[i][c] == 1
    Use: model.Add(sum(...) == 1)
    """

    logger.info(f"Adding exactly one cell per residue for {len(X)} residues")
    console.print(f"[dim]Adding exactly one cell per residue for {len(X)} residues...[/dim]")

    for i in range(len(X)):
        model.Add(sum(X[i][c] for c in cells) == 1)
        console.print(f"[dim]Added exactly one cell per residue for residue {i}...[/dim]")


def add_at_most_one_residue_per_cell(model: cp_model.CpModel,
                                     X: Dict[Residue, Dict[Cell, cp_model.IntVar]],
                                     cells: Sequence[Cell]) -> None:
    """
    For each cell c: sum_i X[i][c] <= 1
    Use: model.Add(sum(...) <= 1)
    """
    logger.info(f"Adding at most one residue per cell for {len(X)} residues and {len(cells)} cells")
    console.print(f"[dim]Adding at most one residue per cell for {len(X)} residues and {len(cells)} cells...[/dim]")
    
    for c in cells:
        model.Add(sum(X[i][c] for i in range(len(X))) <= 1)
        console.print(f"[dim]Added at most one residue per cell for cell {c}...[/dim]")


def add_chain_adjacency_allowed_pairs(model: cp_model.CpModel,
                                      X: Dict[Residue, Dict[Cell, cp_model.IntVar]],
                                      neighbors: Dict[Cell, List[Cell]],
                                      cells: Sequence[Cell]) -> None:
    """
    Enforce for each consecutive pair (i, i+1) that chosen cells are neighbors.
    Two common patterns (choose one and stay consistent):
      A) Implication form:
         For each cell u: (X[i][u] == 1) => sum_{v in N(u)} X[i+1][v] == 1
         Use: OnlyEnforceIf on BoolVars.
      B) Table form (if you create cell choice IntVars):
         Use: model.AddAllowedAssignments([cell_i, cell_j], allowed_pairs)
    Here assume one-hot; implement pattern A.
    """
    logger.info(f"Adding chain adjacency constraints for {len(X)} residues and {len(neighbors)} neighbors")
    console.print(f"[dim]Adding chain adjacency constraints for {len(X)} residues and {len(neighbors)} neighbors...[/dim]")
    
    for i in range(len(X)-1):
        for c in cells:
            model.Add(sum(X[i+1][v] for v in neighbors[c]) == 1).OnlyEnforceIf(X[i][c])
        console.print(f"[dim]Added chain adjacency constraints for residue {i}...[/dim]")


# ---------- Symmetry breaking ----------

def fix_start_and_orientation(model: cp_model.CpModel,
                              X: Dict[Residue, Dict[Cell, cp_model.IntVar]],
                              start_cell: Cell,
                              up_cell: Cell) -> None:
    """
    Fix residue 0 at start_cell and residue 1 at up_cell:
        model.Add(X[0][start_cell] == 1)
        model.Add(X[1][up_cell] == 1)
    Also consider forbidding other cells for residues 0 and 1 for tightness.
    """
    logger.info(f"Fixing residue 0 at {start_cell} and residue 1 at {up_cell}")
    console.print(f"[dim]Fixing residue 0 at {start_cell} and residue 1 at {up_cell}...[/dim]")
    
    # Fix residue 0 at start_cell, zero all others
    for c in X[0]:
        model.Add(X[0][c] == int(c == start_cell))
    
    # Fix residue 1 at up_cell, zero all others
    for c in X[1]:
        model.Add(X[1][c] == int(c == up_cell))
    
    console.print(f"[dim]Fixed residue 0 at {start_cell} and residue 1 at {up_cell}...[/dim]")


# ---------- Contacts (objective features) ----------

def declare_hh_contact_vars(model: cp_model.CpModel,
                            seq: str,
                            X: Dict[Residue, Dict[Cell, cp_model.IntVar]],
                            neighbors: Dict[Cell, List[Cell]],
                            cells: Sequence[Cell]) -> Dict[Tuple[Residue, Residue], cp_model.IntVar]:
    """
    Edge-based contact formulation (memory-lean).

    Builds contact variables per undirected grid edge {u,v} and excludes chain edges via
    directed chain usage indicators B[i,(u,v)]. Returns a dict keyed by normalized
    undirected edges (min(u,v), max(u,v)).
    """
    n = len(seq)
    # Build undirected edge set and directed edges per u
    undirected_edges: List[Tuple[Cell, Cell]] = []
    seen = set()
    for u, nbrs in neighbors.items():
        for v in nbrs:
            key = (u, v) if u < v else (v, u)
            if key not in seen:
                seen.add(key)
                undirected_edges.append(key)

    # Chain edge usage indicators B[i,(u,v)]
    B: Dict[Tuple[int, Cell, Cell], cp_model.IntVar] = {}
    for i in range(n - 1):
        for u in cells:
            for v in neighbors[u]:
                b = model.NewBoolVar(f"B_{i}_{u}_{v}")
                B[(i, u, v)] = b
                model.Add(b <= X[i][u])
                model.Add(b <= X[i+1][v])
            # Tie to chain adjacency: exactly one outgoing if X[i,u]==1
            model.Add(sum(B[(i, u, v)] for v in neighbors[u]) == X[i][u])

    # Occupancy per cell and hydrophobic occupancy flags
    Occ: Dict[Cell, cp_model.IntVar] = {}
    H_at: Dict[Cell, cp_model.IntVar] = {}
    for u in cells:
        occ = model.NewIntVar(0, 1, f"Occ_{u}")
        model.Add(occ == sum(X[i][u] for i in range(n)))
        Occ[u] = occ
        h_occ = model.NewIntVar(0, 1, f"H_at_{u}")
        model.Add(h_occ == sum(X[i][u] for i in range(n) if seq[i] == 'H'))
        H_at[u] = h_occ

    # One contact variable per undirected edge {u,v}
    C_edge: Dict[Tuple[Cell, Cell], cp_model.IntVar] = {}
    for (u, v) in undirected_edges:
        c = model.NewBoolVar(f"Cedge_{u}_{v}")
        # Contacts only if both occupied and both H
        model.Add(c <= Occ[u])
        model.Add(c <= Occ[v])
        model.Add(c <= H_at[u])
        model.Add(c <= H_at[v])
        # Exclude chain edges (both directions)
        for i in range(n - 1):
            model.Add(c <= 1 - B[(i, u, v)])
            model.Add(c <= 1 - B[(i, v, u)])
        key = (u, v) if u < v else (v, u)
        C_edge[key] = c

    logger.info(f"Edge-based contacts: {len(C_edge)} edges, B vars: ~{(n-1)*sum(len(neighbors[u]) for u in cells)}")
    console.print(f"[dim]Edge-based contacts active: {len(C_edge)} edges[/dim]")
    return C_edge


def set_objective_max_contacts(model: cp_model.CpModel,
                               contacts: Iterable[cp_model.IntVar]) -> None:
    """
    model.Maximize(sum(contacts))
    """
    logger.info(f"Setting objective to maximize {len(contacts)} contact variables")
    console.print(f"[dim]Setting objective to maximize {len(contacts)} contact variables...[/dim]")
    
    model.Maximize(sum(contacts))
    console.print(f"[dim]Objective set to maximize {len(contacts)} contact variables...[/dim]")


# ---------- Solve & extract ----------

def configure_solver(solver: cp_model.CpSolver,
                     time_limit_s: float = 10.0,
                     workers: int = 8) -> None:
    """
    Set standard knobs, e.g.:
      solver.parameters.max_time_in_seconds = time_limit_s
      solver.parameters.num_search_workers = workers
    """
    logger.info(f"Setting time limit to {time_limit_s} seconds and {workers} workers")
    console.print(f"[dim]Setting time limit to {time_limit_s} seconds and {workers} workers...[/dim]")
    
    solver.parameters.max_time_in_seconds = time_limit_s
    solver.parameters.num_search_workers = workers
    console.print(f"[dim]Time limit and workers set...[/dim]")


def extract_positions(solver: cp_model.CpSolver,
                      X: Dict[Residue, Dict[Cell, cp_model.IntVar]]) -> Dict[Residue, Cell]:
    """
    Read back the chosen cell for each residue:
      find the unique c with solver.Value(X[i][c]) == 1.
    """
    logger.info(f"Extracting positions for {len(X)} residues")
    console.print(f"[dim]Extracting positions for {len(X)} residues...[/dim]")
    
    return {i: next(c for c in X[i].keys() if solver.Value(X[i][c]) == 1) for i in range(len(X))}


# ---------- Test Driver ----------

@app.command()
def run_hp(
    seq: str = typer.Argument(..., help="Sequence string (e.g., 'HPPH')"),
    L: int = typer.Option(4, "--grid-size", "-L", help="Grid size (L×L)"),
    time_limit: float = typer.Option(5.0, "--time", "-t", help="Solver time limit in seconds"),
    workers: int = typer.Option(8, "--workers", "-w", help="Number of parallel workers")
):
    """
    HP protein folding solver using OR-Tools CP-SAT.
    
    Finds the optimal folding of a protein sequence on a 2D lattice
    by maximizing H-H contacts.
    """
    assert L >= 3, "Use L>=3 so symmetry 'up' is valid"
    model = cp_model.CpModel()

    cells = build_grid(L)
    nbrs = neighbors_map(cells)
    start, up = center_and_up(L)

    # Vars
    X = declare_vars(model, n=len(seq), cells=cells)

    # Core constraints
    add_exactly_one_cell_per_residue(model, X, cells)
    add_at_most_one_residue_per_cell(model, X, cells)
    add_chain_adjacency_allowed_pairs(model, X, nbrs, cells)

    # Symmetry
    fix_start_and_orientation(model, X, start, up)

    # Contacts + objective
    C = declare_hh_contact_vars(model, seq, X, nbrs, cells)
    set_objective_max_contacts(model, C.values())

    # Solve
    solver = cp_model.CpSolver()
    configure_solver(solver, time_limit_s=time_limit, workers=workers)
    status = solver.Solve(model)

    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        console.print("[red]No solution found.[/red]")
        return

    pos = extract_positions(solver, X)
    contacts = sum(int(solver.Value(cij)) for cij in C.values())

    console.print(f"[green]Status: {solver.StatusName(status)}[/green]")
    console.print(f"[cyan]Contacts: {contacts}[/cyan]")
    console.print(f"[dim]Positions: {pos}[/dim]")

    # Pretty print grid
    grid = [["." for _ in range(L)] for __ in range(L)]
    for i, (x, y) in pos.items():
        grid[y][x] = seq[i]
    
    console.print("\n[bold]Grid:[/bold]")
    for row in grid:
        console.print(" ".join(row))


if __name__ == "__main__":
    app()