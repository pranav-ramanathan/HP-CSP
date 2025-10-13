"""
3D HP Lattice Protein Folding Module

Provides 3D grid helpers, constraint builders, and interactive visualization
for the HP protein folding problem on a cubic lattice.
"""

from typing import Dict, List, Tuple, Sequence, Iterable, Optional
import logging
from ortools.sat.python import cp_model
from rich.console import Console

logger = logging.getLogger(__name__)
console = Console()

Cell3D = Tuple[int, int, int]
Residue = int


# ---------- 3D Grid Helpers ----------

def build_grid_3d(L: int) -> List[Cell3D]:
    """
    Return all cells of an L*L*L grid as (x, y, z), in [0..L-1].
    
    Args:
        L: Grid dimension (creates L*L*L cube)
        
    Returns:
        List of all (x, y, z) cell coordinates
    """
    logger.info(f"Building 3D grid for {L}*{L}*{L}")
    console.print(f"[dim]Building 3D grid for {L}*{L}*{L}...[/dim]")
    
    return [(x, y, z) for z in range(L) for y in range(L) for x in range(L)]


def neighbors_map_3d(cells: Sequence[Cell3D]) -> Dict[Cell3D, List[Cell3D]]:
    """
    Return a 6-neighborhood map for 3D: for each cell u, list its Manhattan neighbors v.
    Uses offsets (±1,0,0), (0,±1,0), (0,0,±1).
    
    Args:
        cells: List of all cells in the grid
        
    Returns:
        Dictionary mapping each cell to its list of neighboring cells
    """
    logger.info(f"Building 3D neighborhood map for {len(cells)} cells")
    console.print(f"[dim]Building 3D neighborhood map for {len(cells)} cells...[/dim]")
    
    cell_set = set(cells)
    offsets = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    return {
        u: [(u[0] + dx, u[1] + dy, u[2] + dz) 
            for dx, dy, dz in offsets 
            if (u[0] + dx, u[1] + dy, u[2] + dz) in cell_set]
        for u in cells
    }


def center_and_up_3d(L: int) -> Tuple[Cell3D, Cell3D]:
    """
    Return (center_cell, up_cell) for 3D symmetry fixing.
    Center is at (L//2, L//2, L//2).
    'up_cell' is along the z-axis: (cx, cy, cz-1) if in-bounds, else (cx, cy, cz+1).
    
    Args:
        L: Grid dimension
        
    Returns:
        Tuple of (center_cell, up_cell) coordinates
    """
    logger.info(f"Finding center and up cell for {L}*{L}*{L} grid")
    console.print(f"[dim]Finding center and up cell for {L}*{L}*{L} grid...[/dim]")
    
    cx, cy, cz = L // 2, L // 2, L // 2
    center_cell = (cx, cy, cz)
    # "up" means decreasing z (upward along z-axis)
    up_cell = (cx, cy, cz - 1) if cz > 0 else (cx, cy, cz + 1)
    
    return center_cell, up_cell


# ---------- 3D Variable Declaration ----------

def declare_vars_3d(model: cp_model.CpModel, n: int, cells: Sequence[Cell3D]) -> Dict[Residue, Dict[Cell3D, cp_model.IntVar]]:
    """
    Create one-hot Bool vars X[i][c] ∈ {0,1} meaning residue i is at cell c.
    
    Args:
        model: CP-SAT model
        n: Number of residues
        cells: List of all cells in the grid
        
    Returns:
        Dictionary mapping residue -> cell -> BoolVar
    """
    logger.info(f"Declaring {n} one-hot Bool vars for {len(cells)} cells (3D)")
    console.print(f"[dim]Declaring {n} one-hot Bool vars for {len(cells)} cells (3D)...[/dim]")
    
    return {i: {c: model.NewBoolVar(f"X_{i}_{c}") for c in cells} for i in range(n)}


# ---------- 3D Core Constraints ----------

def add_exactly_one_cell_per_residue_3d(model: cp_model.CpModel, 
                                         X: Dict[Residue, Dict[Cell3D, cp_model.IntVar]], 
                                         cells: Sequence[Cell3D]) -> None:
    """
    For each residue i: sum_c X[i][c] == 1
    
    Args:
        model: CP-SAT model
        X: Placement variables
        cells: List of all cells
    """
    logger.info(f"Adding exactly one cell per residue for {len(X)} residues (3D)")
    console.print(f"[dim]Adding exactly one cell per residue for {len(X)} residues (3D)...[/dim]")
    
    for i in range(len(X)):
        model.Add(sum(X[i][c] for c in cells) == 1)


def add_at_most_one_residue_per_cell_3d(model: cp_model.CpModel,
                                          X: Dict[Residue, Dict[Cell3D, cp_model.IntVar]],
                                          cells: Sequence[Cell3D]) -> None:
    """
    For each cell c: sum_i X[i][c] <= 1
    
    Args:
        model: CP-SAT model
        X: Placement variables
        cells: List of all cells
    """
    logger.info(f"Adding at most one residue per cell for {len(X)} residues and {len(cells)} cells (3D)")
    console.print(f"[dim]Adding at most one residue per cell (3D)...[/dim]")
    
    for c in cells:
        model.Add(sum(X[i][c] for i in range(len(X))) <= 1)


def add_chain_adjacency_allowed_pairs_3d(model: cp_model.CpModel,
                                          X: Dict[Residue, Dict[Cell3D, cp_model.IntVar]],
                                          neighbors: Dict[Cell3D, List[Cell3D]],
                                          cells: Sequence[Cell3D]) -> None:
    """
    Enforce for each consecutive pair (i, i+1) that chosen cells are neighbors.
    Pattern: For each cell u: (X[i][u] == 1) => sum_{v in N(u)} X[i+1][v] == 1
    
    Args:
        model: CP-SAT model
        X: Placement variables
        neighbors: Neighbor map
        cells: List of all cells
    """
    logger.info(f"Adding chain adjacency constraints for {len(X)} residues (3D)")
    console.print(f"[dim]Adding chain adjacency constraints for {len(X)} residues (3D)...[/dim]")
    
    for i in range(len(X) - 1):
        for c in cells:
            model.Add(sum(X[i+1][v] for v in neighbors[c]) == 1).OnlyEnforceIf(X[i][c])


# ---------- 3D Symmetry Breaking ----------

def fix_start_and_orientation_3d(model: cp_model.CpModel,
                                  X: Dict[Residue, Dict[Cell3D, cp_model.IntVar]],
                                  start_cell: Cell3D,
                                  up_cell: Cell3D) -> None:
    """
    Fix residue 0 at start_cell and residue 1 at up_cell.
    Zero out all other one-hot variables for residues 0 and 1 for tight propagation.
    
    Args:
        model: CP-SAT model
        X: Placement variables
        start_cell: Cell for residue 0
        up_cell: Cell for residue 1
    """
    logger.info(f"Fixing residue 0 at {start_cell} and residue 1 at {up_cell} (3D)")
    console.print(f"[dim]Fixing residue 0 at {start_cell} and residue 1 at {up_cell} (3D)...[/dim]")
    
    # Fix residue 0 at start_cell, zero all others
    for c in X[0]:
        model.Add(X[0][c] == int(c == start_cell))
    
    # Fix residue 1 at up_cell, zero all others
    for c in X[1]:
        model.Add(X[1][c] == int(c == up_cell))


# ---------- 3D Contacts ----------

def declare_hh_contact_vars_3d(model: cp_model.CpModel,
                                seq: str,
                                X: Dict[Residue, Dict[Cell3D, cp_model.IntVar]],
                                neighbors: Dict[Cell3D, List[Cell3D]],
                                cells: Sequence[Cell3D]) -> Dict[Tuple[Cell3D, Cell3D], cp_model.IntVar]:
    """
    For each non-consecutive H-H pair (i, j) with |i-j|>1,
    create Bool C[i,j] meaning 'i and j occupy adjacent cells'.
    Links contact vars to placements using reified-AND equivalence.
    
    Args:
        model: CP-SAT model
        seq: Protein sequence string
        X: Placement variables
        neighbors: Neighbor map
        cells: List of all cells
        
    Returns:
        Dictionary mapping (i, j) pairs to contact BoolVars
    """
    n = len(seq)
    # Build undirected edges
    undirected_edges: List[Tuple[Cell3D, Cell3D]] = []
    seen = set()
    for u, nbrs in neighbors.items():
        for v in nbrs:
            key = (u, v) if u < v else (v, u)
            if key not in seen:
                seen.add(key)
                undirected_edges.append(key)

    # Chain edge usage indicators B[i,(u,v)]
    B: Dict[Tuple[int, Cell3D, Cell3D], cp_model.IntVar] = {}
    for i in range(n - 1):
        for u in cells:
            for v in neighbors[u]:
                b = model.NewBoolVar(f"B3_{i}_{u}_{v}")
                B[(i, u, v)] = b
                model.Add(b <= X[i][u])
                model.Add(b <= X[i+1][v])
            model.Add(sum(B[(i, u, v)] for v in neighbors[u]) == X[i][u])

    # Occupancy and H flags
    Occ: Dict[Cell3D, cp_model.IntVar] = {}
    H_at: Dict[Cell3D, cp_model.IntVar] = {}
    for u in cells:
        occ = model.NewIntVar(0, 1, f"Occ3_{u}")
        model.Add(occ == sum(X[i][u] for i in range(n)))
        Occ[u] = occ
        h_occ = model.NewIntVar(0, 1, f"H3_at_{u}")
        model.Add(h_occ == sum(X[i][u] for i in range(n) if seq[i] == 'H'))
        H_at[u] = h_occ

    # Contacts per undirected edge
    C_edge: Dict[Tuple[Cell3D, Cell3D], cp_model.IntVar] = {}
    for (u, v) in undirected_edges:
        c = model.NewBoolVar(f"C3edge_{u}_{v}")
        model.Add(c <= Occ[u])
        model.Add(c <= Occ[v])
        model.Add(c <= H_at[u])
        model.Add(c <= H_at[v])
        for i in range(n - 1):
            model.Add(c <= 1 - B[(i, u, v)])
            model.Add(c <= 1 - B[(i, v, u)])
        key = (u, v) if u < v else (v, u)
        C_edge[key] = c

    logger.info(f"Edge-based contacts (3D): {len(C_edge)} edges")
    console.print(f"[dim]Edge-based contacts active (3D): {len(C_edge)} edges[/dim]")
    return C_edge


# ---------- 3D Visualization ----------

def render_3d_html(seq: str,
                   positions: Dict[Residue, Cell3D],
                   contacts_hh: Iterable[Tuple[Residue, Residue]],
                   out_html: str = "folding_3d.html",
                   show_contacts: bool = True,
                   energy_epsHH: Optional[int] = None) -> str:
    """
    Produce a standalone, interactive HTML with pan/zoom/rotate using Plotly.
    
    Features:
    - Markers as spheres: H = green (#39b54a), P = gray (#b0b0b0)
    - Chain edges as sticks connecting (i, i+1)
    - Optional thin, semi-transparent sticks for H-H contacts
    - Equal aspect ratio, sensible camera, hover text showing index and letter
    
    Args:
        seq: Protein sequence string
        positions: Dictionary mapping residue index to (x, y, z) coordinates
        contacts_hh: Iterable of (i, j) H-H contact pairs
        out_html: Output HTML file path
        show_contacts: Whether to render contact edges
        
    Returns:
        Path to the saved HTML file
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        console.print("[red]Error: plotly not installed. Run: uv add plotly[/red]")
        raise
    
    logger.info(f"Rendering 3D visualization to {out_html}")
    console.print(f"[dim]Rendering 3D visualization to {out_html}...[/dim]")
    
    # Extract coordinates
    n = len(seq)
    xs = [positions[i][0] for i in range(n)]
    ys = [positions[i][1] for i in range(n)]
    zs = [positions[i][2] for i in range(n)]
    
    # Colors: H=green, P=gray
    colors = ['#39b54a' if seq[i] == 'H' else '#b0b0b0' for i in range(n)]
    hover_text = [f"Residue {i}: {seq[i]}" for i in range(n)]
    
    # Create figure
    fig = go.Figure()
    
    # Add chain edges (backbone)
    for i in range(n - 1):
        fig.add_trace(go.Scatter3d(
            x=[xs[i], xs[i+1]],
            y=[ys[i], ys[i+1]],
            z=[zs[i], zs[i+1]],
            mode='lines',
            line=dict(color='#404040', width=6),
            showlegend=False,
            hoverinfo='skip'
        ))
    
    # Add H-H contact edges (if requested)
    if show_contacts:
        for i, j in contacts_hh:
            fig.add_trace(go.Scatter3d(
                x=[xs[i], xs[j]],
                y=[ys[i], ys[j]],
                z=[zs[i], zs[j]],
                mode='lines',
                line=dict(color='rgba(255, 200, 0, 0.3)', width=3, dash='dot'),
                showlegend=False,
                hoverinfo='skip'
            ))
    
    # Add residue markers (spheres)
    fig.add_trace(go.Scatter3d(
        x=xs,
        y=ys,
        z=zs,
        mode='markers',
        marker=dict(
            size=12,
            color=colors,
            line=dict(color='white', width=1),
            opacity=0.95
        ),
        text=hover_text,
        hoverinfo='text',
        showlegend=False
    ))
    
    # Legend entries (H/P colors)
    legend_items = [
        dict(name="H (Hydrophobic)", marker=dict(color="#39b54a"), mode="markers"),
        dict(name="P (Polar)", marker=dict(color="#b0b0b0"), mode="markers"),
    ]

    # Update layout (with optional energy in title)
    fig.update_layout(
        title=dict(
            text=(f"HP Protein Folding: {seq}"
                  + (f" — Energy E/ε_HH: {energy_epsHH}" if energy_epsHH is not None else "")),
            x=0.5,
            xanchor='center'
        ),
        scene=dict(
            xaxis=dict(title='X', showgrid=True, gridcolor='#e0e0e0'),
            yaxis=dict(title='Y', showgrid=True, gridcolor='#e0e0e0'),
            zaxis=dict(title='Z', showgrid=True, gridcolor='#e0e0e0'),
            aspectmode='data',
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=1.5)
            ),
            bgcolor='#fafafa'
        ),
        margin=dict(l=0, r=0, b=0, t=40),
        hovermode='closest',
        paper_bgcolor='white'
    )

    # Add a simple legend panel using annotations
    fig.add_annotation(
        x=0.01, y=0.98, xref='paper', yref='paper', showarrow=False,
        text="<b>Legend</b><br>H: <span style='color:#39b54a;'>green</span><br>P: <span style='color:#b0b0b0;'>gray</span>",
        align='left', bgcolor='rgba(255,255,255,0.7)', bordercolor='#cccccc', borderwidth=1
    )
    
    # Write standalone HTML
    fig.write_html(out_html, include_plotlyjs='cdn')
    
    console.print(f"[green]✓ 3D visualization saved to: {out_html}[/green]")
    return out_html


def print_3d_slices(seq: str, positions: Dict[Residue, Cell3D], L: int) -> None:
    """
    Print Z-slices of the 3D folding for console debugging.
    
    Args:
        seq: Protein sequence string
        positions: Dictionary mapping residue index to (x, y, z) coordinates
        L: Grid dimension
    """
    console.print("\n[bold]3D Grid (Z-slices):[/bold]")
    
    # Group positions by z-level
    by_z = {}
    for i, (x, y, z) in positions.items():
        if z not in by_z:
            by_z[z] = {}
        by_z[z][(x, y)] = seq[i]
    
    # Print each z-slice
    for z in sorted(by_z.keys()):
        console.print(f"\n[cyan]Z = {z}:[/cyan]")
        grid = [["." for _ in range(L)] for _ in range(L)]
        for (x, y), letter in by_z[z].items():
            grid[y][x] = letter
        for row in grid:
            console.print(" ".join(row))

