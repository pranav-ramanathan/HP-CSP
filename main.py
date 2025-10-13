"""
HP Protein Folding Solver - Unified CLI

Supports both 2D and 3D lattice folding with optional interactive visualization.
"""

from typing import Literal, Optional, List
import logging
from pathlib import Path
from ortools.sat.python import cp_model
from rich.console import Console
import typer

# Import modules
import hp2d
import hp3d
import hp_eval
import viz_native
import hp_parse

logger = logging.getLogger(__name__)
console = Console()
app = typer.Typer(rich_markup_mode="rich")


@app.command(name="run-hp")
def run_hp(
    seq: str = typer.Argument(..., help="Sequence string (plain or compressed, e.g., 'HPH2P10' or 'HPH₂P₁₀')"),
    L: int = typer.Option(4, "--grid-size", "-L", help="Grid size (L×L for 2D, L×L×L for 3D)"),
    time_limit: float = typer.Option(5.0, "--time", "-t", help="Solver time limit in seconds"),
    workers: int = typer.Option(8, "--workers", "-w", help="Number of parallel workers"),
    dim: Literal[2, 3] = typer.Option(2, "--dim", "-d", help="Dimension: 2 for 2D lattice, 3 for 3D lattice"),
    viz: Literal["none", "3d", "native"] = typer.Option("none", "--viz", "-v", help="Visualization: 'none', '3d' (Plotly HTML), or 'native' (desktop window)"),
    snap: Optional[str] = typer.Option(None, "--snap", help="Save snapshot to PNG file (offscreen rendering)"),
    strict: bool = typer.Option(False, "--strict", help="Error on any char outside [H P 0-9 Unicode subscripts whitespace]"),
    print_expanded: bool = typer.Option(False, "--print-expanded", help="Print the expanded sequence (clipped to 200 chars)"),
    expect_length: Optional[int] = typer.Option(None, "--expect-length", help="Assert expanded length equals this value")
):
    """
    [bold cyan]HP protein folding solver using OR-Tools CP-SAT.[/bold cyan]
    
    Finds the optimal folding of a protein sequence on a 2D or 3D lattice
    by maximizing H-H contacts. Reports energy in ε_HH units.
    
    [yellow]Examples:[/yellow]
      • 2D folding (default):
        [dim]$ python main.py run-hp HPPH -L 4[/dim]
        
      • 3D folding:
        [dim]$ python main.py run-hp HPPH -L 4 --dim 3[/dim]
        
      • 3D with Plotly HTML:
        [dim]$ python main.py run-hp HPPH -L 6 --dim 3 --viz 3d[/dim]
        
      • 3D with native viewer:
        [dim]$ python main.py run-hp HPPH -L 6 --dim 3 --viz native[/dim]
        
      • Save snapshot:
        [dim]$ python main.py run-hp HPPH -L 6 --dim 3 --snap out/fold.png[/dim]
    """
    assert L >= 3, "Use L>=3 so symmetry 'up' is valid"

    # Parse compressed/plain sequence
    try:
        expanded, warnings = hp_parse.expand_hp_sequence(seq, strict=strict)
    except hp_parse.ParseError as e:
        console.print(f"[red]Parse error:[/red] {e}")
        raise typer.Exit(1)

    if expect_length is not None and len(expanded) != expect_length:
        console.print(f"[yellow]Warning: Expanded length {len(expanded)} != expected {expect_length}[/yellow]")

    # Pretty print summary
    frac = hp_parse.hydrophobic_fraction(expanded)
    console.print(f"[dim]Parsed sequence (compressed): {seq[:40]}{'...' if len(seq) > 40 else ''}[/dim]")
    console.print(f"[dim]Expanded length: {len(expanded)} | Hydrophobic fraction: {frac:.2f}[/dim]")
    if print_expanded:
        clip = expanded if len(expanded) <= 200 else expanded[:200] + '...'
        console.print(f"[dim]Expanded (clip): {clip}[/dim]")

    # Replace input with expanded for solver
    seq = expanded
    
    if dim == 2:
        # ========== 2D MODE (UNCHANGED) ==========
        run_hp_2d(seq, L, time_limit, workers, viz, snap)
    elif dim == 3:
        # ========== 3D MODE ==========
        run_hp_3d(seq, L, time_limit, workers, viz, snap)
    else:
        console.print(f"[red]Error: Invalid dimension {dim}. Use 2 or 3.[/red]")
        raise typer.Exit(1)


def run_hp_2d(seq: str, L: int, time_limit: float, workers: int, viz: str, snap: Optional[str]):
    """
    Run 2D HP folding with optional visualization.
    """
    model = cp_model.CpModel()

    cells = hp2d.build_grid(L)
    nbrs = hp2d.neighbors_map(cells)
    start, up = hp2d.center_and_up(L)

    # Vars
    X = hp2d.declare_vars(model, n=len(seq), cells=cells)

    # Core constraints
    hp2d.add_exactly_one_cell_per_residue(model, X, cells)
    hp2d.add_at_most_one_residue_per_cell(model, X, cells)
    hp2d.add_chain_adjacency_allowed_pairs(model, X, nbrs, cells)

    # Symmetry
    hp2d.fix_start_and_orientation(model, X, start, up)

    # Contacts + objective (edge-based)
    C_edges = hp2d.declare_hh_contact_vars(model, seq, X, nbrs, cells)
    hp2d.set_objective_max_contacts(model, C_edges.values())

    # Solve
    solver = cp_model.CpSolver()
    hp2d.configure_solver(solver, time_limit_s=time_limit, workers=workers)
    status = solver.Solve(model)
    console.print(f"[dim]Solve time: {solver.WallTime():.3f}s[/dim]")

    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        console.print("[red]No solution found.[/red]")
        return

    pos = hp2d.extract_positions(solver, X)
    # Derive realized contacts from edge variables
    contacts_realized = [e for e, var in C_edges.items() if solver.Value(var) == 1]
    num_contacts = len(contacts_realized)
    energy = hp_eval.compute_energy_contacts(contacts_realized)

    console.print(f"[green]Status: {solver.StatusName(status)}[/green]")
    console.print(f"[bold cyan]Energy E/ε_HH: {energy}[/bold cyan]")
    console.print(f"[cyan]Contacts: {num_contacts}[/cyan]")
    console.print(f"[dim]Positions: {pos}[/dim]")

    # Pretty print grid
    grid = [["." for _ in range(L)] for __ in range(L)]
    for i, (x, y) in pos.items():
        grid[y][x] = seq[i]
    
    console.print("\n[bold]Grid:[/bold]")
    for row in grid:
        console.print(" ".join(row))
    
    # Auto-save PNG for 2D viz unless overridden by --snap
    auto_png = f"folding_2d_{seq[:10]}.png"
    
    # Handle visualization for 2D (convert to 3D with z=0)
    if viz == "native" or snap:
        pos_3d = {i: (x, y, 0) for i, (x, y) in pos.items()}
        if snap:
            viz_native.save_native_3d(seq, pos_3d, snap, contacts_realized, show_contacts=True, energy_epsHH=energy)
        else:
            # Save default PNG for 2D runs
            viz_native.save_native_3d(seq, pos_3d, auto_png, contacts_realized, show_contacts=True, energy_epsHH=energy)
        if viz == "native":
            viz_native.render_native_3d(seq, pos_3d, contacts_realized, show_contacts=True, energy_epsHH=energy)


def run_hp_3d(seq: str, L: int, time_limit: float, workers: int, viz: str, snap: Optional[str]):
    """
    Run 3D HP folding with optional visualization.
    """
    model = cp_model.CpModel()

    cells = hp3d.build_grid_3d(L)
    nbrs = hp3d.neighbors_map_3d(cells)
    start, up = hp3d.center_and_up_3d(L)

    # Vars
    X = hp3d.declare_vars_3d(model, n=len(seq), cells=cells)

    # Core constraints
    hp3d.add_exactly_one_cell_per_residue_3d(model, X, cells)
    hp3d.add_at_most_one_residue_per_cell_3d(model, X, cells)
    hp3d.add_chain_adjacency_allowed_pairs_3d(model, X, nbrs, cells)

    # Symmetry
    hp3d.fix_start_and_orientation_3d(model, X, start, up)

    # Contacts + objective (edge-based)
    C_edges = hp3d.declare_hh_contact_vars_3d(model, seq, X, nbrs, cells)
    hp2d.set_objective_max_contacts(model, C_edges.values())  # Same as 2D

    # Solve
    solver = cp_model.CpSolver()
    hp2d.configure_solver(solver, time_limit_s=time_limit, workers=workers)  # Same as 2D
    status = solver.Solve(model)
    console.print(f"[dim]Solve time: {solver.WallTime():.3f}s[/dim]")

    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        console.print("[red]No solution found.[/red]")
        return

    # Extract positions (3D)
    pos = {i: next(c for c in X[i].keys() if solver.Value(X[i][c]) == 1) for i in range(len(X))}
    
    # Extract realized H-H contacts from edges
    realized_contacts = [e for e, cij in C_edges.items() if solver.Value(cij) == 1]
    contacts_count = len(realized_contacts)
    energy = hp_eval.compute_energy_contacts([]) - 0  # compute from count
    energy = -contacts_count

    console.print(f"[green]Status: {solver.StatusName(status)}[/green]")
    console.print(f"[bold cyan]Energy E/ε_HH: {energy}[/bold cyan]")
    console.print(f"[cyan]Contacts: {contacts_count}[/cyan]")
    console.print(f"[dim]Positions: {pos}[/dim]")

    # Print Z-slices for debugging
    hp3d.print_3d_slices(seq, pos, L)

    # Visualization
    if viz == "3d":
        hp3d.render_3d_html(
            seq=seq,
            positions=pos,
            contacts_hh=[],
            out_html=f"folding_3d_{seq[:10]}.html",
            show_contacts=True,
            energy_epsHH=energy
        )
    elif viz == "native":
        viz_native.render_native_3d(
            seq=seq,
            positions=pos,
            contacts_hh=[],
            show_contacts=True,
            energy_epsHH=energy
        )
    
    # Handle snapshot
    if snap:
        viz_native.save_native_3d(
            seq=seq,
            positions=pos,
            path_png=snap,
            contacts_hh=[],
            show_contacts=True,
            energy_epsHH=energy
        )


@app.command()
def batch(
    sequences: List[str] = typer.Option([], "--seq", "-s", help="Sequence(s) to evaluate (can be repeated)"),
    file: Optional[str] = typer.Option(None, "--file", "-f", help="File with sequences (one per line)"),
    L: int = typer.Option(6, "--grid-size", "-L", help="Grid size"),
    dim: Literal[2, 3] = typer.Option(3, "--dim", "-d", help="Dimension: 2 or 3"),
    time_limit: float = typer.Option(120.0, "--time", "-t", help="Solver time limit per run (seconds)"),
    workers: int = typer.Option(8, "--workers", "-w", help="Number of parallel workers"),
    runs: int = typer.Option(1, "--runs", "-r", help="Number of runs per sequence"),
    csv_path: Optional[str] = typer.Option(None, "--csv", help="Output CSV file path"),
    json_path: Optional[str] = typer.Option(None, "--json", help="Output JSON file path"),
    latex_path: Optional[str] = typer.Option(None, "--latex", help="Output LaTeX file path"),
    seeds: List[int] = typer.Option([], "--seed", help="Random seeds (one per run, can be repeated)")
):
    """
    [bold cyan]Batch evaluation of multiple HP sequences.[/bold cyan]
    
    Evaluates multiple sequences and reports energy in ε_HH units.
    Optionally exports results to CSV, JSON, and LaTeX formats.
    
    [yellow]Examples:[/yellow]
      • From file:
        [dim]$ python main.py batch --file sequences.txt -L 6 --dim 3 --time 120[/dim]
        
      • Multiple sequences:
        [dim]$ python main.py batch -s HPPH -s HHPHPH -L 6 --dim 3[/dim]
        
      • With exports:
        [dim]$ python main.py batch --file seqs.txt -L 6 --dim 3 --runs 3 \\
          --csv out/results.csv --json out/results.json --latex out/results.tex[/dim]
    """
    # Load sequences
    all_sequences = list(sequences)
    
    if file:
        file_path = Path(file)
        if not file_path.exists():
            console.print(f"[red]Error: File not found: {file}[/red]")
            raise typer.Exit(1)
        
        with open(file_path) as f:
            file_sequences = [line.strip() for line in f if line.strip() and not line.startswith('#')]
            all_sequences.extend(file_sequences)
    
    if not all_sequences:
        console.print("[red]Error: No sequences provided. Use --seq or --file.[/red]")
        raise typer.Exit(1)
    
    # Validate sequences
    for seq in all_sequences:
        if not all(c in 'HP' for c in seq):
            console.print(f"[red]Error: Invalid sequence '{seq}'. Must contain only H and P.[/red]")
            raise typer.Exit(1)
    
    # Run batch evaluation
    results = hp_eval.evaluate_batch(
        sequences=all_sequences,
        L=L,
        dim=dim,
        time_limit_s=time_limit,
        workers=workers,
        seeds=seeds if seeds else None,
        runs_per_seq=runs
    )
    
    # Print summary table
    hp_eval.print_summary_table(results)
    
    # Export results
    if csv_path:
        hp_eval.save_results_csv(results, csv_path)
    
    if json_path:
        hp_eval.save_results_json(results, json_path)
    
    if latex_path:
        latex_code = hp_eval.latex_table(results)
        Path(latex_path).parent.mkdir(parents=True, exist_ok=True)
        with open(latex_path, 'w') as f:
            f.write(latex_code)
        console.print(f"[green]✓ LaTeX table saved to: {latex_path}[/green]")
    
    # Summary statistics
    total_sequences = len(results)
    optimal_count = sum(1 for r in results if r['is_optimal'])
    avg_energy = sum(r['best_energy_epsHH'] for r in results) / total_sequences if total_sequences > 0 else 0
    
    console.print(f"\n[bold]Summary:[/bold]")
    console.print(f"  Total sequences: {total_sequences}")
    console.print(f"  Optimal solutions: {optimal_count} ({optimal_count/total_sequences*100:.1f}%)")
    console.print(f"  Average energy: {avg_energy:.2f} ε_HH")


if __name__ == "__main__":
    app()

