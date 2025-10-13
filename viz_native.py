"""
Native 3D Visualization for HP Protein Folding

Provides interactive desktop 3D viewer using PyVista (preferred) or vedo (fallback).
No web dependencies - creates native application windows.
"""

from typing import Dict, List, Tuple, Optional
import logging
from pathlib import Path

from rich.console import Console

logger = logging.getLogger(__name__)
console = Console()

# Try importing PyVista first, then vedo as fallback
try:
    import pyvista as pv
    BACKEND = "pyvista"
    console.print("[dim]Using PyVista for 3D visualization[/dim]")
except ImportError:
    try:
        import vedo
        BACKEND = "vedo"
        console.print("[dim]Using vedo for 3D visualization (PyVista not found)[/dim]")
    except ImportError:
        BACKEND = None
        console.print("[yellow]Warning: Neither PyVista nor vedo installed. "
                      "3D native visualization unavailable.[/yellow]")


def render_native_3d(
    seq: str,
    positions: Dict[int, Tuple],
    contacts_hh: Optional[List[Tuple[int, int]]] = None,
    title: Optional[str] = None,
    show_contacts: bool = True,
    sphere_radius: float = 0.20,
    tube_radius: float = 0.07,
    background: str = "white",
    energy_epsHH: Optional[int] = None
) -> None:
    """
    Open an interactive native 3D viewer window with ball-and-stick protein model.
    
    Features:
    - Spheres at residue positions (H=green, P=gray)
    - Chain backbone as tubes between consecutive residues
    - Optional H-H contact tubes (semi-transparent)
    - Smooth shading, anti-aliasing, trackball camera
    
    Keyboard shortcuts:
    - Space: Reset camera
    - S: Save screenshot to out/snap_<seq>.png
    - Q: Quit
    
    Args:
        seq: Protein sequence string
        positions: Dictionary mapping residue index to (x, y, z) coordinates
        contacts_hh: Optional list of (i, j) H-H contact pairs
        title: Window title (default: "HP Folding: <seq>")
        show_contacts: Whether to show H-H contact edges
        sphere_radius: Radius of residue spheres
        tube_radius: Radius of chain backbone tubes
        background: Background color
    """
    if BACKEND is None:
        console.print("[red]Error: No 3D visualization backend available. "
                      "Install PyVista: uv add pyvista[/red]")
        return
    
    if BACKEND == "pyvista":
        _render_pyvista(seq, positions, contacts_hh, title, show_contacts,
                        sphere_radius, tube_radius, background, energy_epsHH)
    else:
        _render_vedo(seq, positions, contacts_hh, title, show_contacts,
                     sphere_radius, tube_radius, background, energy_epsHH)


def save_native_3d(
    seq: str,
    positions: Dict[int, Tuple],
    path_png: str,
    contacts_hh: Optional[List[Tuple[int, int]]] = None,
    show_contacts: bool = True,
    sphere_radius: float = 0.20,
    tube_radius: float = 0.07,
    background: str = "white",
    resolution: Tuple[int, int] = (1920, 1080),
    energy_epsHH: Optional[int] = None
) -> str:
    """
    Render offscreen and save a high-quality PNG image.
    
    Args:
        seq: Protein sequence string
        positions: Dictionary mapping residue index to coordinates
        path_png: Output PNG file path
        contacts_hh: Optional list of (i, j) H-H contact pairs
        show_contacts: Whether to show H-H contact edges
        sphere_radius: Radius of residue spheres
        tube_radius: Radius of chain backbone tubes
        background: Background color
        resolution: Image resolution (width, height)
        
    Returns:
        Path to saved PNG file
    """
    if BACKEND is None:
        console.print("[red]Error: No 3D visualization backend available.[/red]")
        return ""
    
    Path(path_png).parent.mkdir(parents=True, exist_ok=True)
    
    if BACKEND == "pyvista":
        return _save_pyvista(seq, positions, path_png, contacts_hh, show_contacts,
                             sphere_radius, tube_radius, background, resolution, energy_epsHH)
    else:
        return _save_vedo(seq, positions, path_png, contacts_hh, show_contacts,
                          sphere_radius, tube_radius, background, resolution, energy_epsHH)


# ========== PyVista Implementation ==========

def _render_pyvista(
    seq: str,
    positions: Dict[int, Tuple],
    contacts_hh: Optional[List[Tuple[int, int]]],
    title: Optional[str],
    show_contacts: bool,
    sphere_radius: float,
    tube_radius: float,
    background: str,
    energy_epsHH: Optional[int]
) -> None:
    """Render using PyVista."""
    import pyvista as pv
    import numpy as np
    
    # Ensure 3D coordinates
    coords = []
    for i in range(len(seq)):
        pos = positions[i]
        if len(pos) == 2:
            coords.append((pos[0], pos[1], 0.0))
        else:
            coords.append(pos)
    coords = np.array(coords, dtype=float)
    
    # Create plotter
    pl = pv.Plotter()
    pl.background_color = background
    
    # Add residue spheres
    for i, (x, y, z) in enumerate(coords):
        color = '#39b54a' if seq[i] == 'H' else '#b0b0b0'
        sphere = pv.Sphere(radius=sphere_radius, center=(x, y, z), theta_resolution=30, phi_resolution=30)
        pl.add_mesh(sphere, color=color, smooth_shading=True, specular=0.5, specular_power=20)
    
    # Add chain backbone tubes
    n = len(coords)
    for i in range(n - 1):
        p0, p1 = coords[i], coords[i + 1]
        line = pv.Line(p0, p1)
        tube = line.tube(radius=tube_radius, n_sides=20)
        pl.add_mesh(tube, color='#7a7a7a', smooth_shading=True)
    
    # Add H-H contact tubes
    if show_contacts and contacts_hh:
        for e in contacts_hh:
            try:
                i, j = e
                if isinstance(i, int) and isinstance(j, int):
                    p0, p1 = coords[i], coords[j]
                else:
                    u, v = e
                    u3 = (u[0], u[1], (u[2] if len(u) > 2 else 0.0))
                    v3 = (v[0], v[1], (v[2] if len(v) > 2 else 0.0))
                    p0, p1 = u3, v3
            except Exception:
                continue
            line = pv.Line(p0, p1)
            tube = line.tube(radius=tube_radius * 0.5, n_sides=12)
            pl.add_mesh(tube, color='#39b54a', opacity=0.35, smooth_shading=True)
    
    # Configure camera and lighting
    pl.enable_anti_aliasing('msaa')
    pl.add_light(pv.Light(position=(10, 10, 10), light_type='scene light'))
    pl.camera.zoom(1.2)
    pl.reset_camera()
    
    # Set title and legend overlay
    window_title = title or f"HP Folding: {seq[:30]}{'...' if len(seq) > 30 else ''}"
    pl.add_text(window_title, position='upper_edge', font_size=10, color='black')
    pl.add_text("H: green    P: gray", position='lower_left', font_size=10, color='black',
                shadow=True, viewport=True)
    if energy_epsHH is not None:
        pl.add_text(f"Energy E/ε_HH: {energy_epsHH}", position='lower_right', font_size=10, color='black',
                    shadow=True, viewport=True)
    
    # Show
    console.print(f"[cyan]Opening native 3D viewer... (Press 's' to save screenshot, 'q' to quit)[/cyan]")
    pl.show()


def _save_pyvista(
    seq: str,
    positions: Dict[int, Tuple],
    path_png: str,
    contacts_hh: Optional[List[Tuple[int, int]]],
    show_contacts: bool,
    sphere_radius: float,
    tube_radius: float,
    background: str,
    resolution: Tuple[int, int],
    energy_epsHH: Optional[int]
) -> str:
    """Save PNG using PyVista offscreen rendering."""
    import pyvista as pv
    import numpy as np
    
    # Ensure 3D coordinates
    coords = []
    for i in range(len(seq)):
        pos = positions[i]
        if len(pos) == 2:
            coords.append((pos[0], pos[1], 0.0))
        else:
            coords.append(pos)
    coords = np.array(coords, dtype=float)
    
    # Create offscreen plotter
    pl = pv.Plotter(off_screen=True, window_size=resolution)
    pl.background_color = background
    
    # Add residue spheres
    for i, (x, y, z) in enumerate(coords):
        color = '#39b54a' if seq[i] == 'H' else '#b0b0b0'
        sphere = pv.Sphere(radius=sphere_radius, center=(x, y, z), theta_resolution=30, phi_resolution=30)
        pl.add_mesh(sphere, color=color, smooth_shading=True, specular=0.5, specular_power=20)
    
    # Add chain backbone tubes
    n = len(coords)
    for i in range(n - 1):
        p0, p1 = coords[i], coords[i + 1]
        line = pv.Line(p0, p1)
        tube = line.tube(radius=tube_radius, n_sides=20)
        pl.add_mesh(tube, color='#7a7a7a', smooth_shading=True)
    
    # Add H-H contact tubes
    if show_contacts and contacts_hh:
        for i, j in contacts_hh:
            p0, p1 = coords[i], coords[j]
            line = pv.Line(p0, p1)
            tube = line.tube(radius=tube_radius * 0.5, n_sides=12)
            pl.add_mesh(tube, color='#39b54a', opacity=0.35, smooth_shading=True)
    
    # Configure camera and lighting
    pl.enable_anti_aliasing('msaa')
    pl.add_light(pv.Light(position=(10, 10, 10), light_type='scene light'))
    pl.camera.zoom(1.2)
    pl.reset_camera()
    
    # Overlays
    pl.add_text("H: green    P: gray", position='lower_left', font_size=10, color='black',
                shadow=True, viewport=True)
    if energy_epsHH is not None:
        pl.add_text(f"Energy E/ε_HH: {energy_epsHH}", position='lower_right', font_size=10, color='black',
                    shadow=True, viewport=True)

    # Save screenshot
    pl.screenshot(path_png, transparent_background=False)
    
    console.print(f"[green]✓ Screenshot saved to: {path_png}[/green]")
    return path_png


# ========== Vedo Implementation ==========

def _render_vedo(
    seq: str,
    positions: Dict[int, Tuple],
    contacts_hh: Optional[List[Tuple[int, int]]],
    title: Optional[str],
    show_contacts: bool,
    sphere_radius: float,
    tube_radius: float,
    background: str,
    energy_epsHH: Optional[int]
) -> None:
    """Render using vedo."""
    from vedo import Plotter, Sphere, Tube, Points, Text2D
    import numpy as np
    
    # Ensure 3D coordinates
    coords = []
    for i in range(len(seq)):
        pos = positions[i]
        if len(pos) == 2:
            coords.append([pos[0], pos[1], 0.0])
        else:
            coords.append(list(pos))
    coords = np.array(coords, dtype=float)
    
    # Create objects
    objects = []
    
    # Add residue spheres
    for i, (x, y, z) in enumerate(coords):
        color = '#39b54a' if seq[i] == 'H' else '#b0b0b0'
        sphere = Sphere(pos=(x, y, z), r=sphere_radius, c=color)
        sphere.phong()
        objects.append(sphere)
    
    # Add chain backbone tubes
    n = len(coords)
    for i in range(n - 1):
        p0, p1 = coords[i], coords[i + 1]
        tube = Tube([p0, p1], r=tube_radius, c='#7a7a7a')
        tube.phong()
        objects.append(tube)
    
    # Add H-H contact tubes
    if show_contacts and contacts_hh:
        for e in contacts_hh:
            try:
                i, j = e
                if isinstance(i, int) and isinstance(j, int):
                    p0, p1 = coords[i], coords[j]
                else:
                    u, v = e
                    u3 = (u[0], u[1], (u[2] if len(u) > 2 else 0.0))
                    v3 = (v[0], v[1], (v[2] if len(v) > 2 else 0.0))
                    p0, p1 = u3, v3
            except Exception:
                continue
            tube = Tube([p0, p1], r=tube_radius * 0.5, c='#39b54a', alpha=0.35)
            tube.phong()
            objects.append(tube)
    
    # Create plotter and show
    window_title = title or f"HP Folding: {seq[:30]}{'...' if len(seq) > 30 else ''}"
    plt = Plotter(bg=background, title=window_title, axes=0)
    legend = Text2D("H: green    P: gray", pos='bottom-left', s=0.8, c='black', bg='white', alpha=0.7)
    objects.append(legend)
    if energy_epsHH is not None:
        energy_text = Text2D(f"Energy E/ε_HH: {energy_epsHH}", pos='bottom-right', s=0.8, c='black', bg='white', alpha=0.7)
        objects.append(energy_text)
    
    console.print(f"[cyan]Opening native 3D viewer... (Press 'q' to quit)[/cyan]")
    plt.show(*objects, viewup="z")


def _save_vedo(
    seq: str,
    positions: Dict[int, Tuple],
    path_png: str,
    contacts_hh: Optional[List[Tuple[int, int]]],
    show_contacts: bool,
    sphere_radius: float,
    tube_radius: float,
    background: str,
    resolution: Tuple[int, int],
    energy_epsHH: Optional[int]
) -> str:
    """Save PNG using vedo offscreen rendering."""
    from vedo import Plotter, Sphere, Tube
    import numpy as np
    
    # Ensure 3D coordinates
    coords = []
    for i in range(len(seq)):
        pos = positions[i]
        if len(pos) == 2:
            coords.append([pos[0], pos[1], 0.0])
        else:
            coords.append(list(pos))
    coords = np.array(coords, dtype=float)
    
    # Create objects
    objects = []
    
    # Add residue spheres
    for i, (x, y, z) in enumerate(coords):
        color = '#39b54a' if seq[i] == 'H' else '#b0b0b0'
        sphere = Sphere(pos=(x, y, z), r=sphere_radius, c=color)
        sphere.phong()
        objects.append(sphere)
    
    # Add chain backbone tubes
    n = len(coords)
    for i in range(n - 1):
        p0, p1 = coords[i], coords[i + 1]
        tube = Tube([p0, p1], r=tube_radius, c='#7a7a7a')
        tube.phong()
        objects.append(tube)
    
    # Add H-H contact tubes
    if show_contacts and contacts_hh:
        for e in contacts_hh:
            try:
                i, j = e
                if isinstance(i, int) and isinstance(j, int):
                    p0, p1 = coords[i], coords[j]
                else:
                    u, v = e
                    u3 = (u[0], u[1], (u[2] if len(u) > 2 else 0.0))
                    v3 = (v[0], v[1], (v[2] if len(v) > 2 else 0.0))
                    p0, p1 = u3, v3
            except Exception:
                continue
            tube = Tube([p0, p1], r=tube_radius * 0.5, c='#39b54a', alpha=0.35)
            tube.phong()
            objects.append(tube)
    
    # Create offscreen plotter
    plt = Plotter(bg=background, offscreen=True, size=resolution, axes=0)
    # Overlays
    from vedo import Text2D as _Text2D
    objects.append(_Text2D("H: green    P: gray", pos='bottom-left', s=0.8, c='black', bg='white', alpha=0.7))
    if energy_epsHH is not None:
        objects.append(_Text2D(f"Energy E/ε_HH: {energy_epsHH}", pos='bottom-right', s=0.8, c='black', bg='white', alpha=0.7))
    plt.show(*objects, viewup="z")
    plt.screenshot(path_png)
    plt.close()
    
    console.print(f"[green]✓ Screenshot saved to: {path_png}[/green]")
    return path_png

