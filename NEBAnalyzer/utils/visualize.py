import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import nglview as nv


def show_atoms(atoms) -> nv.NGLWidget:
    
    try:
        widget = nv.show_ase(atoms)
    except Exception:
        widget = nv.show_asetraj(atoms)
    widget.clear_representations()
    widget.add_ball_and_stick()
    widget.add_label(labelType='atomindex', color='blue')
    widget.add_unitcell()
    widget.camera = 'orthographic'
    
    return widget


def show_E(E_list: list, coords: list) -> None:
    
    coords_spline = np.arange(0, coords[-1], 0.001)
    E_spline = CubicSpline(coords, E_list, bc_type='natural')(coords_spline)
    E_mark = [(np.abs(coords_spline - val)).argmin() for val in coords]

    fig, ax = plt.subplots(1, 1, figsize=[6, 4])
    ax.plot(coords_spline, E_spline, 'o--', color='darkred', linewidth=2, markersize=12, markevery=E_mark)
    ax.set_ylabel('Energy (kcal/mol)')
    ax.set_xlabel('Reaction coordinate (Å)')
    fig.tight_layout()
    fig.set_dpi(100)


