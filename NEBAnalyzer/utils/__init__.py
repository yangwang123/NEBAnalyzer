import ase
from ase.constraints import FixAtoms
from ase.geometry import wrap_positions
from typing import List


def center_atoms(atoms: List(ase.Atoms), center_indices: list=slice(None), n_image: int=0) -> List(ase.Atoms):
    
    a = atoms[n_image][center_indices]
    mean_xyz = a.get_center_of_mass()
    
    cell = atoms[n_image].get_cell()
    center = cell[0]/2 + cell[1]/2 + cell[2]/2
    shift = center - mean_xyz
    
    constraint = FixAtoms(indices=[])
    shifted_atoms = []
    for a in atoms:
        a = a.copy()
        a.set_constraint(constraint)
        pos = a.get_positions()
        pos[:, :] += shift
        pos = wrap_positions(pos, cell=cell, center=(0.5, 0.5, 0.5))
        a.set_positions(pos)
        shifted_atoms.append(a)
        
    return shifted_atoms


