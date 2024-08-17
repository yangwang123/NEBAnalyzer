import numpy as np
import nglview as nv
from typing import Tuple, List

from .utils.visualize import show_atoms, show_E
from .utils.units import eV_to_kcal_per_mol 


class Analyzer:
    '''
    Template analyer for various QM codes

    Parameters:
    -----------
    ddir: str
        Directory of the NEB calculation
    read_neb: bool
        Whether to read the NEB calculation
    '''
    def __init__(self, ddir: str, read_neb: bool=True):
        self.ddir = ddir
        self.get_n_images()
        if read_neb:
            self.get_E_ini()
            self.get_E_fin()
            self.get_E_all()
            self.get_dists_all()
            self.get_forces()
     
    def get_n_images(self) -> None:
        '''Get the number of images in the NEB calculation; override this method'''
        return
        
    def get_E_ini(self) -> None:
        '''Get energy of the first image; override this method'''
        return
        
    def get_E_fin(self) -> None:
        '''Get energy of the last image; override this method'''
        return
        
    def get_E_all(self) -> None:
        '''Get energy of all images; override this method'''
        return
    
    def get_dists_all(self) -> None:
        '''Get distances of all images; override this method'''
        return
    
    def get_forces(self) -> None:
        '''Get forces of all images; override this method'''
        return
    
    def get_pathway(self, ndx: int=-1, initial: bool=False) -> List:
        '''Get pathway for visualization; override this method'''
        return []
    
    def get_min_force(self) -> Tuple[int, float]:
        '''Get minimum force and its index

        Returns:
        --------
        ndx: int
            Index of the minimum force
        f: float
            Magnitude of the minimum force
        '''
        ndx = np.argmin(np.max(self.forces, axis=-1))
        f = self.forces[ndx, :]
        
        return ndx, f
    
    def get_E(self, ndx: int=-1, ref: str='min') -> Tuple[list, list]:
        '''Get energy and coordinates for plotting
        
        Parameters:
        -----------
        ndx: int
            Pathway index
        ref: str
            Reference energy
        
        Returns:
        --------
        E_list: list
            List of energies
        coords: list
            List of coordinates
        '''
        E_list = [self.E_ini] + list(self.E_all[ndx, :]) + [self.E_fin]
        if ref == 'min':
            E_ref = np.min(E_list)
        elif ref == 'fin':
            E_ref = E_list[-1]
        else:
            E_ref = E_list[0]
        E_list = np.array(E_list) - E_ref
        E_list = list(E_list*eV_to_kcal_per_mol)
        
        dists = list(self.dists_all[ndx, :])
        coords = [sum(dists[:i]) for i in range(len(dists)+1)]
        
        return E_list, coords
    
    
    def plot_E(self, ndx: int=-1, ref: str='min') -> None:
        '''Plot energy profile

        Parameters:
        -----------
        ndx: int
            Pathway index
        ref: str
            Reference energy to set zero
        '''
        E_list, coords = self.get_E(ndx=ndx, ref=ref)
        return show_E(E_list, coords)
        
        
    def vis_pathway(self, ndx: int=-1, initial: bool=False) -> nv.NGLWidget:
        '''Visualize pathway
        
        Parameters:
        -----------
        ndx: int
            Pathway index
        initial: bool
            Whether to show the initial state
        '''
        atoms = self.get_pathway(ndx, initial)
        return show_atoms(atoms)

