import os
import numpy as np
import ase
from ase import io
from typing import List

from . import Analyzer


class VaspAnalyzer(Analyzer):
    '''
    Analyzer for VASP NEB calculations

    Parameters:
    -----------
    ddir: str
        Directory of the NEB calculation
    read_neb: bool
        Whether to read the NEB calculation
    '''
    def __init__(self, ddir: str, read_neb: bool=True):
        super().__init__(ddir, read_neb)
        
    def get_n_images(self) -> None:
        '''Get the number of images in the NEB calculation'''
        nums = set(str(i) for i in range(10))
        self.n_images = len([d for d in os.listdir(self.ddir) if len(d) == 2
                             and os.path.isdir(os.path.join(self.ddir, d))
                             and d[0] in nums and d[1] in nums])
    
    @staticmethod
    def get_E_image(file: str) -> List[float]:
        '''Get the total energy of each image
        
        Parameters:
        -----------
        file: str
            OUTCAR file of the image
        
        Returns:
        --------
        Es: list
            Total energies of images
        '''
        try:
            atoms = io.read(file, format='vasp-out', index=slice(None))
            return [a.get_total_energy() for a in atoms]
        except Exception:
            Es = []
            with open(file, 'r') as f:
                lines = f.readlines()
            for line in lines:
                if 'energy  without entropy=' in line:
                    line = line.split()
                    Es.append(float(line[-1]))
            return Es
        
    def get_E_ini(self) -> None:
        '''Get energy of the first image'''
        file = os.path.join(self.ddir, f'00/OUTCAR')
        self.E_ini = self.get_E_image(file)[-1]
        
    def get_E_fin(self) -> None:
        '''Get energy of the last image'''
        file = os.path.join(self.ddir, f'{self.n_images-1:02d}/OUTCAR')
        self.E_fin = self.get_E_image(file)[-1]
        
    def get_E_all(self) -> None:
        '''Get energy of intermediate images'''
        E_all = []
        for i in range(self.n_images):
            if i == 0 or i == self.n_images-1:
                continue
            file = os.path.join(self.ddir, f'{i:02d}/OUTCAR')
            E_all.append(self.get_E_image(file))
        self.E_all = np.array(E_all).T
    
    @staticmethod
    def get_dists_image(file: str, prev: bool=False) -> List[float]:
        '''Get the distance between the current image and the previous or next image
        
        Parameters:
        -----------
        file: str
            OUTCAR file of the image
        prev: bool
            Whether to get the distance to the previous image
        
        Returns:
        --------
        dists: list
            Distance to the previous or next image
        '''
        dists = []
        with open(file, 'r') as f:
            for line in f:
                if 'NEB: distance to prev, next image, angle between' in line:
                    line = line.split()
                    if prev:
                        d = float(line[8])
                    else:
                        d = float(line[9])
                    dists.append(d)
        return dists
    
    def get_dists_all(self) -> None:
        '''Get distances between images'''
        dists_all = []
        for i in range(self.n_images):
            if i == 0 or i == self.n_images-1:
                continue
            file = os.path.join(self.ddir, f'{i:02d}/OUTCAR')
            if i == 1:
                dists_all.append(self.get_dists_image(file, prev=True))
            dists_all.append(self.get_dists_image(file, prev=False))
        
        self.dists_all = np.array(dists_all).T
    
    @staticmethod
    def get_forces_image(file: str) -> List[float]:
        '''Get the maximum force of an image
        
        Parameters:
        -----------
        file: str
            OUTCAR file of the image
        
        Returns:
        --------
        forces: list
            Maximum force at each iteration
        '''
        with open(file, 'r') as f:
            lines = f.readlines()
        forces = []
        for line in lines:
            if 'FORCES: max atom, RMS' in line:
                forces.append(float(line.strip().split()[4]))
        return forces
    
    def get_forces(self) -> None:
        '''Get maximum forces of intermediate images'''
        forces_list = []
        for i in range(self.n_images):
            if i == 0 or i == self.n_images-1:
                continue
            file = os.path.join(self.ddir, f'{i:02d}/OUTCAR')
            forces_list.append(self.get_forces_image(file))
        self.forces = np.array(forces_list).T
    
    def get_pathway(self, ndx: int=-1) -> List[ase.Atoms]:
        '''Get the pathway of the NEB calculation
        
        Parameters:
        -----------
        ndx: int
            Pathway index
            
        Returns:
        --------
        atoms: list
            List of ASE atoms objects
        '''
        atoms = []
        for i in range(self.n_images):
            if i == 0 or i == self.n_images-1 or ndx == 0:
                file = os.path.join(self.ddir, f'{i:02d}/POSCAR')
                a = io.read(file, format='vasp')
            elif ndx == -1:
                file = os.path.join(self.ddir, f'{i:02d}/CONTCAR')
                a = io.read(file, format='vasp')
            else:
                file = os.path.join(self.ddir, f'{i:02d}/OUTCAR')
                a = io.read(file, format='vasp-out', index=ndx)
            atoms.append(a)
        return atoms

