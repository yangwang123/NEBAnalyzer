import os
import numpy as np
from ase import io

from .analysis import Analyzer


class VaspAnalyzer(Analyzer):
    
    
    def __init__(self, ddir: str):

        super().__init__(ddir)
        
    
    def get_n_images(self) -> None:
        
        nums = set(str(i) for i in range(10))
        self.n_images = len([d for d in os.listdir(self.ddir) if len(d) == 2
                             and os.path.isdir(os.path.join(self.ddir, d))
                             and d[0] in nums and d[1] in nums])
    
    
    @staticmethod
    def get_E_image(file: str) -> list:
        
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
        
        file = os.path.join(self.ddir, f'00/OUTCAR')
        self.E_ini = self.get_E_image(file)[-1]
        
        
    def get_E_fin(self) -> None:
        
        file = os.path.join(self.ddir, f'{self.n_images-1:02d}/OUTCAR')
        self.E_fin = self.get_E_image(file)[-1]
        
    
    def get_E_all(self) -> None:
        
        E_all = []
        for i in range(self.n_images):
            if i == 0 or i == self.n_images-1:
                continue
            file = os.path.join(self.ddir, f'{i:02d}/OUTCAR')
            E_all.append(self.get_E_image(file))
        
        self.E_all = np.array(E_all).T
    
        
    @staticmethod
    def get_dists_image(file: str, prev: bool=False) -> list:
        
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
    def get_forces_image(file) -> list:
        
        with open(file, 'r') as f:
            lines = f.readlines()
        forces = []
        for line in lines:
            if 'FORCES: max atom, RMS' in line:
                forces.append(float(line.strip().split()[4]))
        
        return forces
    
    
    def get_forces(self) -> None:
        
        forces_list = []
        for i in range(self.n_images):
            if i == 0 or i == self.n_images-1:
                continue
            file = os.path.join(self.ddir, f'{i:02d}/OUTCAR')
            forces_list.append(self.get_forces_image(file))
        
        self.forces = np.array(forces_list).T
    
    
    def get_pathway(self, ndx: int=-1, initial: bool=False) -> list:
        
        atoms = []
        for i in range(self.n_images):
            if i == 0 or i == self.n_images-1 or initial:
                file = os.path.join(self.ddir, f'{i:02d}/POSCAR')
                a = io.read(file, format='vasp')
            else:
                file = os.path.join(self.ddir, f'{i:02d}/OUTCAR')
                a = io.read(file, format='vasp-out', index=ndx)
            atoms.append(a)
            
        return atoms


