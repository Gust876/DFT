from abc import ABC, abstractmethod
from rdkit.Chem import Mol
 
 
class ForceField(ABC):
    
    @abstractmethod
    def optimize(self, mol: Mol):
        pass