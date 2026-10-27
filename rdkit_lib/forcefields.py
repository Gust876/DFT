from rdkit.Chem import Mol, AllChem
from interfaces.forcefield_interface import ForceField
 
 
class MMFF94Strategy(ForceField):
    name = "MMFF94"
 
    def optimize(self, mol: Mol):
        return AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, maxIters=500)
 
 
class UFFStrategy(ForceField):
    name = "UFF"
    
    def optimize(self, mol: Mol):
        return AllChem.UFFOptimizeMoleculeConfs(mol, numThreads=0)