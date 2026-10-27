from rdkit import Chem
from rdkit.Chem import AllChem


class ConformerGenerator:
    def __init__(self, mol, num_confs=100, prune_rms=0.5):
        self.mol = mol
        self.num_confs = num_confs
        self.prune_rms = prune_rms

    def _prepare_3d(self):
        mol_with_hs = Chem.AddHs(self.mol)
        params = AllChem.ETKDGv3()
        params.numThreads = 0
        params.pruneRmsThresh = self.prune_rms
        
        AllChem.EmbedMultipleConfs(mol_with_hs, numConfs=self.num_confs, params=params)
        return mol_with_hs

    def run_optimization(self, strategy):
        mol_3d = self._prepare_3d()
        res = strategy.optimize(mol_3d)
        return mol_3d, res