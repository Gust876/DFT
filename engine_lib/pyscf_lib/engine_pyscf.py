from engine_lib.engine_interface import EngineDFT
from pyscf.gto import Mole


class EnginePySCF(EngineDFT):
    def __init__(self, mol: Mole, xc: str):
        self.mol = mol
        self.xc = xc
    
    def dft_method(self):
        mf = self.mol.KS().density_fit()
        mf.xc = self.xc
        
        mf.max_memory = 6000
        mf.conv_tol = 1e-6

        return mf