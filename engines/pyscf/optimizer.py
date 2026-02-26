from pyscf.geomopt.geometric_solver import optimize
from interfaces.engine_interface import EngineDFT
from interfaces.optimizer_interface import FactoryOptimizer
from interfaces.optimizer_interface import OptimizerProduct
from engines.pyscf.engine_pyscf import EnginePySCF
from pyscf.gto import Mole


class FactoryOptimizerPySCF(FactoryOptimizer):

    def __init__(self, mol: Mole, xc: str):
        self.mol = mol
        self.xc = xc

    def _initialize_engine(self):
        return EnginePySCF(mol=self.mol, xc=self.xc)

    def factory_method(self):
        return OptimizerPySCF(
            engine=self._initialize_engine()
            )
    

class OptimizerPySCF(OptimizerProduct):
    def __init__(self, engine: EngineDFT):
        self._engine = engine

    def opt_geometry(self, maxsteps: int, verbose: int):
        try:
            mf = self._engine.dft_method()

            opt_mol = optimize(
                mf,
                maxsteps = maxsteps,
                verbose = verbose
                )
            xyz_string = opt_mol.tostring(format="xyz")

            return {
                "converged": True,
                "xyz_data": xyz_string
            }
        
        except Exception as error:
            return {
                "converged": False,
                "xyz_data": None,
                "error": str(error)
            }
        