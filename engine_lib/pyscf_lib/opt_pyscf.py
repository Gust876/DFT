from pyscf.geomopt.geometric_solver import optimize
from engine_lib.engine_interface import EngineDFT


class OptimizerPySCF:
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
            
            xyz_string = opt_mol.tostring(format='xyz')
            return {
                "converged": True,
                "opt_mol": opt_mol,
                "xyz_data": xyz_string
            }
        
        except Exception as error:
            return {
                "converged": False,
                "xyz_data": None,
                "error": str(error)
            }

