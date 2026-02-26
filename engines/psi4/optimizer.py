from interfaces.optimizer_interface import FactoryOptimizer, OptimizerProduct
from engines.psi4.mol import FactoryMolPsi4 
from engines.utils import factory_mol
from pathlib import Path
import psi4


class FactoryOptimizerPsi4(FactoryOptimizer):

    def factory_method(self):
        return OptimizerPsi4()
    

class OptimizerPsi4(OptimizerProduct):
     
     def opt_geometry(self, xyz_file: Path):
        psi4.set_memory("500 MB")
        
        try:
            construct_mol = factory_mol(FactoryMolPsi4())
            mol = construct_mol.create_mol(
                xyz_file=xyz_file.read_text()
            )

            psi4.set_options({"reference": "rhf"})
            psi4.optimize("scf/cc-pvdz", molecule=mol)

            xyz_string = f"{mol.natom()} \n {mol.save_string_xyz()}"
            return {
                "converged": True,
                "xyz_data": xyz_string
            }

        except Exception as error:
            return {
                "converged": False,
                "error": str(error)
            }
        