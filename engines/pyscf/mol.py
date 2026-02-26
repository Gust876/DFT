from interfaces.molecule_interface import MolProduct, MolFactory
from pyscf.gto import Mole
from pyscf import gto


class FactoryMolPySCF(MolFactory):
    
    def factory_method(self):
        return MolPySCF()


class MolPySCF(MolProduct):

    def create_mol(self, xyz_file: str, basis: str) -> Mole:

        atom_coordinates = "\n".join(xyz_file.strip().splitlines()[2:])

        mol = gto.M(
            atom = atom_coordinates,
            basis = basis
        )
        return mol
    