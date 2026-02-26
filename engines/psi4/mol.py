from interfaces.molecule_interface import MolProduct, MolFactory
import psi4


class FactoryMolPsi4(MolFactory):

    def factory_method(self):
        return MolPsi4()


class MolPsi4(MolProduct):

    def create_mol(self, xyz_file):
        mol = psi4.geometry(xyz_file)
        return mol
    