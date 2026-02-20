from pyscf import gto
from pyscf.gto import Mole


def create_molecule(xyz_file: str, basis: str) -> Mole:
    mol = gto.M(
        atom = xyz_file,
        basis = basis
    )
    return mol
