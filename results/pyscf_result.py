from interfaces.result_interface import EletronicStructure
from pyscf.tools import molden, cubegen
from pathlib import Path


class PySCFResult(EletronicStructure):

    def __init__(self, frequencies, mol, mf):
        self._frequencies = frequencies
        self.mol = mol
        self.mf = mf

    def write_molden(self, path: Path):

        with open(path, "w") as file:
            molden.header(self.mol, file)
            molden.orbital_coeff(
                self.mol,
                file,
                mo_coeff=self.mf.mo_coeff,
                ene=self.mf.mo_energy,
                occ=self.mf.mo_occ
            )

    def write_density_cube(self, path: Path):
        cubegen.density(
            self.mol,
            str(path),
            self.mf.make_rdm1()
        )

    @property
    def frequencies(self):
        return self._frequencies
