from interfaces.frequency_interface import FactoryFrequency, FrequencyProduct
from results.psi4_result import Psi4Result
from engines.psi4.mol import FactoryMolPsi4
from engines.utils import factory_mol
from pathlib import Path
import psi4


class FactoryFrequencyPsi4(FactoryFrequency):

    def __init__(self, xyz_file: Path):
        self.xyz_file = xyz_file

    def factory_method(self):
        return FrequencyPsi4(xyz_file=self.xyz_file)
    

class FrequencyPsi4(FrequencyProduct):

    def __init__(self, xyz_file: Path):
        self.xyz_file = xyz_file

    def vibrational_frequency(self):
        psi4.set_memory("500 MB")

        construct_mol = factory_mol(FactoryMolPsi4())
        mol = construct_mol.create_mol(
            xyz_file=self.xyz_file.read_text()
        )

        scf_e, scf_wfn = psi4.frequency(
            "scf/cc-pvdz",
            molecule=mol,
            return_wfn=True
        )
        frequencies = scf_wfn.frequencies().to_array()

        return Psi4Result(
            frequencies=frequencies,
            wfn=scf_wfn
        )
        