from interfaces.frequency_interface import FactoryFrequency
from interfaces.frequency_interface import FrequencyProduct
from interfaces.result_interface import EletronicStructure
from interfaces.engine_interface import EngineDFT
from engines.pyscf.engine_pyscf import EnginePySCF
from results.pyscf_result import PySCFResult
from pyscf.hessian import thermo
from pyscf.gto import Mole


class FactoryFrequencyPySCF(FactoryFrequency):

    def __init__(self, mol: Mole, xc: str):
        self.mol = mol
        self.xc = xc

    def _initialize_engine(self):
        return EnginePySCF(mol=self.mol, xc=self.xc)

    def factory_method(self):
        return FrequencyPySCF(
            engine=self._initialize_engine()
        )


class FrequencyPySCF(FrequencyProduct):

    def __init__(self, engine: EngineDFT) -> EletronicStructure:
        self._engine = engine

    def vibrational_frequency(self):
        mf = self._engine.dft_method()
        mf.run()
        
        hessian = mf.Hessian().kernel()
        freq_info = thermo.harmonic_analysis(mf.mol, hessian)
        frequencies = freq_info["freq_wavenumber"]

        return PySCFResult(frequencies, mf.mol, mf)