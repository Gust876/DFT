from interfaces.strategy_interface import EngineStrategy
from engines.utils import factory_mol, factory_optimizer, factory_frequency
from engines.pyscf.optimizer import FactoryOptimizerPySCF
from engines.pyscf.frequency import FactoryFrequencyPySCF
from engines.pyscf.mol import FactoryMolPySCF
from engines.psi4.optimizer import FactoryOptimizerPsi4
from engines.psi4.frequency import FactoryFrequencyPsi4
from engines.psi4.mol import FactoryMolPsi4
from pathlib import Path


class PySCFStrategy(EngineStrategy):

    def construct_mol(self, xyz_file: Path, basis: str):
        construct_mol = factory_mol(FactoryMolPySCF())
        mol = construct_mol.create_mol(
            xyz_file=xyz_file.read_text(),
            basis=basis
        )
        return mol
    
    def optimizer(self,  xyz_file: Path, **kwargs):
        mol = self.construct_mol(xyz_file, kwargs['basis'])

        optimizer = factory_optimizer(
            FactoryOptimizerPySCF(
                mol=mol,
                xc=kwargs['xc']
            )
        )
        dict_optimize = optimizer.opt_geometry(
            maxsteps=200,
            verbose=0
        )
        return dict_optimize    
    
    def frequency(self, xyz_file: Path, **kwargs):
        mol = self.construct_mol(xyz_file, kwargs['basis'])

        freq_obj = factory_frequency(
            FactoryFrequencyPySCF(
                mol=mol,
                xc=kwargs['xc']
            )
        )
        return freq_obj.vibrational_frequency()
    

class Psi4Strategy(EngineStrategy):

    def construct_mol(self, xyz_file):
        construct_mol = factory_mol(FactoryMolPsi4())
        mol = construct_mol.create_mol(xyz_file)
        return mol
    
    def optimizer(self, xyz_file: Path, **kwargs):
        optimizer = factory_optimizer(FactoryOptimizerPsi4())
        dict_optimizer = optimizer.opt_geometry(
            xyz_file=xyz_file,
            xc=kwargs['xc'],
            basis=kwargs['basis']
        )
        return dict_optimizer
    
    def frequency(self, xyz_file, **kwargs):
        freq_obj = factory_frequency(
            FactoryFrequencyPsi4(
                xyz_file=xyz_file,
                xc=kwargs['xc'],
                basis=kwargs['basis']
            )
        )
        return freq_obj.vibrational_frequency()
    