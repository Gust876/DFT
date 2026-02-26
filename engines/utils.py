from interfaces.optimizer_interface import FactoryOptimizer, OptimizerProduct
from interfaces.frequency_interface import FactoryFrequency, FrequencyProduct
from interfaces.molecule_interface import MolFactory, MolProduct


def factory_mol(factory: MolFactory) -> MolProduct:
    return factory.factory_method()

def factory_optimizer(factory: FactoryOptimizer) -> OptimizerProduct:
    return factory.factory_method()

def factory_frequency(factory: FactoryFrequency) -> FrequencyProduct:
    return factory.factory_method()
