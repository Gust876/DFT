from interfaces.optimizer_interface import FactoryOptimizer, OptimizerProduct
from interfaces.frequency_interface import FactoryFrequency, FrequencyProduct
from interfaces.molecule_interface import MolFactory, MolProduct


def factory_mol(factory: MolFactory) -> MolProduct:
    '''
    Instancia o produto de molécula a partir de uma factory (padrão Factory Method).
 
    Args:
        factory (MolFactory): Factory concreta da engine desejada
            (ex: FactoryMolPySCF, FactoryMolPsi4).
 
    Returns:
        MolProduct: Objeto responsável por construir a molécula no formato da engine.
    '''
    return factory.factory_method()

def factory_optimizer(factory: FactoryOptimizer) -> OptimizerProduct:
    '''
    Instancia o otimizador geométrico a partir de uma factory (padrão Factory Method).
 
    Args:
        factory (FactoryOptimizer): Factory concreta da engine desejada
            (ex: FactoryOptimizerPySCF, FactoryOptimizerPsi4).
 
    Returns:
        OptimizerProduct: Objeto responsável por realizar a otimização geométrica.
    '''
    return factory.factory_method()

def factory_frequency(factory: FactoryFrequency) -> FrequencyProduct:
    '''
    Instancia o calculador de frequências a partir de uma factory (padrão Factory Method).
 
    Args:
        factory (FactoryFrequency): Factory concreta da engine desejada
            (ex: FactoryFrequencyPySCF, FactoryFrequencyPsi4).
 
    Returns:
        FrequencyProduct: Objeto responsável por calcular as frequências vibracionais.
    '''
    return factory.factory_method()
