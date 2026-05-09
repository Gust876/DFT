from abc import ABC, abstractmethod


class FactoryOptimizer(ABC):
    '''
    Interface abstrata para factories de otimizadores geométricos (padrão Factory Method).

    Cada engine implementa sua própria factory para construir o
    otimizador de geometria molecular correspondente.
    '''

    @abstractmethod
    def factory_method(self):
        '''
        Cria e retorna uma instância do produto otimizador.

        Returns:
            OptimizerProduct: Objeto responsável por realizar a otimização geométrica.
        '''
        pass


class OptimizerProduct(ABC):
    '''
    Interface abstrata para produtos de otimização geométrica.

    Define o contrato para a execução da otimização geométrica
    de moléculas usando diferentes engines de cálculo DFT.
    '''

    @abstractmethod
    def opt_geometry(self):
        '''
        Executa a otimização geométrica da molécula.

        Returns:
            dict: Dicionário contendo:
                - 'converged' (bool): True se a otimização convergiu.
                - 'xyz_data' (str): Geometria otimizada em formato XYZ.
                - 'error' (str, opcional): Mensagem de erro caso não haja convergência.
        '''
        pass

