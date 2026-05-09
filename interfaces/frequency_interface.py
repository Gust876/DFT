from abc import ABC, abstractmethod


class FactoryFrequency(ABC):
    '''
    Interface abstrata para factories de cálculo de frequências vibracionais (padrão Factory Method).

    Cada engine implementa sua própria factory para construir o
    calculador de frequências vibracionais correspondente.
    '''

    @abstractmethod
    def factory_method(self):
        '''
        Cria e retorna uma instância do produto de frequência.

        Returns:
            FrequencyProduct: Objeto responsável por calcular as frequências vibracionais.
        '''
        pass


class FrequencyProduct(ABC):
    '''
    Interface abstrata para produtos de cálculo de frequências vibracionais.

    Define o contrato para o cálculo de frequências vibracionais
    de moléculas em estado fundamental usando diferentes engines DFT.
    A presença exclusiva de frequências reais (positivas) confirma
    que a geometria corresponde a um mínimo de energia (estado fundamental).
    '''

    @abstractmethod
    def vibrational_frequency(self):
        '''
        Executa o cálculo de frequências vibracionais.

        Returns:
            EletronicStructure: Objeto de resultado contendo frequências e
                informações necessárias para geração de arquivos de saída.
        '''
        pass