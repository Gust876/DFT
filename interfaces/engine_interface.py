from abc import ABC, abstractmethod

class EngineDFT(ABC):
    ''''
    Interface abstrata para engines de cálculo DFT.

    Define o contrato que todas as engines de cálculo de teoria
    do funcional da densidade (DFT) devem seguir. Implementações
    concretas devem fornecer o método DFT configurado e pronto
    para execução.
    '''

    @abstractmethod
    def dft_method(self):
        '''
        Retorna o objeto de método DFT configurado para a engine.

        Returns:
            Objeto de método DFT específico da engine (ex: pyscf.dft.RKS).
        '''
        pass
