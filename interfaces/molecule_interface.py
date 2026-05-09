from abc import ABC, abstractmethod


class MolFactory(ABC):
    '''
    Interface abstrata para factories de moléculas (padrão Factory Method).

    Cada engine implementa sua própria factory para construir objetos
    de molécula no formato esperado pela engine correspondente
    '''

    @abstractmethod
    def factory_method(self):
        '''
        Cria e retorna uma instância do produto de molécula.

        Returns:
            MolProduct: Objeto responsável por construir a molécula.
        '''
        pass


class MolProduct(ABC):
    '''
    Interface abstrata para produtos de molécula.

    Define o contrato para a criação de objetos de molécula
    a partir de arquivos de coordenadas no formato XYZ.
    '''

    @abstractmethod
    def create_mol(self):
        '''
        Constrói e retorna o objeto de molécula da engine.

        Args:
            *args: Argumentos posicionais específicos da engine.
            **kwargs: Argumentos nomeados específicos da engine.

        Returns:
            Objeto de molécula no formato da engine (ex: pyscf.gto.Mole, psi4.Molecule).
        '''
        pass