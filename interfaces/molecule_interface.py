from abc import ABC, abstractmethod


class MolFactory(ABC):

    @abstractmethod
    def factory_method(self):
        pass


class MolProduct(ABC):

    @abstractmethod
    def create_mol(self):
        pass