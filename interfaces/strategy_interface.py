from abc import ABC, abstractmethod


class EngineStrategy(ABC):

    @abstractmethod
    def construct_mol(self, *args):
        pass

    @abstractmethod
    def optimizer(self, xyz_file, **kwargs):
        pass

    @abstractmethod
    def frequency(self, xyz_file, **kwargs):
        pass
    