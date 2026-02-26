from abc import ABC, abstractmethod


class EletronicStructure(ABC):

    @abstractmethod
    def write_molden(self):
        pass

    @abstractmethod
    def write_density_cube(self):
        pass
