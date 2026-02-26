from abc import ABC, abstractmethod


class FactoryOptimizer(ABC):

    @abstractmethod
    def factory_method(self):
        pass


class OptimizerProduct(ABC):

    @abstractmethod
    def opt_geometry(self):
        pass

