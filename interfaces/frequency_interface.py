from abc import ABC, abstractmethod


class FactoryFrequency(ABC):

    @abstractmethod
    def factory_method(self):
        pass


class FrequencyProduct(ABC):

    @abstractmethod
    def vibrational_frequency(self):
        pass