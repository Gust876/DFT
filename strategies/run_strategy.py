from interfaces.strategy_interface import EngineStrategy


class SetOptimizerStrategy:
    def __init__(self, strategy: EngineStrategy):
        self._strategy = strategy
    
    def construct_mol(self, *args):
        return self._strategy.construct_mol(*args)
    
    def optimizer(self, xyz_file, **kwargs):
        return self._strategy.optimizer(xyz_file, **kwargs)
    

class SetFrequencyStrategy:
    def __init__(self, strategy: EngineStrategy):
        self._strategy = strategy
    
    def construct_mol(self, *args):
        return self._strategy.construct_mol(*args)
    
    def frequency(self, xyz_file, **kwargs):
        return self._strategy.frequency(xyz_file, **kwargs)
    