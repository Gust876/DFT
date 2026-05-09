from interfaces.strategy_interface import EngineStrategy


class SetOptimizerStrategy:
    '''
    Contexto do padrão Strategy para execução de otimização geométrica.

    Encapsula uma estratégia de engine (PySCF ou Psi4) e delega
    a execução da otimização geométrica à estratégia configurada,
    mantendo os workflows completamente desacoplados das implementações
    específicas de cada engine.

    Attributes:
        _strategy (EngineStrategy): Estratégia de engine configurada.
    '''
    def __init__(self, strategy: EngineStrategy):
        '''
        Inicializa o contexto com a estratégia de engine desejada.

        Args:
            strategy (EngineStrategy): Instância da estratégia de engine
                (ex: PySCFStrategy, Psi4Strategy).
        '''
        self._strategy = strategy
    
    def construct_mol(self, *args):
        '''
        Delega a construção da molécula à estratégia configurada.

        Args:
            *args: Argumentos passados diretamente à estratégia.

        Returns:
            Objeto de molécula no formato da engine correspondente.
        '''
        return self._strategy.construct_mol(*args)
    
    def optimizer(self, xyz_file, **kwargs):
        '''
        Delega a otimização geométrica à estratégia configurada.

        Args:
            xyz_file (Path): Caminho para o arquivo XYZ com a geometria inicial.
            **kwargs: Parâmetros do cálculo (basis, xc).

        Returns:
            dict: Resultado da otimização com chaves 'converged', 'xyz_data' e 'error'.
        '''
        return self._strategy.optimizer(xyz_file, **kwargs)
    

class SetFrequencyStrategy:
    '''
    Contexto do padrão Strategy para execução de cálculo de frequências vibracionais.

    Encapsula uma estratégia de engine (PySCF ou Psi4) e delega
    o cálculo de frequências à estratégia configurada, permitindo
    que o workflow de frequências opere independentemente da engine escolhida.

    Attributes:
        _strategy (EngineStrategy): Estratégia de engine configurada.
    '''
    def __init__(self, strategy: EngineStrategy):
        '''
        Inicializa o contexto com a estratégia de engine desejada.

        Args:
            strategy (EngineStrategy): Instância da estratégia de engine
                (ex: PySCFStrategy, Psi4Strategy).
        '''
        self._strategy = strategy
    
    def construct_mol(self, *args):
        '''
        Delega a construção da molécula à estratégia configurada.

        Args:
            *args: Argumentos passados diretamente à estratégia.

        Returns:
            Objeto de molécula no formato da engine correspondente.
        '''
        return self._strategy.construct_mol(*args)
    
    def frequency(self, xyz_file, **kwargs):
        '''
        Delega o cálculo de frequências vibracionais à estratégia configurada.

        Args:
            xyz_file (Path): Caminho para o arquivo XYZ com a geometria otimizada.
            **kwargs: Parâmetros do cálculo (basis, xc).

        Returns:
            EletronicStructure: Resultado com frequências e estrutura eletrônica.
        '''
        return self._strategy.frequency(xyz_file, **kwargs)
    