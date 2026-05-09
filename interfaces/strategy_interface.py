from abc import ABC, abstractmethod


class EngineStrategy(ABC):
    '''
    Interface abstrata para estratégias de execução de engines DFT (padrão Strategy).

    Define o contrato que todas as estratégias de engine devem seguir,
    permitindo que os workflows sejam completamente desacoplados das
    implementações específicas de PySCF e Psi4. Novas engines podem
    ser integradas ao pipeline implementando esta interface.
    '''

    @abstractmethod
    def construct_mol(self, *args):
        '''
        Constrói o objeto de molécula específico da engine.

        Args:
            *args: Argumentos necessários para construção da molécula
                (ex: caminho do arquivo XYZ, basis set).

        Returns:
            Objeto de molécula no formato da engine correspondente.
        '''
        pass

    @abstractmethod
    def optimizer(self, xyz_file, **kwargs):
        '''
        Executa a otimização geométrica usando a engine configurada.

        Args:
            xyz_file (Path): Caminho para o arquivo XYZ com a geometria inicial.
            **kwargs: Parâmetros do cálculo:
                - basis (str): Conjunto de funções de base (ex: 'cc-pvdz').
                - xc (str): Funcional de troca-correlação ou método (ex: 'm06-2x').

        Returns:
            dict: Resultado da otimização com chaves 'converged', 'xyz_data' e 'error'.
        '''
        pass

    @abstractmethod
    def frequency(self, xyz_file, **kwargs):
        '''
        Executa o cálculo de frequências vibracionais usando a engine configurada.

        Args:
            xyz_file (Path): Caminho para o arquivo XYZ com a geometria otimizada.
            **kwargs: Parâmetros do cálculo:
                - basis (str): Conjunto de funções de base.
                - xc (str): Funcional de troca-correlação ou método.

        Returns:
            EletronicStructure: Resultado do cálculo com frequências e estrutura eletrônica.
        '''
        pass
    