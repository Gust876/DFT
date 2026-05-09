from interfaces.strategy_interface import EngineStrategy
from engines.utils import factory_mol, factory_optimizer, factory_frequency
from engines.pyscf.optimizer import FactoryOptimizerPySCF
from engines.pyscf.frequency import FactoryFrequencyPySCF
from engines.pyscf.mol import FactoryMolPySCF
from engines.psi4.optimizer import FactoryOptimizerPsi4
from engines.psi4.frequency import FactoryFrequencyPsi4
from engines.psi4.mol import FactoryMolPsi4
from pathlib import Path


class PySCFStrategy(EngineStrategy):
    '''
    Estratégia de execução DFT utilizando a engine PySCF.

    Implementa o padrão Strategy para a engine PySCF, orquestrando
    a construção da molécula e a execução dos cálculos de otimização
    geométrica e frequências vibracionais. A construção explícita do
    objeto Mole é necessária pois o PySCF recebe a molécula como
    argumento dos calculadores, diferentemente do Psi4.
    '''

    def construct_mol(self, xyz_file: Path, basis: str):
        '''
        Constrói o objeto Mole do PySCF a partir de um arquivo XYZ.

        Args:
            xyz_file (Path): Caminho para o arquivo XYZ com coordenadas atômicas.
            basis (str): Conjunto de funções de base (ex: 'cc-pvdz', '6-311g**').

        Returns:
            pyscf.gto.Mole: Objeto de molécula configurado para uso no PySCF.
        '''
        construct_mol = factory_mol(FactoryMolPySCF())
        mol = construct_mol.create_mol(
            xyz_file=xyz_file.read_text(),
            basis=basis
        )
        return mol
    
    def optimizer(self,  xyz_file: Path, **kwargs):
        '''
        Executa a otimização geométrica com PySCF.

        Constrói a molécula, inicializa o otimizador geométrico via
        factory e executa a otimização usando o solver geométrico
        do PySCF (geomeTRIC).

        Args:
            xyz_file (Path): Caminho para o arquivo XYZ com a geometria inicial.
            **kwargs:
                basis (str): Conjunto de funções de base.
                xc (str): Funcional de troca-correlação (ex: 'm06-2x', 'b3lyp').

        Returns:
            dict: Dicionário com chaves 'converged' (bool), 'xyz_data' (str)
                e 'error' (str, opcional).
        '''
        mol = self.construct_mol(xyz_file, kwargs['basis'])

        optimizer = factory_optimizer(
            FactoryOptimizerPySCF(
                mol=mol,
                xc=kwargs['xc']
            )
        )
        dict_optimize = optimizer.opt_geometry(
            maxsteps=200,
            verbose=0
        )
        return dict_optimize    
    
    def frequency(self, xyz_file: Path, **kwargs):
        '''
        Executa o cálculo de frequências vibracionais com PySCF.

        Constrói a molécula, inicializa o calculador de frequências
        via factory e executa o cálculo da Hessiana e análise harmônica.

        Args:
            xyz_file (Path): Caminho para o arquivo XYZ com a geometria otimizada.
            **kwargs:
                basis (str): Conjunto de funções de base.
                xc (str): Funcional de troca-correlação.

        Returns:
            PySCFResult: Objeto com frequências vibracionais e estrutura eletrônica.
        '''
        mol = self.construct_mol(xyz_file, kwargs['basis'])

        freq_obj = factory_frequency(
            FactoryFrequencyPySCF(
                mol=mol,
                xc=kwargs['xc']
            )
        )
        return freq_obj.vibrational_frequency()
    

class Psi4Strategy(EngineStrategy):
    '''
    Estratégia de execução DFT utilizando a engine Psi4.

    Implementa o padrão Strategy para a engine Psi4, orquestrando
    a execução dos cálculos de otimização geométrica e frequências
    vibracionais. O Psi4 gerencia internamente a leitura do arquivo
    XYZ, de modo que a construção explícita da molécula não é
    necessária nos métodos de alto nível.
    '''

    def construct_mol(self, xyz_file):
        '''
        Constrói o objeto de molécula do Psi4 a partir de um arquivo XYZ.

        Args:
            xyz_file (str | Path): Conteúdo ou caminho do arquivo XYZ.

        Returns:
            psi4.core.Molecule: Objeto de molécula configurado para uso no Psi4.
        '''
        construct_mol = factory_mol(FactoryMolPsi4())
        mol = construct_mol.create_mol(xyz_file)
        return mol
    
    def optimizer(self, xyz_file: Path, **kwargs):
        '''
        Executa a otimização geométrica com Psi4.

        Args:
            xyz_file (Path): Caminho para o arquivo XYZ com a geometria inicial.
            **kwargs:
                basis (str): Conjunto de funções de base.
                xc (str): Método de cálculo (ex: 'b3lyp', 'm06-2x').

        Returns:
            dict: Dicionário com chaves 'converged' (bool), 'xyz_data' (str)
                e 'error' (str, opcional).
        '''
        optimizer = factory_optimizer(FactoryOptimizerPsi4())
        dict_optimizer = optimizer.opt_geometry(
            xyz_file=xyz_file,
            xc=kwargs['xc'],
            basis=kwargs['basis']
        )
        return dict_optimizer
    
    def frequency(self, xyz_file, **kwargs):
        '''
        Executa o cálculo de frequências vibracionais com Psi4.

        Args:
            xyz_file (Path): Caminho para o arquivo XYZ com a geometria otimizada.
            **kwargs:
                basis (str): Conjunto de funções de base.
                xc (str): Método de cálculo.

        Returns:
            Psi4Result: Objeto com frequências vibracionais e estrutura eletrônica.
        '''
        freq_obj = factory_frequency(
            FactoryFrequencyPsi4(
                xyz_file=xyz_file,
                xc=kwargs['xc'],
                basis=kwargs['basis']
            )
        )
        return freq_obj.vibrational_frequency()
    