from interfaces.frequency_interface import FactoryFrequency, FrequencyProduct
from results.psi4_result import Psi4Result
from engines.psi4.mol import FactoryMolPsi4
from engines.utils import factory_mol
from pathlib import Path
import psi4


class FactoryFrequencyPsi4(FactoryFrequency):
    '''
    Factory de calculadores de frequências vibracionais para Psi4 (padrão Factory Method).
 
    Recebe os parâmetros do cálculo e cria a instância do calculador
    de frequências para a engine Psi4. Diferentemente do PySCF, o
    gerenciamento da molécula é feito internamente pelo calculador.
 
    Attributes:
        xyz_file (Path): Caminho para o arquivo XYZ com a geometria otimizada.
        xc (str): Método de cálculo.
        basis (str): Conjunto de funções de base.
    '''

    def __init__(self, xyz_file: Path, xc: str, basis: str):
        '''
        Inicializa a factory com os parâmetros do cálculo.
 
        Args:
            xyz_file (Path): Caminho para o arquivo XYZ com a geometria otimizada.
            xc (str): Método de cálculo (ex: 'b3lyp', 'm06-2x').
            basis (str): Conjunto de funções de base (ex: 'cc-pvdz').
        '''
        self.xyz_file = xyz_file
        self.xc = xc
        self.basis = basis

    def factory_method(self):
        '''
        Cria e retorna uma instância de FrequencyPsi4.
 
        Returns:
            FrequencyPsi4: Calculador de frequências configurado para a engine Psi4.
        '''
        return FrequencyPsi4(
            xyz_file=self.xyz_file,
            xc=self.xc,
            basis=self.basis
        )
    

class FrequencyPsi4(FrequencyProduct):
    '''
    Calculador de frequências vibracionais para a engine Psi4.
 
    Executa o cálculo de frequências vibracionais usando o módulo
    de frequências do Psi4, retornando os resultados encapsulados
    em um Psi4Result com acesso à função de onda para geração
    dos arquivos de saída.
 
    Attributes:
        xyz_file (Path): Caminho para o arquivo XYZ com a geometria otimizada.
        xc (str): Método de cálculo.
        basis (str): Conjunto de funções de base.
    '''

    def __init__(self, xyz_file: Path, xc: str, basis: str):
        '''
        Inicializa o calculador com os parâmetros do cálculo.
 
        Args:
            xyz_file (Path): Caminho para o arquivo XYZ com a geometria otimizada.
            xc (str): Método de cálculo (ex: 'b3lyp', 'm06-2x').
            basis (str): Conjunto de funções de base (ex: 'cc-pvdz').
        '''
        self.xyz_file = xyz_file
        self.xc = xc
        self.basis = basis

    def vibrational_frequency(self):
        '''
        Executa o cálculo de frequências vibracionais com Psi4.
 
        Constrói a molécula a partir do arquivo XYZ, executa o cálculo
        de frequências vibracionais e retorna o resultado encapsulado
        com acesso à função de onda para exportação dos arquivos de saída.
 
        Returns:
            Psi4Result: Objeto contendo as frequências vibracionais e a
                função de onda necessária para geração dos arquivos .molden e .cube.
        '''
        psi4.set_memory("500 MB")

        construct_mol = factory_mol(FactoryMolPsi4())
        mol = construct_mol.create_mol(
            xyz_file=self.xyz_file.read_text()
        )

        scf_e, scf_wfn = psi4.frequency(
            f"{self.xc}/{self.basis}",
            molecule=mol,
            return_wfn=True
        )
        frequencies = scf_wfn.frequencies().to_array()

        return Psi4Result(
            frequencies=frequencies,
            wfn=scf_wfn
        )
        