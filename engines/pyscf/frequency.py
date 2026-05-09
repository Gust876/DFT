from interfaces.frequency_interface import FactoryFrequency
from interfaces.frequency_interface import FrequencyProduct
from interfaces.result_interface import EletronicStructure
from interfaces.engine_interface import EngineDFT
from engines.pyscf.engine_pyscf import EnginePySCF
from results.pyscf_result import PySCFResult
from pyscf.hessian import thermo
from pyscf.gto import Mole


class FactoryFrequencyPySCF(FactoryFrequency):
    '''
    Factory de calculadores de frequências vibracionais para PySCF (padrão Factory Method).
 
    Recebe a molécula e o funcional configurados e cria a instância
    do calculador de frequências, inicializando a engine DFT internamente.
 
    Attributes:
        mol (pyscf.gto.Mole): Objeto de molécula do PySCF.
        xc (str): Funcional de troca-correlação.
    '''

    def __init__(self, mol: Mole, xc: str):
        '''
        Inicializa a factory com os parâmetros do cálculo.
 
        Args:
            mol (pyscf.gto.Mole): Objeto de molécula já construído.
            xc (str): Funcional de troca-correlação (ex: 'm06-2x').
        '''
        self.mol = mol
        self.xc = xc

    def _initialize_engine(self):
        '''
        Inicializa e retorna a engine DFT configurada.
 
        Returns:
            EnginePySCF: Engine DFT configurada com a molécula e o funcional.
        '''
        return EnginePySCF(mol=self.mol, xc=self.xc)

    def factory_method(self):
        '''
        Cria e retorna uma instância de FrequencyPySCF.
 
        Returns:
            FrequencyPySCF: Calculador de frequências configurado com a engine PySCF.
        '''
        return FrequencyPySCF(
            engine=self._initialize_engine()
        )


class FrequencyPySCF(FrequencyProduct):
    '''
    Calculador de frequências vibracionais para a engine PySCF.
 
    Executa o cálculo SCF, computa a Hessiana analítica e realiza
    a análise harmônica para obtenção das frequências vibracionais,
    retornando os resultados encapsulados em um PySCFResult.
 
    Attributes:
        _engine (EngineDFT): Engine DFT configurada para o cálculo.
    '''

    def __init__(self, engine: EngineDFT) -> EletronicStructure:
        '''
        Inicializa o calculador com a engine DFT configurada.
 
        Args:
            engine (EngineDFT): Instância da engine DFT (EnginePySCF).
        '''
        self._engine = engine

    def vibrational_frequency(self):
        '''
        Executa o cálculo de frequências vibracionais com PySCF.
 
        Realiza o cálculo SCF, computa a Hessiana analítica via
        método de diferenças finitas e executa a análise harmônica
        para obter as frequências vibracionais em cm⁻¹.
 
        Returns:
            PySCFResult: Objeto contendo as frequências vibracionais e
                os dados necessários para geração dos arquivos .molden e .cube.
        '''
        mf = self._engine.dft_method()
        mf.run()
        
        hessian = mf.Hessian().kernel()
        freq_info = thermo.harmonic_analysis(mf.mol, hessian)
        frequencies = freq_info["freq_wavenumber"]

        return PySCFResult(frequencies, mf.mol, mf)