from interfaces.engine_interface import EngineDFT
from pyscf.gto import Mole


class EnginePySCF(EngineDFT):
    '''
    Engine DFT concreta para cálculos com PySCF.
 
    Configura e retorna um objeto de método DFT usando a abordagem
    Kohn-Sham com ajuste de densidade (density fitting), adequada
    para cálculos de custo computacional reduzido mantendo boa
    precisão para funcionais híbridos e meta-GGA.
 
    Attributes:
        mol (pyscf.gto.Mole): Objeto de molécula do PySCF.
        xc (str): Funcional de troca-correlação (ex: 'm06-2x', 'b3lyp').
    '''
    def __init__(self, mol: Mole, xc: str):
        '''
        Inicializa a engine com a molécula e o funcional desejados.
 
        Args:
            mol (pyscf.gto.Mole): Objeto de molécula já construído e configurado.
            xc (str): Funcional de troca-correlação (ex: 'm06-2x', 'b3lyp').
        '''
        self.mol = mol
        self.xc = xc
    
    def dft_method(self):
        '''
        Configura e retorna o objeto de método DFT Kohn-Sham com density fitting.
 
        Aplica density fitting para redução do custo computacional,
        define o funcional de troca-correlação e configura parâmetros
        de convergência e memória.
 
        Returns:
            pyscf.dft.rks.RKS: Objeto de método DFT configurado e pronto
                para execução via .run() ou .kernel().
        '''
        mf = self.mol.KS().density_fit()
        mf.xc = self.xc
        
        mf.max_memory = 6000
        mf.conv_tol = 1e-6

        return mf