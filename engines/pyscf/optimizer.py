from pyscf.geomopt.geometric_solver import optimize
from interfaces.engine_interface import EngineDFT
from interfaces.optimizer_interface import FactoryOptimizer
from interfaces.optimizer_interface import OptimizerProduct
from engines.pyscf.engine_pyscf import EnginePySCF
from pyscf.gto import Mole


class FactoryOptimizerPySCF(FactoryOptimizer):
    '''
    Factory de otimizadores geométricos para a engine PySCF (padrão Factory Method).
 
    Recebe a molécula e o funcional configurados e cria a instância
    do otimizador geométrico correspondente, inicializando a engine
    DFT internamente.
 
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
        Cria e retorna uma instância de OptimizerPySCF.
 
        Returns:
            OptimizerPySCF: Otimizador geométrico configurado com a engine PySCF.
        '''
        return OptimizerPySCF(
            engine=self._initialize_engine()
            )
    

class OptimizerPySCF(OptimizerProduct):
    '''
    Otimizador geométrico para a engine PySCF.
 
    Executa a otimização da geometria molecular usando o solver
    geométrico geomeTRIC integrado ao PySCF, retornando a geometria
    otimizada em formato XYZ.
 
    Attributes:
        _engine (EngineDFT): Engine DFT configurada para o cálculo.
    '''
    def __init__(self, engine: EngineDFT):
        '''
        Inicializa o otimizador com a engine DFT configurada.
 
        Args:
            engine (EngineDFT): Instância da engine DFT (EnginePySCF).
        '''
        self._engine = engine

    def opt_geometry(self, maxsteps: int, verbose: int):
        '''
        Executa a otimização geométrica da molécula com PySCF.
 
        Utiliza o solver geomeTRIC para minimizar a energia em relação
        às coordenadas atômicas, retornando a geometria otimizada em
        formato XYZ caso a convergência seja atingida.
 
        Args:
            maxsteps (int): Número máximo de passos de otimização.
            verbose (int): Nível de verbosidade do otimizador (0 = silencioso).
 
        Returns:
            dict: Dicionário contendo:
                - 'converged' (bool): True se a otimização convergiu.
                - 'xyz_data' (str): Geometria otimizada em formato XYZ.
                - 'error' (str, opcional): Mensagem de erro se não convergiu.
        '''
        try:
            mf = self._engine.dft_method()

            opt_mol = optimize(
                mf,
                maxsteps = maxsteps,
                verbose = verbose
                )
            xyz_string = opt_mol.tostring(format="xyz")

            return {
                "converged": True,
                "xyz_data": xyz_string
            }
        
        except Exception as error:
            return {
                "converged": False,
                "xyz_data": None,
                "error": str(error)
            }
        