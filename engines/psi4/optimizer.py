from interfaces.optimizer_interface import FactoryOptimizer, OptimizerProduct
from engines.psi4.mol import FactoryMolPsi4 
from engines.utils import factory_mol
from pathlib import Path
import psi4


class FactoryOptimizerPsi4(FactoryOptimizer):
    '''
    Factory de otimizadores geométricos para a engine Psi4 (padrão Factory Method).
 
    Cria instâncias de OptimizerPsi4, responsáveis por executar a
    otimização geométrica via Psi4. Diferentemente do PySCF, o Psi4
    gerencia internamente a molécula durante a otimização, de modo
    que esta factory não precisa receber parâmetros de inicialização
    '''

    def factory_method(self):
        '''
        Cria e retorna uma instância de OptimizerPsi4.
 
        Returns:
            OptimizerPsi4: Otimizador geométrico configurado para a engine Psi4.
        '''
        return OptimizerPsi4()
    

class OptimizerPsi4(OptimizerProduct):
     '''
     Otimizador geométrico para a engine Psi4.
 
    Executa a otimização da geometria molecular usando o otimizador
    interno do Psi4, retornando a geometria otimizada em formato XYZ.
    O Psi4 gerencia internamente a construção e atualização da
    molécula durante o processo de otimização.
     '''
     
     def opt_geometry(self, xyz_file: Path, xc: str, basis: str):
        '''
        Executa a otimização geométrica da molécula com Psi4.
 
        Constrói a molécula a partir do arquivo XYZ, configura a
        referência RHF e executa a otimização com o método e basis
        set especificados.
 
        Args:
            xyz_file (Path): Caminho para o arquivo XYZ com a geometria inicial.
            xc (str): Método de cálculo (ex: 'b3lyp', 'm06-2x').
            basis (str): Conjunto de funções de base (ex: 'cc-pvdz', '6-311g**').
 
        Returns:
            dict: Dicionário contendo:
                - 'converged' (bool): True se a otimização convergiu.
                - 'xyz_data' (str): Geometria otimizada em formato XYZ.
                - 'error' (str, opcional): Mensagem de erro se não convergiu.
        '''
        psi4.set_memory("500 MB")
        
        try:
            construct_mol = factory_mol(FactoryMolPsi4())
            mol = construct_mol.create_mol(
                xyz_file=xyz_file.read_text()
            )

            psi4.set_options({"reference": "rhf"})
            psi4.optimize(f"{xc}/{basis}", molecule=mol)

            xyz_string = f"{mol.natom()} \n {mol.save_string_xyz()}"
            return {
                "converged": True,
                "xyz_data": xyz_string
            }

        except Exception as error:
            return {
                "converged": False,
                "error": str(error)
            }
        