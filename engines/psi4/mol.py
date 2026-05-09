from interfaces.molecule_interface import MolProduct, MolFactory
import psi4


class FactoryMolPsi4(MolFactory):
    '''
    Factory de moléculas para a engine Psi4 (padrão Factory Method).
 
    Cria instâncias de MolPsi4, responsáveis por construir objetos
    psi4.core.Molecule a partir do conteúdo de arquivos XYZ.
    '''

    def factory_method(self):
        '''
        Cria e retorna uma instância de MolPsi4.
 
        Returns:
            MolPsi4: Objeto responsável por construir a molécula no formato Psi4.
        '''
        return MolPsi4()


class MolPsi4(MolProduct):
    '''
    Produto de molécula para a engine Psi4.
 
    Constrói objetos psi4.core.Molecule a partir do conteúdo de
    arquivos XYZ, utilizando a interface de geometria do Psi4.
    Diferentemente do PySCF, o Psi4 lê a molécula diretamente
    do conteúdo do arquivo, sem necessidade de parsing manual.
    '''

    def create_mol(self, xyz_file):
        '''
        Constrói e retorna um objeto Molecule do Psi4.
 
        Args:
            xyz_file (str): Conteúdo completo do arquivo XYZ como string,
                incluindo o número de átomos e o comentário nas primeiras linhas.
 
        Returns:
            psi4.core.Molecule: Objeto de molécula configurado para uso no Psi4.
        '''
        mol = psi4.geometry(xyz_file)
        return mol
    