from interfaces.molecule_interface import MolProduct, MolFactory
from pyscf.gto import Mole
from pyscf import gto


class FactoryMolPySCF(MolFactory):
    '''
    Factory de moléculas para a engine PySCF (padrão Factory Method).
 
    Cria instâncias de MolPySCF, responsáveis por construir objetos
    pyscf.gto.Mole a partir de arquivos de coordenadas no formato XYZ.
    '''
    
    def factory_method(self):
        '''
        Cria e retorna uma instância de MolPySCF.
 
        Returns:
            MolPySCF: Objeto responsável por construir a molécula no formato PySCF.
        '''
        return MolPySCF()


class MolPySCF(MolProduct):
    '''
    Produto de molécula para a engine PySCF.
 
    Constrói objetos pyscf.gto.Mole a partir do conteúdo de arquivos
    XYZ, extraindo as coordenadas atômicas e configurando o basis set
    e demais parâmetros necessários para os cálculos DFT.
    '''

    def create_mol(self, xyz_file: str, basis: str) -> Mole:
        '''
        Constrói e retorna um objeto Mole do PySCF.
 
        Lê as coordenadas atômicas do conteúdo do arquivo XYZ
        (ignorando as duas primeiras linhas de cabeçalho), configura
        o basis set e compila a molécula.
 
        Args:
            xyz_file (str): Conteúdo completo do arquivo XYZ como string.
            basis (str): Conjunto de funções de base (ex: 'cc-pvdz', '6-311g**').
 
        Returns:
            pyscf.gto.Mole: Objeto de molécula compilado e pronto para cálculos.
        '''

        atom_coordinates = "\n".join(xyz_file.strip().splitlines()[2:])

        mol = gto.Mole()
        mol.atom = atom_coordinates
        mol.basis = basis
        mol.max_memory = 2000
        mol.build()

        return mol
    