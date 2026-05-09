from interfaces.result_interface import EletronicStructure
from pyscf.tools import molden, cubegen
from pathlib import Path


class PySCFResult(EletronicStructure):
    '''
    Resultado de cálculo DFT de frequências vibracionais via PySCF.

    Encapsula os dados produzidos após o cálculo da Hessiana e análise
    harmônica, fornecendo métodos para exportar a estrutura eletrônica
    nos formatos padrão de química computacional (.molden e .cube).

    Attributes:
        _frequencies (ndarray): Array com as frequências vibracionais
            em cm⁻¹. Frequências puramente reais indicam estado fundamental.
        mol (pyscf.gto.Mole): Objeto de molécula do PySCF.
        mf (pyscf.dft.RKS): Objeto de campo médio com os resultados SCF.
    '''

    def __init__(self, frequencies, mol, mf):
        '''
        Inicializa o resultado com os dados do cálculo PySCF.

        Args:
            frequencies (ndarray): Array de frequências vibracionais em cm⁻¹.
            mol (pyscf.gto.Mole): Objeto de molécula do PySCF.
            mf (pyscf.dft.RKS): Objeto de campo médio com resultados SCF.
        '''
        self._frequencies = frequencies
        self.mol = mol
        self.mf = mf

    def write_molden(self, path: Path):
        '''
        Exporta os orbitais moleculares no formato Molden.

        Gera um arquivo .molden contendo o cabeçalho da molécula e
        os coeficientes dos orbitais moleculares, energias e números
        de ocupação, compatível com visualizadores como Avogadro e VESTA.

        Args:
            path (Path): Caminho completo do arquivo de saída (.molden).
        '''

        with open(path, "w") as file:
            molden.header(self.mol, file)
            molden.orbital_coeff(
                self.mol,
                file,
                mo_coeff=self.mf.mo_coeff,
                ene=self.mf.mo_energy,
                occ=self.mf.mo_occ
            )

    def write_density_cube(self, path: Path):
        '''
        Exporta a densidade eletrônica no formato Gaussian Cube.

        Gera um arquivo .cube com a densidade eletrônica calculada a
        partir da matriz densidade de um elétron (RDM1), compatível
        com visualizadores volumétricos como VMD e VESTA.

        Args:
            path (Path): Caminho completo do arquivo de saída (.cube).
        '''
        cubegen.density(
            self.mol,
            str(path),
            self.mf.make_rdm1()
        )

    @property
    def frequencies(self):
        '''
        Array de frequências vibracionais em cm⁻¹.

        Returns:
            ndarray: Frequências vibracionais. Valores puramente reais
                indicam que a geometria corresponde a um mínimo de energia.
        '''
        return self._frequencies
