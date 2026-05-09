from interfaces.result_interface import EletronicStructure
from pathlib import Path
import psi4


class Psi4Result(EletronicStructure):
    '''
    Resultado de cálculo DFT de frequências vibracionais via Psi4.

    Encapsula os dados produzidos após o cálculo de frequências no Psi4,
    fornecendo métodos para exportar a estrutura eletrônica nos formatos
    padrão de química computacional (.molden e .cube).

    Attributes:
        _frequencies (ndarray): Array com as frequências vibracionais em cm⁻¹.
        wfn (psi4.core.Wavefunction): Objeto de função de onda do Psi4
            contendo todos os dados eletrônicos do cálculo.
    '''

    def __init__(self, frequencies, wfn):
        '''
        Inicializa o resultado com os dados do cálculo Psi4.

        Args:
            frequencies (ndarray): Array de frequências vibracionais em cm⁻¹.
            wfn (psi4.core.Wavefunction): Função de onda resultante do cálculo.
        '''
        self._frequencies = frequencies
        self.wfn = wfn   
    
    def write_molden(self, path: Path):
        '''
        Exporta os orbitais moleculares no formato Molden via Psi4.

        Args:
            path (Path): Caminho completo do arquivo de saída (.molden).
        '''
        psi4.molden(
            wfn=self.wfn,
            filename=path
        )
    
    def write_density_cube(self, path: Path):
        '''
        Exporta a densidade eletrônica no formato Gaussian Cube via Psi4.

        Utiliza o módulo cubeprop do Psi4 para gerar o arquivo .cube
        da densidade eletrônica. O arquivo gerado pelo Psi4 é então
        renomeado para o caminho de saída especificado.

        Args:
            path (Path): Caminho completo do arquivo de saída (.cube).

        Note:
            O Psi4 gera o arquivo .cube no diretório de trabalho corrente.
            Esta implementação localiza e renomeia o arquivo automaticamente.
        '''
        psi4.set_options(
            {
                "cubeprop_tasks": ["density"]
            }
        )
        psi4.cubeprop(wfn=self.wfn)

        cube_file = list(Path(".").glob("*.cube"))
        cube_file[0].rename(path)

    @property
    def frequencies(self):
        '''
        Array de frequências vibracionais em cm⁻¹.

        Returns:
            ndarray: Frequências vibracionais. Valores puramente reais
                indicam que a geometria corresponde a um mínimo de energia.
        '''
        return self._frequencies
