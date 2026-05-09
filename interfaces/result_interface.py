from abc import ABC, abstractmethod


class EletronicStructure(ABC):
    '''
    Interface abstrata para resultados de estrutura eletrônica.

    Define o contrato para objetos que encapsulam os resultados de
    cálculos DFT de frequências vibracionais, fornecendo métodos
    para exportar a estrutura eletrônica nos formatos padrão
    utilizados em química computacional.
    '''

    @abstractmethod
    def write_molden(self):
        '''
        Exporta a estrutura eletrônica no formato Molden.

        O formato Molden é amplamente utilizado para visualização
        de orbitais moleculares em softwares como Molden, Avogadro
        e VESTA.

        Args:
            path (Path): Caminho completo do arquivo de saída (.molden).
        '''
        pass

    @abstractmethod
    def write_density_cube(self):
        '''
        Exporta a densidade eletrônica no formato Gaussian Cube.

        O formato Cube permite a visualização volumétrica da densidade
        eletrônica em softwares como VMD, VESTA e Avogadro.

        Args:
            path (Path): Caminho completo do arquivo de saída (.cube).
        '''
        pass
