from interfaces.result_interface import EletronicStructure
from pathlib import Path
import psi4


class Psi4Result(EletronicStructure):

    def __init__(self, frequencies, wfn):
        self._frequencies = frequencies
        self.wfn = wfn   
    
    def write_molden(self, path: Path):
        psi4.molden(
            wfn=self.wfn,
            filename=path
        )
    
    def write_density_cube(self, path: Path):
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
        return self._frequencies
