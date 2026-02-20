from pyscf.tools import molden, cubegen
from zipfile import ZIP_DEFLATED, ZipFile
from pathlib import Path

def file_molden(
            mol,
            file_name,
            DATA_DIR: Path,
            coeff_mo,
            energy_mo,
            occ_mo
            ):
    
    path_dir = DATA_DIR / f"{file_name}.molden"
    
    with open(path_dir, "w") as f1:
        molden.header(mol, f1)
        molden.orbital_coeff(
                mol,
                f1,
                coeff_mo,
                ene=energy_mo,
                occ=occ_mo
                )
        
    return None

def file_cubegen(mol, dm, file_name, DATA_DIR):
    path_dir = DATA_DIR / f"{file_name}_den.cube"
    cubegen.density(mol, str(path_dir), dm)

    return None

def zip_files(DATA_DIR: Path, file_name):
    molden_file = DATA_DIR / f"{file_name}.molden"
    cube_file = DATA_DIR / f"{file_name}_den.cube"
    zip_path = DATA_DIR / f"{file_name}.zip"

    with ZipFile(zip_path, "w", compression=ZIP_DEFLATED) as myzip:

        if molden_file.exists():
            myzip.write(molden_file, arcname=molden_file.name)

        if cube_file.exists():
            myzip.write(cube_file, arcname=cube_file.name)

    molden_file.unlink()
    cube_file.unlink()

    return None
