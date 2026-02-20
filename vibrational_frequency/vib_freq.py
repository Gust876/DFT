from vibrational_frequency.generate_files import file_molden, file_cubegen, zip_files
from engine_lib.pyscf_lib.utils_pyscf import create_molecule
from engine_lib.pyscf_lib.engine_pyscf import EnginePySCF
from pyscf.hessian import thermo
from numpy import isrealobj
from pathlib import Path

OPT_DIR = Path("./xyz_opt")
DATA_DIR = Path("./data_dir")
DATA_DIR.mkdir(parents=True, exist_ok=True)

for xyz_file in OPT_DIR.iterdir():
    if not xyz_file.is_file() or xyz_file.suffix != ".xyz":
        continue

    file_name = xyz_file.stem
    data_file = DATA_DIR / f"{file_name}.zip"

    if not data_file.exists():

        mol = create_molecule(xyz_file=str(xyz_file), basis="cc-pvdz")
        engine = EnginePySCF(mol=mol, xc="b3lyp")
        mf = engine.dft_method()
        mf.run()

        hessian = mf.Hessian().kernel()
        freq_info = thermo.harmonic_analysis(mf.mol, hessian)


        if isrealobj(freq_info["freq_wavenumber"]):

            file_molden(
                mol=mol,
                file_name=file_name,
                DATA_DIR=DATA_DIR,
                coeff_mo=mf.mo_coeff,
                energy_mo=mf.mo_energy,
                occ_mo=mf.mo_occ
            )

            file_cubegen(
                mol=mol,
                dm=mf.make_rdm1(),
                file_name=file_name,
                DATA_DIR=DATA_DIR
            )

            del engine, hessian, mol, mf

            zip_files(DATA_DIR, file_name)


# print(freq_info["freq_wavenumber"])
# print(freq_info["norm_mode"])
# print(freq_info["norm_mode"][0])
# print()
# print(mol.atom_coords())
    