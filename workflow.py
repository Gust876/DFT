from engine_lib.pyscf_lib.opt_pyscf import OptimizerPySCF
from engine_lib.pyscf_lib.engine_pyscf import EnginePySCF
from engine_lib.pyscf_lib.utils_pyscf import create_molecule
from pathlib import Path
import logging

XYZ_DIR = Path("./xyz_semi_opt")
OPT_DIR = Path("./xyz_opt")
OPT_DIR.mkdir(exist_ok=True)

logging.basicConfig(
    filename='optimization_process.log',
    filemode='a',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)

for xyz_file in XYZ_DIR.iterdir():

    if not xyz_file.is_file() or xyz_file.suffix != ".xyz":
        continue

    try:
        opt_xyz_file = OPT_DIR / f"opt_{xyz_file.name}"
        if not opt_xyz_file.exists():

            mol = create_molecule(xyz_file=str(xyz_file), basis="cc-pvdz")

            engine_pyscf =  EnginePySCF(mol=mol, xc="b3lyp")
            opt_geom = OptimizerPySCF(engine=engine_pyscf)
            dict_optimize = opt_geom.opt_geometry(maxsteps=200, verbose=0)

            converged = dict_optimize["converged"]
            if converged:

                opt_xyz_file.write_text(dict_optimize["xyz_data"])
                logging.info(f"{xyz_file.name} convergida")
                
            else:
                logging.warning(f"{xyz_file.name} - error: {dict_optimize['error']}")

    except Exception as error:
        logging.error(f"{xyz_file.name} falhou: {str(error)}")
        continue
