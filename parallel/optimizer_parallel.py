from engines.pyscf.optimizer import FactoryOptimizerPySCF
from engines.utils import factory_mol, factory_optimizer
from engines.pyscf.mol import FactoryMolPySCF
from pathlib import Path
import logging

LOGS_DIR = Path("./logs")
LOGS_DIR.mkdir(exist_ok=True)

log_file = LOGS_DIR / "frequencies_process.log"

logging.basicConfig(
    filename=str(log_file),
    filemode='a',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)

def geometry_optimizer(xyz_file: Path, OPT_DIR: Path) -> None:

    try:
        opt_xyz_file = OPT_DIR / f"opt_{xyz_file.name}"
        if not opt_xyz_file.exists():

            construct_mol = factory_mol(FactoryMolPySCF())
            mol = construct_mol.create_mol(
                xyz_file=xyz_file.read_text(),
                basis="cc-pvdz"
            )

            optimizer = factory_optimizer(
                FactoryOptimizerPySCF(
                    mol=mol,
                    xc="b3lyp"
                )
            )
            dict_optimize = optimizer.opt_geometry(maxsteps=200, verbose=0)

            converged = dict_optimize["converged"]
            if converged:

                opt_xyz_file.write_text(dict_optimize["xyz_data"])
                logging.info(f"{xyz_file.name} convergida")
                
            else:
                logging.warning(f"{xyz_file.name} - error: {dict_optimize["error"]}")

    except Exception as error:
        logging.error(f"{xyz_file.name} falhou: {str(error)}")
        