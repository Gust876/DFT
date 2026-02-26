from parallel.optimizer_parallel import geometry_optimizer
from joblib import Parallel, delayed
from pathlib import Path
import logging

XYZ_DIR = Path("./xyz_semi_opt")
OPT_DIR = Path("./xyz_opt")
OPT_DIR.mkdir(exist_ok=True)

LOGS_DIR = Path("./logs")
LOGS_DIR.mkdir(parents=True, exist_ok=True)

log_file = LOGS_DIR / "parallel_optimization.log"

logging.basicConfig(
    filename=str(log_file),
    filemode='a',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)


try:
    Parallel(n_jobs=8)(
        delayed(geometry_optimizer)
        (
            xyz_file=xyz_file,
            OPT_DIR=OPT_DIR
        )
        for xyz_file in XYZ_DIR.iterdir()
    )

except Exception as error:
    logging.error(f"erro ao paralelizar: {str(error)}")
