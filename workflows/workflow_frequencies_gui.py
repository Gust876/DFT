from parallel.frequency_parallel import frequency
from utils.prompt import strategy_from_env
from joblib import Parallel, delayed
from pathlib import Path
import os
import logging

OPT_DIR    = Path("./xyz_opt")
OUTPUT_DIR = Path("./zip_dir")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

LOGS_DIR = Path("./logs")
LOGS_DIR.mkdir(parents=True, exist_ok=True)

log_file = LOGS_DIR / "parallel_frequencies.log"

logging.basicConfig(
    filename=str(log_file),
    filemode='a',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)

try:
    n_jobs   = int(os.environ.get("DFT_NJOBS", -1))
    strategy = strategy_from_env()

    Parallel(n_jobs=n_jobs)(
        delayed(frequency)(
            xyz_file=xyz_file,
            OUTPUT_DIR=OUTPUT_DIR,
            **strategy
        )
        for xyz_file in OPT_DIR.iterdir()
    )

except Exception as error:
    logging.error(f"erro ao paralelizar: {str(error)}")
