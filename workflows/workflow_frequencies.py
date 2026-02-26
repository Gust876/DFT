from parallel.frequency_parallel import frequency
from joblib import Parallel, delayed
from pathlib import Path
import logging

OPT_DIR = Path("./xyz_opt")
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
    Parallel(n_jobs=-1)(
        delayed(frequency)
        (
            xyz_file=xyz_file,
            OUTPUT_DIR=OUTPUT_DIR
        )
        for xyz_file in OPT_DIR.iterdir()
    )

except Exception as error:
    logging.error(f"erro ao paralelizar: {str(error)}")
