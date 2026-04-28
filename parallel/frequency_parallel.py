from strategies.run_strategy import SetFrequencyStrategy
from numpy import isrealobj
from pathlib import Path
import logging
from zipfile import ZipFile, ZIP_DEFLATED

LOGS_DIR = Path("./logs")
LOGS_DIR.mkdir(parents=True, exist_ok=True)

log_file = LOGS_DIR / "frequencies_process.log"

logging.basicConfig(
    filename=str(log_file),
    filemode='a',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)

def frequency(xyz_file: Path, OUTPUT_DIR: Path, **kwargs) -> None:
    file_name = xyz_file.stem
    data_file = OUTPUT_DIR / f"{file_name}.zip"

    if not data_file.exists():

        try:
            set_strategy = SetFrequencyStrategy(kwargs['engine'])

            result = set_strategy.frequency(
                xyz_file=xyz_file,
                basis=kwargs['basis'],
                xc=kwargs['xc']
            )
            frequencies = result.frequencies

            if isrealobj(frequencies):

                result.write_molden(OUTPUT_DIR / f"{file_name}.molden")
                result.write_density_cube(OUTPUT_DIR / f"{file_name}.cube")
                zip_files(OUTPUT_DIR, file_name)

                for ext in ("*.cube", "*.xyz", "*.dat"):
                    for file in Path(".").glob(ext):
                        file.unlink()
                
                logging.info(f"{xyz_file.name} concluída")

            else:
                logging.info(f"{xyz_file.name} possui frequências imaginárias")

        except Exception as error:
            logging.error(f"{xyz_file.name} falhou: {str(error)}")

def zip_files(
        output_dir: Path,
        file_name: str,
        remove_originals: bool = True
    ) -> None:
    
    molden_file = output_dir / f"{file_name}.molden"
    cube_file = output_dir / f"{file_name}.cube"
    zip_path = output_dir / f"{file_name}.zip"

    with ZipFile(zip_path, "w", compression=ZIP_DEFLATED) as myzip:

        if molden_file.exists():
            myzip.write(molden_file, arcname=molden_file.name)

        if cube_file.exists():
            myzip.write(cube_file, arcname=cube_file.name)

    if remove_originals:
        if molden_file.exists():
            molden_file.unlink()

        if cube_file.exists():
            cube_file.unlink()