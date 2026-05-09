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
    '''
    Executa o cálculo de frequências vibracionais de uma molécula de forma isolada.

    Função projetada para execução paralela via Joblib. Cada chamada processa
    uma geometria otimizada, verificando se o resultado já existe antes de iniciar
    o cálculo (retomada automática). Apenas moléculas com todas as frequências
    reais (estado fundamental confirmado) têm seus arquivos de saída gerados.

    Args:
        xyz_file (Path): Caminho para o arquivo XYZ com a geometria otimizada.
        OUTPUT_DIR (Path): Diretório de saída para os arquivos comprimidos (.zip).
        **kwargs:
            engine (EngineStrategy): Estratégia de engine configurada
                (PySCFStrategy ou Psi4Strategy).
            basis (str): Conjunto de funções de base (ex: 'cc-pvdz').
            xc (str): Funcional de troca-correlação ou método (ex: 'm06-2x').

    Returns:
        None

    Side effects:
        - Gera arquivo .molden com orbitais moleculares em OUTPUT_DIR.
        - Gera arquivo .cube com densidade eletrônica em OUTPUT_DIR.
        - Comprime ambos em um arquivo .zip e remove os originais.
        - Remove arquivos temporários (.cube, .xyz, .dat) do diretório de trabalho.
        - Registra resultado (sucesso, aviso de frequência imaginária ou erro) no log.

    Note:
        Moléculas com frequências imaginárias (complexas) indicam que a geometria
        não corresponde a um mínimo de energia e, portanto, não são processadas.
    '''
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
    '''
    Comprime os arquivos de saída (.molden e .cube) em um arquivo .zip.

    Args:
        output_dir (Path): Diretório onde os arquivos de saída estão localizados
            e onde o .zip será criado.
        file_name (str): Nome base dos arquivos (sem extensão).
        remove_originals (bool): Se True, remove os arquivos originais após a
            compressão. Padrão: True.

    Returns:
        None

    Side effects:
        - Cria arquivo .zip contendo os arquivos .molden e .cube em output_dir.
        - Remove os arquivos originais se remove_originals for True. 
    '''
    
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