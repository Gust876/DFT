from pathlib import Path
import logging
from strategies.run_strategy import SetOptimizerStrategy

LOGS_DIR = Path("./logs")
LOGS_DIR.mkdir(exist_ok=True)

log_file = LOGS_DIR / "parallel_optimization.log"

logging.basicConfig(
    filename=str(log_file),
    filemode='a',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)

def geometry_optimizer(
        xyz_file: Path,
        OPT_DIR: Path,
        **kwargs
    ) -> None:
    '''
    Executa a otimização geométrica de uma molécula de forma isolada.

    Função projetada para execução paralela via Joblib. Cada chamada
    processa um único arquivo XYZ, verificando se o resultado já existe
    antes de iniciar o cálculo (retomada automática de execuções interrompidas).
    O resultado convergido é salvo como arquivo XYZ na pasta de saída.

    Args:
        xyz_file (Path): Caminho para o arquivo XYZ com a geometria inicial
            semi-otimizada.
        OPT_DIR (Path): Diretório de saída para as geometrias otimizadas.
        **kwargs:
            engine (EngineStrategy): Estratégia de engine configurada
                (PySCFStrategy ou Psi4Strategy).
            basis (str): Conjunto de funções de base (ex: 'cc-pvdz').
            xc (str): Funcional de troca-correlação ou método (ex: 'm06-2x').

    Returns:
        None

    Side effects:
        - Cria arquivo XYZ com prefixo 'opt_' em OPT_DIR se a otimização convergir.
        - Registra resultado (convergência, aviso ou erro) no arquivo de log.
    '''

    try:
        opt_xyz_file = OPT_DIR / f"opt_{xyz_file.name}"
        if not opt_xyz_file.exists():

            set_strategy = SetOptimizerStrategy(kwargs['engine'])

            dict_optimize = set_strategy.optimizer(
                xyz_file=xyz_file,
                basis=kwargs['basis'],
                xc=kwargs['xc']
            )

            converged = dict_optimize['converged']
            if converged:

                opt_xyz_file.write_text(dict_optimize['xyz_data'])
                logging.info(f"{xyz_file.name} convergida")
                
            else:
                logging.warning(f"{xyz_file.name} - error: {dict_optimize['error']}")

    except Exception as error:
        logging.error(f"{xyz_file.name} falhou: {str(error)}")
        