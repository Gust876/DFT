from strategies.strategies import PySCFStrategy, Psi4Strategy
import os


def strategy_from_env() -> dict:
    '''
    Constrói o dicionário de estratégia a partir de variáveis de ambiente.

    Utilizado pela interface gráfica (facade.py) para passar os parâmetros
    de execução ao workflow sem depender de input interativo do terminal.
    As variáveis de ambiente são definidas pelo processo pai (Streamlit)
    antes de invocar o subprocesso do workflow.

    Environment Variables:
        DFT_ENGINE (str): Engine de cálculo ('pyscf' ou 'psi4').
        DFT_BASIS (str): Conjunto de funções de base (ex: 'cc-pvdz').
        DFT_XC (str): Funcional de troca-correlação ou método (ex: 'm06-2x').

    Returns:
        dict: Dicionário com chaves:
            - 'engine' (EngineStrategy): Instância da estratégia da engine.
            - 'basis' (str): Conjunto de funções de base.
            - 'xc' (str): Funcional de troca-correlação ou método.

    Raises:
        ValueError: Se DFT_ENGINE não for 'pyscf' ou 'psi4', ou se
            DFT_BASIS ou DFT_XC estiverem vazios.
    '''

    engine_choices = os.environ.get("DFT_ENGINE", "").strip().lower()
    basis = os.environ.get("DFT_BASIS", "").strip()
    xc = os.environ.get("DFT_XC", "").strip()

    if engine_choices not in ["pyscf", "psi4"]:
        raise ValueError(
            f"Engine inválida: '{engine_choices}'. Defina DFT_ENGINE como'pyscf' ou 'psi4'."
        )
    if not basis:
        raise ValueError(
            "Basis set não definido. Defina a variável de ambiente DFT_BASIS."
        )
    if not xc:
        raise ValueError(
            "Funcional/método não definido. Defina a variável de ambiente DFT_XC."
        )

    if engine_choices == "pyscf":
        return {
            'engine': PySCFStrategy(),
            'basis': basis,
            'xc': xc
        }
    
    elif engine_choices == "psi4":
        return {
            'engine': Psi4Strategy(),
            'basis': basis,
            'xc': xc
        }


def prompt_strategy(workflow: str = "optimizer") -> dict:
    '''
    Coleta os parâmetros de execução interativamente via terminal.

    Utilizado pelos workflows de linha de comando (workflow_optimizer.py
    e workflow_frequencies.py). Solicita ao usuário a escolha da engine,
    o conjunto de funções de base e o funcional de troca-correlação.

    Args:
        workflow (str): Identificador do workflow sendo executado.
            Atualmente não utilizado na lógica interna, reservado
            para extensões futuras. Padrão: 'optimizer'.

    Returns:
        dict: Dicionário com chaves:
            - 'engine' (EngineStrategy): Instância da estratégia da engine.
            - 'basis' (str): Conjunto de funções de base.
            - 'xc' (str): Funcional de troca-correlação ou método.

    Raises:
        ValueError: Se a engine fornecida não for 'pyscf' ou 'psi4'.
    '''

    engine_choices = input(
        'Escolha uma engine (pyscf/psi4): '
    ).strip().lower()

    if engine_choices not in ["pyscf", "psi4"]:
        raise ValueError(
            f"Engine inválida: '{engine_choices}'. Escolha 'pyscf' ou 'psi4'."
        )
    basis = input('Defina um basis set: ').strip()

    if engine_choices == "pyscf":
        xc = input('Defina um xc: ').strip()
        return {
            'engine': PySCFStrategy(),
            'basis': basis,
            'xc': xc
        }
    
    elif engine_choices == "psi4":
        method = input('Defina um método: ').strip()
        return {
            'engine': Psi4Strategy(),
            'basis': basis,
            'xc': method
        }
    