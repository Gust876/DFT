from strategies.strategies import PySCFStrategy, Psi4Strategy
import os


def strategy_from_env() -> dict:

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
    