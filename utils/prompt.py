from strategies.strategies import PySCFStrategy, Psi4Strategy

def prompt_strategy() -> dict:

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
    