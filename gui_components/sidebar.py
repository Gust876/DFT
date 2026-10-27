"""Componente de sidebar da interface gráfica DFT."""
import streamlit as st
from pathlib import Path


def count_files(directory: Path, suffix: str = "") -> int:
    """Conta os arquivos em um diretório, opcionalmente filtrando por extensão.

    Args:
        directory (Path): Diretório a ser inspecionado.
        suffix (str): Extensão para filtrar (ex: '.xyz'). Vazio = todos os arquivos.

    Returns:
        int: Número de arquivos encontrados.
    """
    if not directory.exists():
        return 0
    return len([f for f in directory.iterdir()
                if f.is_file() and (not suffix or f.suffix == suffix)])


def render_sidebar(xyz_dir: Path, opt_dir: Path, output_dir: Path) -> dict:
    """Renderiza a sidebar de configuração e retorna os parâmetros escolhidos.

    Args:
        xyz_dir (Path): Diretório de arquivos XYZ de entrada.
        opt_dir (Path): Diretório de geometrias otimizadas.
        output_dir (Path): Diretório de resultados finais.

    Returns:
        dict: Configurações escolhidas ('engine', 'basis', 'xc', 'n_jobs',
            'engine_label') e contagens de arquivos ('n_input', 'n_opt', 'n_output').
    """
    with st.sidebar:
        st.markdown('<div class="main-title">DFT</div>', unsafe_allow_html=True)
        st.markdown('<div class="main-subtitle">Automation Pipeline</div>', unsafe_allow_html=True)

        engine = st.selectbox(
            "Engine", ["pyscf", "psi4"],
            format_func=lambda x: "PySCF" if x == "pyscf" else "Psi4"
        )
        basis = st.text_input("Basis set", value="cc-pvdz", help="ex: cc-pvdz, 6-311g**")
        xc_label = "Funcional" if engine == "pyscf" else "Método"
        xc = st.text_input(xc_label, value="m06-2x", help="ex: b3lyp, m06-2x, pbe0")
        n_jobs = st.slider("Núcleos", min_value=1, max_value=8, value=4,
                           help="Moléculas calculadas em paralelo")

    engine_label = "PySCF" if engine == "pyscf" else "Psi4"
    return {
        "engine": engine, "basis": basis, "xc": xc,
        "n_jobs": n_jobs, "engine_label": engine_label,
        "n_input": count_files(xyz_dir, ".xyz"),
        "n_opt": count_files(opt_dir, ".xyz"),
        "n_output": count_files(output_dir, ".zip"),
    }