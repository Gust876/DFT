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

    Exibe os controles de seleção de engine, basis set, funcional e
    número de núcleos, além dos contadores de arquivos por diretório.

    Args:
        xyz_dir (Path): Diretório de arquivos XYZ de entrada.
        opt_dir (Path): Diretório de geometrias otimizadas.
        output_dir (Path): Diretório de resultados finais.

    Returns:
        dict: Dicionário com as configurações escolhidas:
            - 'engine' (str): 'pyscf' ou 'psi4'.
            - 'basis' (str): Conjunto de funções de base.
            - 'xc' (str): Funcional de troca-correlação ou método.
            - 'n_jobs' (int): Número de núcleos para paralelização.
            - 'engine_label' (str): Nome formatado da engine para exibição.
            - 'n_input' (int): Número de arquivos de entrada.
            - 'n_opt' (int): Número de geometrias otimizadas.
            - 'n_output' (int): Número de resultados gerados.
    """
    with st.sidebar:
        st.markdown('<div class="main-title">⚛ DFT</div>', unsafe_allow_html=True)
        st.markdown('<div class="main-subtitle">Automation Pipeline</div>', unsafe_allow_html=True)
        st.markdown("---")

        st.markdown('<div class="section-title">Configuração</div>', unsafe_allow_html=True)

        engine = st.selectbox(
            "Engine", ["pyscf", "psi4"],
            format_func=lambda x: "PySCF" if x == "pyscf" else "Psi4"
        )
        basis = st.text_input("Basis Set", value="cc-pvdz", placeholder="ex: cc-pvdz, 6-311g**")
        xc_label = "Funcional XC" if engine == "pyscf" else "Método"
        xc = st.text_input(xc_label, value="m06-2x", placeholder="ex: b3lyp, m06-2x, pbe0")
        n_jobs = st.slider("Núcleos paralelos", min_value=1, max_value=8, value=4)

        st.markdown("---")
        st.markdown('<div class="section-title">Status</div>', unsafe_allow_html=True)

        n_input  = count_files(xyz_dir,    ".xyz")
        n_opt    = count_files(opt_dir,    ".xyz")
        n_output = count_files(output_dir, ".zip")

        col1, col2, col3 = st.columns(3)
        for col, val, label in zip([col1, col2, col3],
                                   [n_input, n_opt, n_output],
                                   ["Input", "Otim.", "Saída"]):
            with col:
                st.markdown(
                    f'<div class="metric-box">'
                    f'<div class="metric-value">{val}</div>'
                    f'<div class="metric-label">{label}</div>'
                    f'</div>',
                    unsafe_allow_html=True
                )

        engine_badge = "badge-pyscf" if engine == "pyscf" else "badge-psi4"
        engine_label = "PySCF" if engine == "pyscf" else "Psi4"
        st.markdown(f"""
        <br>
        <div style="font-size:0.78rem; color:#8b949e; font-family:'IBM Plex Mono',monospace;">
            Engine: <span class="badge {engine_badge}">{engine_label}</span><br><br>
            Basis: <span style="color:#e6edf3">{basis or '—'}</span><br><br>
            XC/Método: <span style="color:#e6edf3">{xc or '—'}</span>
        </div>
        """, unsafe_allow_html=True)

    return {
        "engine": engine, "basis": basis, "xc": xc,
        "n_jobs": n_jobs, "engine_label": engine_label,
        "n_input": n_input, "n_opt": n_opt, "n_output": n_output
    }
