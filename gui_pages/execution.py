"""Página de execução dos workflows da interface gráfica DFT."""
import streamlit as st
import subprocess
import sys
from pathlib import Path


def run_workflow(workflow: str, engine: str, basis: str, xc: str,
                 n_jobs: int, base_dir: Path) -> tuple[bool, str]:
    """Executa o workflow como subprocesso com os parâmetros via variáveis de ambiente.

    Args:
        workflow (str): Identificador do workflow ('optimizer' ou 'frequencies').
        engine (str): Engine de cálculo ('pyscf' ou 'psi4').
        basis (str): Conjunto de funções de base.
        xc (str): Funcional de troca-correlação ou método.
        n_jobs (int): Número de núcleos para paralelização.
        base_dir (Path): Diretório raiz do projeto.

    Returns:
        tuple[bool, str]: (sucesso, mensagem de saída ou erro).
    """
    import os
    env = os.environ.copy()
    env["PYTHONPATH"] = str(base_dir)
    env["DFT_ENGINE"] = engine
    env["DFT_BASIS"]  = basis
    env["DFT_XC"]     = xc
    env["DFT_NJOBS"]  = str(n_jobs)

    script = base_dir / "workflows" / f"workflow_{workflow}_gui.py"

    try:
        result = subprocess.run(
            [sys.executable, str(script)],
            capture_output=True, text=True, env=env
        )
        return (True, result.stdout) if result.returncode == 0 else (False, result.stderr)
    except Exception as e:
        return False, str(e)


def _node(num: int, label: str) -> None:
    """Etapa do pipeline: contagem de arquivos e rótulo."""
    st.markdown(
        f'<div class="flow-node">'
        f'<div class="flow-num">{num}</div><div class="flow-label">{label}</div></div>',
        unsafe_allow_html=True,
    )


def _run(workflow: str, count: int, config: dict, base_dir: Path) -> None:
    """Valida, executa o workflow e notifica o resultado."""
    if not config["basis"].strip() or not config["xc"].strip():
        st.toast("Preencha basis set e funcional na barra lateral")
    elif count == 0:
        st.toast("Nenhuma molécula nesta etapa — use Carregar ou Criar")
    else:
        with st.spinner(f"{count} molécula(s) · {config['engine_label']}"):
            ok, output = run_workflow(
                workflow, config["engine"], config["basis"],
                config["xc"], config["n_jobs"], base_dir
            )
        if ok:
            st.toast("Concluído")
            st.rerun()
        else:
            st.error(output)


def render_execucao(config: dict, base_dir: Path) -> None:
    """Renderiza o pipeline: etapas com contadores e botões entre elas.

    Args:
        config (dict): Configurações retornadas pela sidebar.
        base_dir (Path): Diretório raiz do projeto.
    """
    cols = st.columns(5, vertical_alignment="center", gap="small")
    with cols[0]:
        _node(config["n_input"], "Iniciais")
    with cols[1]:
        opt = st.button("Otimizar →", key="btn_opt", type="primary")
    with cols[2]:
        _node(config["n_opt"], "Otimizadas")
    with cols[3]:
        freq = st.button("Frequências →", key="btn_freq", type="primary")
    with cols[4]:
        _node(config["n_output"], "Resultados")

    if opt:
        _run("optimizer", config["n_input"], config, base_dir)
    if freq:
        _run("frequencies", config["n_opt"], config, base_dir)
