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


def render_execucao(config: dict, base_dir: Path) -> None:
    """Renderiza a aba de execução dos workflows.

    Exibe os cards de otimização geométrica e frequências vibracionais
    com botões de execução e feedback de resultado.

    Args:
        config (dict): Configurações retornadas pela sidebar
            (engine, basis, xc, n_jobs, engine_label, n_input, n_opt).
        base_dir (Path): Diretório raiz do projeto.
    """
    st.markdown("### Workflows")
    st.markdown("Configure os parâmetros na barra lateral e execute o workflow desejado.")
    st.markdown("")

    col_opt, col_freq = st.columns(2)

    with col_opt:
        st.markdown("""
        <div class="section-card">
            <div class="section-title">Otimização Geométrica</div>
            <div style="font-size:0.82rem; color:#8b949e; margin-bottom:1rem; line-height:1.6;">
                Lê os arquivos <span style="color:#e6edf3; font-family:'IBM Plex Mono',monospace;">.xyz</span>
                de <code>xyz_semi_opt/</code>, realiza a otimização geométrica
                e salva os resultados em <code>xyz_opt/</code>.
            </div>
        </div>
        """, unsafe_allow_html=True)

        if st.button("Executar Otimização", key="btn_opt"):
            if not config["basis"].strip() or not config["xc"].strip():
                st.error("Preencha o Basis Set e o Funcional/Método antes de executar.")
            elif config["n_input"] == 0:
                st.warning("Nenhum arquivo .xyz encontrado em xyz_semi_opt/")
            else:
                with st.spinner(f"Otimizando {config['n_input']} molécula(s) com {config['engine_label']}..."):
                    ok, output = run_workflow(
                        "optimizer", config["engine"], config["basis"],
                        config["xc"], config["n_jobs"], base_dir
                    )
                    if ok:
                        st.success("Otimização concluída. Veja os resultados na aba Arquivos.")
                    else:
                        st.error(f"Erro:\n{output}")

    with col_freq:
        st.markdown("""
        <div class="section-card">
            <div class="section-title">Frequências Vibracionais</div>
            <div style="font-size:0.82rem; color:#8b949e; margin-bottom:1rem; line-height:1.6;">
                Lê as geometrias otimizadas de <code>xyz_opt/</code>, calcula as
                frequências vibracionais e gera os arquivos
                <span style="color:#e6edf3; font-family:'IBM Plex Mono',monospace;">.cube</span> e
                <span style="color:#e6edf3; font-family:'IBM Plex Mono',monospace;">.molden</span>
                em <code>zip_dir/</code>.
            </div>
        </div>
        """, unsafe_allow_html=True)

        if st.button("Executar Frequências", key="btn_freq"):
            if not config["basis"].strip() or not config["xc"].strip():
                st.error("Preencha o Basis Set e o Funcional/Método antes de executar.")
            elif config["n_opt"] == 0:
                st.warning("Nenhuma geometria otimizada encontrada em xyz_opt/")
            else:
                with st.spinner(f"Calculando frequências de {config['n_opt']} molécula(s)..."):
                    ok, output = run_workflow(
                        "frequencies", config["engine"], config["basis"],
                        config["xc"], config["n_jobs"], base_dir
                    )
                    if ok:
                        st.success("Frequências concluídas. Veja os resultados na aba Arquivos.")
                    else:
                        st.error(f"Erro:\n{output}")

    st.markdown("---")
    st.markdown("### Como usar")
    st.markdown("""
    1. **Configure** a engine, basis set e funcional na barra lateral
    2. **Coloque** os arquivos `.xyz` semi-otimizados em `xyz_semi_opt/`
    3. **Execute** a Otimização Geométrica primeiro
    4. **Execute** o cálculo de Frequências após a otimização concluir
    5. **Confira** os resultados em `zip_dir/` na aba Arquivos
    """)
