"""Página de visualização de logs da interface gráfica DFT."""
import streamlit as st
from html import escape
from pathlib import Path


def read_log(log_file: Path) -> str:
    """Lê o conteúdo de um arquivo de log.

    Args:
        log_file (Path): Caminho para o arquivo de log.

    Returns:
        str: Conteúdo do log ou "—" se vazio/inexistente.
    """
    if not log_file.exists():
        return "—"
    content = log_file.read_text()
    return content if content.strip() else "—"


def render_logs(logs_dir: Path) -> None:
    """Renderiza os logs de execução em sub-abas.

    Args:
        logs_dir (Path): Diretório onde os arquivos de log estão armazenados.
    """
    sources = [
        ("Otimização", "parallel_optimization.log"),
        ("Frequências", "parallel_frequencies.log"),
        ("Processo", "frequencies_process.log"),
        ("Confôrmeros", "conformers_process.log"),
    ]
    for tab, (_, filename) in zip(st.tabs([label for label, _ in sources]), sources):
        with tab:
            log = escape(read_log(logs_dir / filename))
            st.markdown(f'<div class="log-box">{log}</div>', unsafe_allow_html=True)

    if st.button("↻", key="refresh_logs", help="Atualizar"):
        st.rerun()