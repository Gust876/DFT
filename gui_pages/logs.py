"""Página de visualização de logs da interface gráfica DFT."""
import streamlit as st
from pathlib import Path


def read_log(log_file: Path) -> str:
    """Lê o conteúdo de um arquivo de log.

    Args:
        log_file (Path): Caminho para o arquivo de log.

    Returns:
        str: Conteúdo do log ou mensagem padrão se vazio/inexistente.
    """
    if not log_file.exists():
        return "Nenhum log disponível ainda."
    content = log_file.read_text()
    return content if content.strip() else "Nenhum log disponível ainda."


def render_logs(logs_dir: Path) -> None:
    """Renderiza a aba de visualização de logs de execução.

    Exibe os logs de otimização, frequências e processamento
    lado a lado com botão de atualização.

    Args:
        logs_dir (Path): Diretório onde os arquivos de log estão armazenados.
    """
    st.markdown("### Logs de execução")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown('<div class="section-title">Otimização</div>', unsafe_allow_html=True)
        log = read_log(logs_dir / "parallel_optimization.log")
        st.markdown(f'<div class="log-box">{log}</div>', unsafe_allow_html=True)

    with col2:
        st.markdown('<div class="section-title">Frequências</div>', unsafe_allow_html=True)
        log = read_log(logs_dir / "parallel_frequencies.log")
        st.markdown(f'<div class="log-box">{log}</div>', unsafe_allow_html=True)

    st.markdown("")
    st.markdown('<div class="section-title">Processamento de frequências</div>', unsafe_allow_html=True)
    log = read_log(logs_dir / "frequencies_process.log")
    st.markdown(f'<div class="log-box">{log}</div>', unsafe_allow_html=True)

    if st.button("Atualizar logs"):
        st.rerun()
