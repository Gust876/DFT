"""Página de listagem de arquivos da interface gráfica DFT."""
import streamlit as st
from html import escape
from pathlib import Path


def list_files(directory: Path) -> list[Path]:
    """Lista e ordena os arquivos de um diretório.

    Args:
        directory (Path): Diretório a ser listado.

    Returns:
        list[Path]: Lista ordenada de arquivos (excluindo subdiretórios).
    """
    if not directory.exists():
        return []
    return sorted([f for f in directory.iterdir() if f.is_file()])


def render_file_column(title: str, files: list[Path]) -> None:
    """Renderiza uma coluna de arquivos com pasta e contador.

    Args:
        title (str): Nome da pasta.
        files (list[Path]): Arquivos a exibir.
    """
    st.markdown(
        f'<div class="section-title">{escape(title)} '
        f'<span class="badge badge-pyscf">{len(files)}</span></div>',
        unsafe_allow_html=True,
    )
    if not files:
        st.markdown('<div class="file-empty">—</div>', unsafe_allow_html=True)
        return
    for f in files:
        st.markdown(
            f'<div class="file-item"><span class="file-name">{escape(f.name)}</span>'
            f'<span style="margin-left:auto">{f.stat().st_size/1024:.1f} KB</span></div>',
            unsafe_allow_html=True,
        )


def render_arquivos(xyz_dir: Path, opt_dir: Path, output_dir: Path) -> None:
    """Renderiza a aba de arquivos: entrada, otimizados e resultados.

    Args:
        xyz_dir (Path): Diretório de arquivos XYZ de entrada.
        opt_dir (Path): Diretório de geometrias otimizadas.
        output_dir (Path): Diretório de resultados finais (.zip).
    """
    col_a, col_b, col_c = st.columns(3, gap="large")
    with col_a:
        render_file_column("xyz_semi_opt", list_files(xyz_dir))
    with col_b:
        render_file_column("xyz_opt", list_files(opt_dir))
    with col_c:
        render_file_column("zip_dir", list_files(output_dir))

    if st.button("↻", key="refresh_files", help="Atualizar"):
        st.rerun()
