"""Página de listagem de arquivos da interface gráfica DFT."""
import streamlit as st
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
    """Renderiza uma coluna de listagem de arquivos.

    Args:
        title (str): Título da coluna.
        files (list[Path]): Lista de arquivos a exibir.
        icon (str): Ícone a exibir ao lado de cada arquivo.
    """
    st.markdown(f'<div class="section-title">{title}</div>', unsafe_allow_html=True)
    if files:
        for f in files:
            size = f.stat().st_size
            st.markdown(f"""
            <div class="file-item">
                <span class="file-name">{f.name}</span>
                <span style="margin-left:auto">{size/1024:.1f} KB</span>
            </div>
            """, unsafe_allow_html=True)
    else:
        st.markdown(
            '<div style="color:#8b949e; font-size:0.8rem; '
            'font-family:\'IBM Plex Mono\',monospace;">Nenhum arquivo</div>',
            unsafe_allow_html=True
        )


def render_arquivos(xyz_dir: Path, opt_dir: Path, output_dir: Path) -> None:
    """Renderiza a aba de listagem de arquivos do projeto.

    Exibe três colunas com os arquivos de entrada, geometrias otimizadas
    e resultados finais, com tamanho de cada arquivo e botão de atualização.

    Args:
        xyz_dir (Path): Diretório de arquivos XYZ de entrada.
        opt_dir (Path): Diretório de geometrias otimizadas.
        output_dir (Path): Diretório de resultados finais (.zip).
    """
    st.markdown("### Arquivos do projeto")

    col_a, col_b, col_c = st.columns(3)

    with col_a:
        render_file_column("Input (xyz_semi_opt)", list_files(xyz_dir))
    with col_b:
        render_file_column("Otimizados (xyz_opt)", list_files(opt_dir))
    with col_c:
        render_file_column("Resultados (zip_dir)", list_files(output_dir))

    if st.button("Atualizar lista de arquivos"):
        st.rerun()
