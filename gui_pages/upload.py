"""Página de carregamento de moléculas (.xyz) a partir do computador do usuário."""
import streamlit as st
from pathlib import Path
 
from gui_pages.xyz_utils import (
    validate_xyz, normalize_xyz, sanitize_filename, save_xyz,
    choose_destination, flash, show_flash, page_title, chip,
)
 
 
def render_carregar(xyz_dir: Path, opt_dir: Path) -> None:
    """Renderiza a aba de carregamento de arquivos .xyz.
 
    Args:
        xyz_dir (Path): Diretório de geometrias iniciais (xyz_semi_opt).
        opt_dir (Path): Diretório de geometrias otimizadas (xyz_opt).
    """
    page_title("⬆", "Carregar moléculas")
    show_flash()
 
    version = st.session_state.setdefault("uploader_version", 0)
    files = st.file_uploader(
        "Arquivos .xyz", type=["xyz"], accept_multiple_files=True,
        key=f"uploader_{version}", label_visibility="collapsed",
    )
    if not files:
        return
 
    valid, chips = [], []
    for f in files:
        try:
            content = f.getvalue().decode("utf-8")
        except UnicodeDecodeError:
            chips.append(chip(f.name, False, "não é texto UTF-8"))
            continue
        ok, msg = validate_xyz(content)
        chips.append(chip(f.name, ok, msg))
        if ok:
            valid.append((sanitize_filename(f.name), normalize_xyz(content)))
 
    st.markdown(f'<div class="chip-list">{"".join(chips)}</div>', unsafe_allow_html=True)
 
    if not valid:
        return
 
    col_dest, col_over, col_btn = st.columns([4, 2, 3], vertical_alignment="center")
    with col_dest:
        dest_dir = choose_destination(xyz_dir, opt_dir, key="upload_dest")
    with col_over:
        overwrite = st.toggle("Sobrescrever", key="upload_overwrite")
    with col_btn:
        clicked = st.button(f"Salvar {len(valid)}", type="primary", key="btn_save_upload")
 
    if clicked:
        saved = [n for n, c in valid if save_xyz(c, n, dest_dir, overwrite)]
        skipped = [n for n, _ in valid if n not in saved]
 
        msg = f"{len(saved)} salvo(s) em {dest_dir.name}/"
        if skipped:
            msg += f" · {len(skipped)} já existia(m)"
        flash("success" if saved else "warning", msg)
 
        st.session_state["uploader_version"] += 1
        st.rerun()
 