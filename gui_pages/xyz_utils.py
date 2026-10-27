"""Utilitários compartilhados para carregar/criar arquivos .xyz na GUI."""
import re
from html import escape
import streamlit as st
from pathlib import Path

_SYMBOL = re.compile(r"^[A-Za-z]{1,2}$")


def validate_xyz(content: str) -> tuple[bool, str]:
    """Valida o formato XYZ padrão (nº de átomos, comentário, átomos).

    Args:
        content (str): Conteúdo completo do arquivo.

    Returns:
        tuple[bool, str]: (válido, mensagem de erro ou resumo).
    """
    lines = content.strip().splitlines()
    if len(lines) < 3:
        return False, "arquivo muito curto"
    try:
        n_atoms = int(lines[0].strip())
    except ValueError:
        return False, "1ª linha deve ser o nº de átomos"

    atoms = [ln for ln in lines[2:] if ln.strip()]
    if n_atoms <= 0 or len(atoms) != n_atoms:
        return False, f"cabeçalho: {n_atoms} átomos, encontrados: {len(atoms)}"

    for i, ln in enumerate(atoms, start=3):
        parts = ln.split()
        if len(parts) < 4 or not _SYMBOL.match(parts[0]):
            return False, f"linha {i}: esperado 'Elemento x y z'"
        try:
            [float(v) for v in parts[1:4]]
        except ValueError:
            return False, f"linha {i}: coordenadas inválidas"
    return True, f"{n_atoms} átomos"


def normalize_xyz(content: str) -> str:
    """Remove espaços/linhas extras e garante quebra de linha final."""
    return content.strip().replace("\r\n", "\n") + "\n"


def sanitize_filename(name: str) -> str:
    """Gera um nome de arquivo seguro, sempre com extensão .xyz."""
    stem = Path(name).stem
    stem = re.sub(r"[^\w\-.]", "_", stem.strip()) or "molecula"
    return f"{stem}.xyz"


def save_xyz(content: str, filename: str, dest_dir: Path, overwrite: bool = False) -> bool:
    """Salva o conteúdo em dest_dir. Retorna False se existir e overwrite=False."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    target = dest_dir / filename
    if target.exists() and not overwrite:
        return False
    target.write_text(content, encoding="utf-8")
    return True


def page_title(icon: str, title: str) -> None:
    """Título compacto de página."""
    st.markdown(f'<div class="page-title"><span>{icon}</span>{escape(title)}</div>',
                unsafe_allow_html=True)


def chip(name: str, ok: bool, meta: str) -> str:
    """HTML de um item de arquivo com status (✓/✕) e metadado."""
    cls, icon = ("ok", "✓") if ok else ("err", "✕")
    return (f'<div class="chip {cls}"><span class="chip-icon">{icon}</span>'
            f'<span class="chip-name">{escape(name)}</span>'
            f'<span class="chip-meta">{escape(meta)}</span></div>')


def choose_destination(xyz_dir: Path, opt_dir: Path, key: str) -> Path:
    """Seletor compacto de destino (explicação apenas no tooltip)."""
    options = {"Para otimizar": xyz_dir, "Já otimizada": opt_dir}
    label = st.radio(
        "Destino", list(options), key=key, horizontal=True,
        label_visibility="collapsed",
        help="Para otimizar → xyz_semi_opt/\n\nJá otimizada → xyz_opt/ (segue direto para as frequências)",
    )
    return options[label]


def flash(kind: str, message: str) -> None:
    """Guarda uma notificação para ser exibida após o st.rerun()."""
    st.session_state["_flash"] = (kind, message)


def show_flash() -> None:
    """Exibe (e consome) a notificação guardada por flash() como toast."""
    item = st.session_state.pop("_flash", None)
    if item:
        st.toast(item[1])