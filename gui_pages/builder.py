"""Página de criação de estruturas iniciais (.xyz) pela interface gráfica."""
import streamlit as st
import pandas as pd
from html import escape
from pathlib import Path

from gui_pages.xyz_utils import (
    validate_xyz, normalize_xyz, sanitize_filename, save_xyz,
    choose_destination, flash, show_flash, page_title,
)

DEFAULT_ATOMS = pd.DataFrame(
    [
        {"Elemento": "O", "x": 0.000000, "y": 0.000000, "z": 0.117790},
        {"Elemento": "H", "x": 0.000000, "y": 0.755453, "z": -0.471161},
        {"Elemento": "H", "x": 0.000000, "y": -0.755453, "z": -0.471161},
    ]
)


def table_to_xyz(df: pd.DataFrame, comment: str) -> str:
    """Converte a tabela de átomos para o formato XYZ."""
    df = df.dropna(subset=["Elemento", "x", "y", "z"])
    df = df[df["Elemento"].astype(str).str.strip() != ""]
    rows = [
        f"{str(r.Elemento).strip().capitalize():<3s} {r.x:>14.8f} {r.y:>14.8f} {r.z:>14.8f}"
        for r in df.itertuples()
    ]
    return f"{len(rows)}\n{comment}\n" + "\n".join(rows) + "\n"


def smiles_to_xyz(smiles: str, comment: str) -> str:
    """Gera geometria 3D a partir de SMILES (requer RDKit)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("SMILES inválido.")
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        raise ValueError("Não foi possível gerar a geometria 3D.")
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    except Exception:
        AllChem.UFFOptimizeMolecule(mol, maxIters=2000)

    atoms = Chem.MolToXYZBlock(mol).strip().splitlines()[2:]
    return f"{len(atoms)}\n{comment}\n" + "\n".join(atoms) + "\n"


def _preview(content: str, name: str, xyz_dir: Path, opt_dir: Path, key: str) -> None:
    """Card de prévia com validação, destino e botão de salvar."""
    ok, msg = validate_xyz(content)
    badge = "badge-ok" if ok else "badge-psi4"
    st.markdown(
        f'<div class="section-title">Prévia &nbsp;<span class="badge {badge}">{escape(msg)}</span></div>',
        unsafe_allow_html=True,
    )
    st.code(content, language="text")
    if not ok:
        return

    dest_dir = choose_destination(xyz_dir, opt_dir, key=f"{key}_dest")
    overwrite = st.toggle("Sobrescrever", key=f"{key}_overwrite")
    filename = sanitize_filename(name)

    if st.button(f"Salvar {filename}", type="primary", key=f"{key}_save"):
        if save_xyz(normalize_xyz(content), filename, dest_dir, overwrite):
            flash("success", f"{filename} salvo em {dest_dir.name}/")
            st.rerun()
        else:
            st.toast(f"{filename} já existe — ative “Sobrescrever”")


def render_criar(xyz_dir: Path, opt_dir: Path) -> None:
    """Renderiza a aba de criação de estruturas iniciais.

    Args:
        xyz_dir (Path): Diretório de geometrias iniciais (xyz_semi_opt).
        opt_dir (Path): Diretório de geometrias otimizadas (xyz_opt).
    """
    show_flash()

    col_name, col_mode = st.columns([3, 2])
    with col_name:
        name = st.text_input("Nome", value="molecula", key="builder_name")
    with col_mode:
        mode = st.radio("Método", ["Tabela", "SMILES"], horizontal=True, key="builder_mode")
    comment = name.strip() or "molecula"

    st.markdown("")
    col_left, col_right = st.columns([3, 2], gap="large")
    content = None

    with col_left:
        if mode == "Tabela":
            df = st.data_editor(
                DEFAULT_ATOMS, num_rows="dynamic", key="builder_table",
                column_config={
                    "Elemento": st.column_config.TextColumn("Elemento"),
                    "x": st.column_config.NumberColumn("x (Å)", format="%.6f"),
                    "y": st.column_config.NumberColumn("y (Å)", format="%.6f"),
                    "z": st.column_config.NumberColumn("z (Å)", format="%.6f"),
                },
            )
            content = table_to_xyz(df, comment)
        else:
            try:
                import rdkit  # noqa: F401
            except ImportError:
                st.info("Requer RDKit: `conda install -c conda-forge rdkit`")
            else:
                smiles = st.text_input("SMILES", value="CCO", placeholder="ex: CCO, c1ccccc1", key="builder_smiles")
                if st.button("Gerar 3D", key="btn_gen_smiles"):
                    try:
                        st.session_state["smiles_xyz"] = smiles_to_xyz(smiles.strip(), comment)
                    except Exception as e:
                        st.session_state.pop("smiles_xyz", None)
                        st.toast(str(e))
                content = st.session_state.get("smiles_xyz")

    with col_right:
        if content:
            _preview(content, name, xyz_dir, opt_dir, key="table" if mode == "Tabela" else "smiles")
 