"""Página de busca conformacional (RDKit) da interface gráfica DFT."""
import sys
import streamlit as st
import pandas as pd
from pathlib import Path

from gui_pages.xyz_utils import chip, flash, show_flash

MAX_CHIPS = 8


def _results_table(results: list[dict]) -> pd.DataFrame:
    """Converte o retorno de run_conformer_search em tabela de exibição."""
    rows = []
    for r in results:
        if r["status"] == "ok":
            for c in r["conformers"]:
                rows.append({
                    "Molécula": r["id"], "Arquivo": c["file"],
                    "E (kcal/mol)": round(c["energy"], 3),
                    "ΔE (kcal/mol)": round(c["rel_energy"], 3),
                    "Campo de força": r["force_field"], "Status": "ok",
                })
        else:
            rows.append({
                "Molécula": r["id"], "Arquivo": "—", "E (kcal/mol)": None,
                "ΔE (kcal/mol)": None, "Campo de força": "—", "Status": r["message"],
            })
    return pd.DataFrame(rows)


def render_conformeros(xyz_dir: Path) -> None:
    """Renderiza a aba de busca conformacional.

    Os confôrmeros de menor energia são gravados em xyz_semi_opt/,
    prontos para a otimização geométrica DFT.

    Args:
        xyz_dir (Path): Diretório de geometrias iniciais (xyz_semi_opt).
    """
    show_flash()

    try:
        from rdkit import Chem  # noqa: F401
    except ImportError as e:
        st.info(f"RDKit indisponível em `{sys.executable}`: {e}")
        return

    try:
        from rdkit_lib.conformer_pipeline import (
            parse_entries, parse_json_entries, describe_smiles, run_conformer_search,
        )
    except ImportError as e:
        st.error(f"Falha ao importar rdkit_lib: {e}")
        return

    mode = st.radio("Entrada", ["Texto", "JSON"], horizontal=True,
                    label_visibility="collapsed", key="conf_mode")

    entries, error_chip = [], None
    if mode == "Texto":
        text = st.text_area("SMILES", height=150, label_visibility="collapsed",
                            placeholder="nome SMILES  (um por linha)", key="conf_text")
        entries = parse_entries(text)
    else:
        upload = st.file_uploader("JSON", type=["json"], label_visibility="collapsed",
                                  key="conf_json")
        if upload:
            try:
                entries = parse_json_entries(upload.getvalue())
            except ValueError as e:
                error_chip = chip(upload.name, False, str(e))

    valid, chips = [], []
    for entry in entries:
        ok, meta = describe_smiles(entry["canonical_smiles"])
        if ok:
            valid.append(entry)
        if len(chips) < MAX_CHIPS:
            chips.append(chip(entry["id"], ok, meta))
    if error_chip:
        chips.append(error_chip)
    if chips:
        extra = len(entries) - MAX_CHIPS
        more = (f'<div class="chip-meta">+{extra}</div>' if extra > 0 else "")
        st.markdown(f'<div class="chip-list">{"".join(chips)}{more}</div>',
                    unsafe_allow_html=True)

    with st.expander("Parâmetros"):
        c1, c2, c3 = st.columns(3)
        with c1:
            num_confs = st.number_input("Gerar", min_value=10, max_value=1000,
                                        value=100, step=10, key="conf_gen")
        with c2:
            max_confs = st.number_input("Manter", min_value=1, max_value=50,
                                        value=5, step=1, key="conf_keep")
        with c3:
            prune = st.number_input("RMSD (Å)", min_value=0.1, max_value=3.0,
                                    value=0.5, step=0.1, key="conf_rms")

    col_over, col_btn = st.columns([2, 3], vertical_alignment="center")
    with col_over:
        overwrite = st.toggle("Sobrescrever", key="conf_overwrite")
    with col_btn:
        clicked = st.button(f"Gerar confôrmeros ({len(valid)})", type="primary",
                            disabled=not valid, key="btn_conformers")

    if clicked:
        bar = st.progress(0.0)

        def on_progress(done: int, total: int, mol_id: str) -> None:
            bar.progress(done / total, text=mol_id)

        results = run_conformer_search(
            valid, xyz_dir, max_confs=int(max_confs), num_confs=int(num_confs),
            prune_rms=float(prune), overwrite=overwrite, on_progress=on_progress,
        )
        st.session_state["conf_results"] = results

        n_files = sum(len(r["conformers"]) for r in results if r["status"] == "ok")
        flash("success" if n_files else "warning",
              f"{n_files} confôrmero(s) salvo(s) em {xyz_dir.name}/")
        st.rerun()

    results = st.session_state.get("conf_results")
    if results:
        st.dataframe(_results_table(results), hide_index=True)