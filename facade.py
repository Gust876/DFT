"""Ponto de entrada da interface gráfica do pipeline DFT.

Orquestra a renderização da sidebar, das abas de execução,
arquivos e logs, delegando cada responsabilidade ao módulo
correspondente em components/ e pages/.

Usage:
    streamlit run facade.py
"""
import streamlit as st
from pathlib import Path

from gui_components.styles import inject_css
from gui_components.sidebar import render_sidebar
from gui_pages.execution import render_execucao
from gui_pages.files import render_arquivos
from gui_pages.logs import render_logs


st.set_page_config(
    page_title="DFT Automation",
    page_icon="⚛️",
    layout="wide",
    initial_sidebar_state="expanded"
)

inject_css()


BASE_DIR   = Path(__file__).parent
XYZ_DIR    = BASE_DIR / "xyz_semi_opt"
OPT_DIR    = BASE_DIR / "xyz_opt"
OUTPUT_DIR = BASE_DIR / "zip_dir"
LOGS_DIR   = BASE_DIR / "logs"

for d in [XYZ_DIR, OPT_DIR, OUTPUT_DIR, LOGS_DIR]:
    d.mkdir(parents=True, exist_ok=True)


config = render_sidebar(XYZ_DIR, OPT_DIR, OUTPUT_DIR)


tab1, tab2, tab3 = st.tabs(["▶  Execução", "▶  Arquivos", "▶  Logs"])

with tab1:
    render_execucao(config, BASE_DIR)

with tab2:
    render_arquivos(XYZ_DIR, OPT_DIR, OUTPUT_DIR)

with tab3:
    render_logs(LOGS_DIR)
