"""Estilos CSS globais da interface gráfica DFT."""
import streamlit as st

CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@300;400;600&display=swap');

html, body, [class*="css"] { font-family: 'IBM Plex Sans', sans-serif; }
.stApp { background-color: #0d1117; color: #e6edf3; }
.block-container { padding-top: 2rem; max-width: 1100px; }

[data-testid="stSidebar"] { background-color: #161b22; border-right: 1px solid #30363d; }

.main-title {
    font-family: 'IBM Plex Mono', monospace; font-size: 2rem; font-weight: 600;
    color: #58a6ff; letter-spacing: -0.5px; margin-bottom: 0;
}
.main-subtitle {
    font-family: 'IBM Plex Mono', monospace; font-size: 0.85rem; color: #8b949e;
    margin-top: 4px; margin-bottom: 2rem;
}

.section-card {
    background: #161b22; border: 1px solid #30363d;
    border-radius: 8px; padding: 1.5rem; margin-bottom: 1rem;
}
.section-title {
    font-family: 'IBM Plex Mono', monospace; font-size: 0.75rem; font-weight: 600;
    color: #8b949e; text-transform: uppercase; letter-spacing: 1.5px; margin-bottom: 1rem;
}

.badge { display: inline-block; padding: 2px 10px; border-radius: 20px;
    font-family: 'IBM Plex Mono', monospace; font-size: 0.72rem; font-weight: 600; }
.badge-pyscf { background: #1f3a5f; color: #58a6ff; border: 1px solid #58a6ff44; }
.badge-psi4  { background: #2d1f5f; color: #d2a8ff; border: 1px solid #d2a8ff44; }
.badge-ok    { background: #12261e; color: #3fb950; border: 1px solid #3fb95044; }

.log-box {
    background: #0d1117; border: 1px solid #30363d; border-radius: 6px;
    padding: 1rem; font-family: 'IBM Plex Mono', monospace;
    font-size: 0.78rem; color: #8b949e; max-height: 300px;
    overflow-y: auto; white-space: pre-wrap;
}

.metric-box {
    background: #161b22; border: 1px solid #30363d;
    border-radius: 8px; padding: 1rem 1.25rem; text-align: center;
}
.metric-value { font-family: 'IBM Plex Mono', monospace; font-size: 2rem; font-weight: 600; color: #58a6ff; }
.metric-label { font-size: 0.75rem; color: #8b949e; text-transform: uppercase; letter-spacing: 1px; }

/* Botões */
.stButton > button {
    background: #21262d; color: #e6edf3; border: 1px solid #30363d;
    border-radius: 8px; font-family: 'IBM Plex Mono', monospace;
    font-size: 0.85rem; font-weight: 600; padding: 0.55rem 1.5rem;
    transition: all 0.15s ease; width: 100%;
}
.stButton > button:hover { background: #30363d; border-color: #58a6ff; color: #58a6ff; }
.stButton > button[kind="primary"],
.stButton > button[data-testid="stBaseButton-primary"] {
    background: #1f6feb; border-color: #1f6feb; color: #fff;
}
.stButton > button[kind="primary"]:hover,
.stButton > button[data-testid="stBaseButton-primary"]:hover {
    background: #388bfd; border-color: #388bfd; color: #fff;
}

.stSelectbox > div > div,
.stTextInput > div > div > input {
    background: #0d1117 !important; border: 1px solid #30363d !important;
    border-radius: 8px !important; color: #e6edf3 !important;
    font-family: 'IBM Plex Mono', monospace !important;
}

.stTabs [data-baseweb="tab-list"] { background: transparent; gap: 0; border-bottom: 1px solid #30363d; }
.stTabs [data-baseweb="tab"] {
    background: transparent; color: #8b949e;
    font-family: 'IBM Plex Mono', monospace; font-size: 0.85rem;
    padding: 0.6rem 1.25rem; border: none;
}
.stTabs [aria-selected="true"] {
    color: #58a6ff !important; border-bottom: 2px solid #58a6ff !important;
    background: transparent !important;
}

hr { border-color: #30363d; }

.file-item {
    display: flex; align-items: center; gap: 8px; padding: 6px 0;
    border-bottom: 1px solid #21262d; font-family: 'IBM Plex Mono', monospace;
    font-size: 0.8rem; color: #8b949e;
}
.file-item:last-child { border-bottom: none; }
.file-name { color: #e6edf3; }

/* ---------- Novo: títulos de página ---------- */
.page-title {
    font-family: 'IBM Plex Mono', monospace; font-size: 1.15rem; font-weight: 600;
    color: #e6edf3; margin: 0.5rem 0 1.25rem 0; letter-spacing: -0.3px;
}
.page-title span { color: #58a6ff; margin-right: 8px; }

/* ---------- Novo: dropzone ---------- */
[data-testid="stFileUploaderDropzone"] {
    background: #161b22; border: 1.5px dashed #30363d; border-radius: 14px;
    padding: 2.25rem 1rem; justify-content: center; transition: all 0.15s ease;
}
[data-testid="stFileUploaderDropzone"]:hover { border-color: #58a6ff; background: #131c2b; }
[data-testid="stFileUploaderDropzone"] button {
    background: #1f6feb; border: none; color: #fff; border-radius: 8px;
}
[data-testid="stFileUploaderFile"] { display: none; }

/* ---------- Novo: chips de arquivo ---------- */
.chip-list { display: flex; flex-direction: column; gap: 6px; margin: 1rem 0 1.25rem 0; }
.chip {
    display: flex; align-items: center; gap: 12px; padding: 10px 14px;
    background: #161b22; border: 1px solid #30363d; border-radius: 10px;
    font-family: 'IBM Plex Mono', monospace; font-size: 0.8rem;
}
.chip-icon { width: 20px; height: 20px; border-radius: 50%; display: flex;
    align-items: center; justify-content: center; font-size: 0.7rem; font-weight: 700; }
.chip.ok  .chip-icon { background: #12261e; color: #3fb950; }
.chip.err { border-color: #f8514944; }
.chip.err .chip-icon { background: #2d1416; color: #f85149; }
.chip-name { color: #e6edf3; flex: 1; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
.chip-meta { color: #8b949e; font-size: 0.74rem; }
.chip.err .chip-meta { color: #f85149; }

/* ---------- Novo: fluxo do pipeline ---------- */
.flow { display: flex; align-items: center; justify-content: center; gap: 0; margin: 0.5rem 0 2rem 0; }
.flow-node {
    background: #161b22; border: 1px solid #30363d; border-radius: 12px;
    padding: 1rem 1.75rem; text-align: center; min-width: 130px;
}
.flow-num { font-family: 'IBM Plex Mono', monospace; font-size: 1.9rem; font-weight: 600; color: #58a6ff; line-height: 1.1; }
.flow-label { font-size: 0.7rem; color: #8b949e; text-transform: uppercase; letter-spacing: 1.2px; margin-top: 4px; }
.flow-link { flex: 1; min-width: 90px; text-align: center; font-family: 'IBM Plex Mono', monospace;
    font-size: 0.7rem; color: #8b949e; position: relative; padding-bottom: 10px; }
.flow-link i { display: block; height: 1px; background: #30363d; margin-top: 6px; position: relative; }
.flow-link i::after { content: ''; position: absolute; right: 0; top: -3px; width: 7px; height: 7px;
    border-top: 1px solid #58a6ff; border-right: 1px solid #58a6ff; transform: rotate(45deg); }

/* ---------- Novo: cards de ação ---------- */
.action-card {
    background: #161b22; border: 1px solid #30363d; border-radius: 12px;
    padding: 1.25rem 1.5rem; margin-bottom: 0.75rem; display: flex; align-items: center; gap: 14px;
}
.action-icon { font-size: 1.5rem; }
.action-title { font-weight: 600; font-size: 0.95rem; color: #e6edf3; }
.action-sub { font-family: 'IBM Plex Mono', monospace; font-size: 0.72rem; color: #8b949e; margin-top: 2px; }

/* ---------- Novo: radios como "pílulas" ---------- */
div[role="radiogroup"] { gap: 6px; }
div[role="radiogroup"] > label {
    background: #161b22; border: 1px solid #30363d; border-radius: 999px;
    padding: 4px 14px; transition: all 0.15s ease;
}
div[role="radiogroup"] > label:hover { border-color: #58a6ff; }

/* ---------- Minimalismo: dropzone sem texto padrão ---------- */
[data-testid="stFileUploaderDropzoneInstructions"] { display: none; }
[data-testid="stFileUploaderDropzone"] { flex-direction: column; gap: 14px; }
[data-testid="stFileUploaderDropzone"]::before { content: "↑"; font-size: 2.4rem; font-family: 'IBM Plex Mono', monospace; color: #58a6ff; line-height: 1; }

/* ---------- Fluxo com botões integrados ---------- */
.flow-node { min-width: 0; width: 100%; padding: 1.1rem 0.5rem; }
.flow-icon { font-size: 1.3rem; margin-bottom: 2px; }
.file-empty { color: #30363d; font-size: 1.4rem; font-family: 'IBM Plex Mono', monospace; padding: 6px 0; }
</style>
"""


def inject_css() -> None:
    """Injeta os estilos CSS globais na interface Streamlit."""
    st.markdown(CSS, unsafe_allow_html=True)
 