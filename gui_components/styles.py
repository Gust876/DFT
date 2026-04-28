"""Estilos CSS globais da interface gráfica DFT."""
import streamlit as st

CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@300;400;600&display=swap');

html, body, [class*="css"] { font-family: 'IBM Plex Sans', sans-serif; }

.stApp { background-color: #0d1117; color: #e6edf3; }

[data-testid="stSidebar"] {
    background-color: #161b22;
    border-right: 1px solid #30363d;
}

.main-title {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 2rem; font-weight: 600;
    color: #58a6ff; letter-spacing: -0.5px; margin-bottom: 0;
}
.main-subtitle {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 0.85rem; color: #8b949e;
    margin-top: 4px; margin-bottom: 2rem;
}

.section-card {
    background: #161b22; border: 1px solid #30363d;
    border-radius: 8px; padding: 1.5rem; margin-bottom: 1rem;
}
.section-title {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 0.75rem; font-weight: 600; color: #8b949e;
    text-transform: uppercase; letter-spacing: 1.5px; margin-bottom: 1rem;
}

.badge { display: inline-block; padding: 2px 10px; border-radius: 20px;
    font-family: 'IBM Plex Mono', monospace; font-size: 0.72rem; font-weight: 600; }
.badge-pyscf { background: #1f3a5f; color: #58a6ff; border: 1px solid #58a6ff44; }
.badge-psi4  { background: #2d1f5f; color: #d2a8ff; border: 1px solid #d2a8ff44; }

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
.metric-value {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 2rem; font-weight: 600; color: #58a6ff;
}
.metric-label { font-size: 0.75rem; color: #8b949e; text-transform: uppercase; letter-spacing: 1px; }

.stButton > button {
    background: #21262d; color: #e6edf3; border: 1px solid #30363d;
    border-radius: 6px; font-family: 'IBM Plex Mono', monospace;
    font-size: 0.85rem; font-weight: 600; padding: 0.5rem 1.5rem;
    transition: all 0.15s ease; width: 100%;
}
.stButton > button:hover { background: #30363d; border-color: #58a6ff; color: #58a6ff; }

div[data-testid="column"]:first-child .stButton > button {
    background: #1f6feb; border-color: #1f6feb; color: white;
}
div[data-testid="column"]:first-child .stButton > button:hover {
    background: #388bfd; border-color: #388bfd; color: white;
}

.stSelectbox > div > div,
.stTextInput > div > div > input {
    background: #0d1117 !important; border: 1px solid #30363d !important;
    border-radius: 6px !important; color: #e6edf3 !important;
    font-family: 'IBM Plex Mono', monospace !important;
}

.stTabs [data-baseweb="tab-list"] { background: transparent; gap: 0; border-bottom: 1px solid #30363d; }
.stTabs [data-baseweb="tab"] {
    background: transparent; color: #8b949e;
    font-family: 'IBM Plex Mono', monospace; font-size: 0.85rem;
    padding: 0.5rem 1.25rem; border: none;
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
</style>
"""


def inject_css() -> None:
    """Injeta os estilos CSS globais na interface Streamlit."""
    st.markdown(CSS, unsafe_allow_html=True)
