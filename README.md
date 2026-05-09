# DFT Automation Pipeline

Pipeline de automação para cálculos de Teoria do Funcional da Densidade (DFT), desenvolvido como parte de uma Iniciação Científica (PIVIC) no Departamento de Química. O projeto automatiza otimização geométrica e cálculo de frequências vibracionais de moléculas, com suporte às engines **PySCF** e **Psi4**, execução paralela e interface gráfica via Streamlit.

---

## Funcionalidades

- Otimização geométrica de conjuntos de moléculas em paralelo
- Cálculo de frequências vibracionais para validação de estados fundamentais
- Geração automática de arquivos `.molden` e `.cube` para moléculas em estado fundamental
- Compressão dos resultados em arquivos `.zip`
- Suporte às engines PySCF e Psi4, intercambiáveis via linha de comando ou interface gráfica
- Interface gráfica (Streamlit) modularizada com visualização de logs e status dos arquivos em tempo real
- Retomada automática de cálculos interrompidos

---

## Estrutura do projeto

```
DFT/
├── docker/
│   ├── Dockerfile
│   └── .dockerignore
├── engines/
│   ├── pyscf/
│   │   ├── engine_pyscf.py   # Engine DFT Kohn-Sham com density fitting
│   │   ├── mol.py            # Construção da molécula PySCF
│   │   ├── optimizer.py      # Otimização geométrica PySCF (geomeTRIC)
│   │   └── frequency.py      # Frequências vibracionais PySCF (Hessiana)
│   ├── psi4/
│   │   ├── mol.py            # Construção da molécula Psi4
│   │   ├── optimizer.py      # Otimização geométrica Psi4
│   │   └── frequency.py      # Frequências vibracionais Psi4
│   └── utils.py              # Funções factory de instanciação
├── interfaces/
│   ├── engine_interface.py   # ABC para engines DFT
│   ├── molecule_interface.py # ABC para factories de moléculas
│   ├── optimizer_interface.py# ABC para otimizadores
│   ├── frequency_interface.py# ABC para calculadores de frequência
│   ├── result_interface.py   # ABC para resultados
│   └── strategy_interface.py # ABC para estratégias de engine
├── strategies/
│   ├── strategies.py         # PySCFStrategy e Psi4Strategy
│   └── run_strategy.py       # Contextos SetOptimizerStrategy e SetFrequencyStrategy
├── parallel/
│   ├── optimizer_parallel.py # Função geometry_optimizer (execução paralela)
│   └── frequency_parallel.py # Função frequency (execução paralela)
├── results/
│   ├── pyscf_result.py       # PySCFResult: escrita de .molden e .cube
│   └── psi4_result.py        # Psi4Result: escrita de .molden e .cube
├── workflows/
│   ├── workflow_optimizer.py       # Workflow de otimização (terminal)
│   ├── workflow_frequencies.py     # Workflow de frequências (terminal)
│   ├── workflow_optimizer_gui.py   # Workflow de otimização (GUI)
│   └── workflow_frequencies_gui.py # Workflow de frequências (GUI)
├── utils/
│   └── prompt.py             # Coleta de parâmetros (terminal e GUI)
├── gui_components/
│   ├── styles.py             # CSS global da interface gráfica
│   └── sidebar.py            # Sidebar de configuração
├── gui_pages/
│   ├── execution.py          # Aba de execução dos workflows
│   ├── files.py              # Aba de listagem de arquivos
│   └── logs.py               # Aba de visualização de logs
├── facade.py                 # Ponto de entrada da interface gráfica (Streamlit)
├── main.py                   # Ponto de entrada CLI
├── xyz_semi_opt/             # Geometrias de entrada (.xyz semi-otimizados)
├── xyz_opt/                  # Geometrias otimizadas (.xyz)
├── zip_dir/                  # Resultados finais (.zip com .molden e .cube)
└── logs/                     # Logs de execução
```

---

## Pré-requisitos

- Linux (ou Windows com WSL/Dev Container)
- Conda (recomendado)
- PySCF e Psi4 instalados no ambiente

O projeto está configurado para uso com **Dev Containers** no VS Code, o que facilita a execução em Windows, onde PySCF e Psi4 não têm suporte nativo.

---

## Instalação

```bash
# 1. Clone o repositório
git clone https://github.com/Gust876/DFT.git
cd DFT

# 2. Crie o ambiente conda (ou use o Dev Container)
conda env create -f environment.yml
conda activate dft
```

---

## Como usar

### Via linha de comando

```bash
# Otimização geométrica
python main.py optimizer

# Frequências vibracionais
python main.py frequencies
```

O programa solicitará interativamente a escolha da engine, basis set e funcional de troca-correlação.

**Formato de entrada:** coloque os arquivos `.xyz` das moléculas em `xyz_semi_opt/` antes de executar.

### Via interface gráfica

```bash
streamlit run facade.py
```

Acesse `http://localhost:8501` no navegador. Configure a engine, basis set e funcional na barra lateral e clique no botão do workflow desejado.

---

## Fluxo de cálculo

```
xyz_semi_opt/          →   Otimização Geométrica   →   xyz_opt/
(geometrias de entrada)    (PySCF ou Psi4)             (geometrias otimizadas)

xyz_opt/               →   Frequências Vibracionais →   zip_dir/
(geometrias otimizadas)    (PySCF ou Psi4)             (.molden + .cube comprimidos)
```

Moléculas com frequências imaginárias (não pertencentes a um mínimo de energia) são automaticamente excluídas da geração de arquivos de saída.

---

## Engines suportadas

| Engine | Otimização | Frequências | Molden | Cube |
|--------|-----------|-------------|--------|------|
| PySCF  | ✅        | ✅          | ✅     | ✅   |
| Psi4   | ✅        | ✅          | ✅     | ✅   |

---

## Arquitetura

O projeto utiliza os padrões de projeto **Factory Method** e **Strategy**, organizados em camadas bem definidas:

- **interfaces/** — contratos abstratos (ABCs) para todas as entidades do sistema
- **engines/** — implementações concretas de PySCF e Psi4
- **strategies/** — orquestração do fluxo de execução por engine
- **parallel/** — funções otimizadas para execução paralela via Joblib
- **workflows/** — coordenação de alto nível dos cálculos (terminal e GUI)
- **results/** — exportação dos resultados em formatos de química computacional
- **gui_components/** — componentes reutilizáveis da interface gráfica
- **gui_pages/** — páginas da interface gráfica (execução, arquivos, logs)

A adição de uma nova engine (ex: ORCA, Gaussian) requer apenas a criação de novos módulos em `engines/` e uma nova classe em `strategies/strategies.py`, bem como a implementação de uma nova classe em `results/` para encapsular os resultados da engine, sem modificar nenhum outro componente do sistema.

---

## Interface gráfica

A interface gráfica é construída com **Streamlit** e organizada em três abas:

- **Execução** — seleção de engine, basis set e funcional; execução dos workflows com feedback em tempo real
- **Arquivos** — listagem dos arquivos de entrada, geometrias otimizadas e resultados gerados
- **Logs** — visualização dos logs de execução por workflow

A GUI se comunica com os workflows via variáveis de ambiente (`DFT_ENGINE`, `DFT_BASIS`, `DFT_XC`, `DFT_NJOBS`), mantendo os workflows completamente independentes da interface.

---

## Logs

Os logs de execução são salvos automaticamente em `logs/`:

- `parallel_optimization.log` — progresso da otimização geométrica
- `parallel_frequencies.log` — progresso do cálculo de frequências
- `frequencies_process.log` — detalhes do processamento por molécula

---

## Licença

Projeto desenvolvido para fins de pesquisa acadêmica (PIVIC) — Departamento de Química da Universidade Federal da Paraíba.
