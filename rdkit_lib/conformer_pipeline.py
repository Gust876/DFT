"""Busca conformacional (RDKit) integrada ao pipeline DFT.

Reaproveita, sem modificá-los, os módulos do projeto Conformers-RDKit
(ConformerGenerator, estratégias de campo de força, seleção por energia,
escrita de .xyz e serialização). Esta camada apenas reorganiza o fluxo do
`workflow.py` original em funções reutilizáveis pela GUI, gravando os
confôrmeros em `xyz_semi_opt/`, a entrada do pipeline de otimização.
"""
import json
import logging
import re
import shutil
import tempfile
from pathlib import Path
from typing import Callable, Optional

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from rdkit_lib.forcefields import MMFF94Strategy, UFFStrategy
from rdkit_lib.generate_conformers import ConformerGenerator
from rdkit_lib.utils_rdkit import select_low_energy, write_xyz_files
from rdkit_lib.serializer import serializer_json

_UNSAFE = re.compile(r"[^\w\-.]")


def _get_logger(logs_dir: Path) -> logging.Logger:
    """Logger dedicado, gravando em logs/conformers_process.log."""
    logger = logging.getLogger("dft.conformers")
    if not logger.handlers:
        logs_dir.mkdir(parents=True, exist_ok=True)
        handler = logging.FileHandler(logs_dir / "conformers_process.log")
        handler.setFormatter(
            logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        )
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        logger.propagate = False
    return logger


def _unique_ids(entries: list[dict]) -> list[dict]:
    """Sanitiza os ids para uso em nomes de arquivo e remove duplicidades."""
    used = set()
    for entry in entries:
        base = _UNSAFE.sub("_", str(entry["id"]).strip()) or "mol"
        candidate, k = base, 2
        while candidate in used:
            candidate, k = f"{base}_{k}", k + 1
        used.add(candidate)
        entry["id"] = candidate
    return entries


def parse_entries(text: str) -> list[dict]:
    """Interpreta linhas 'nome SMILES' ou apenas 'SMILES' (linhas '#' ignoradas)."""
    entries = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) == 1:
            mol_id, smiles = f"mol_{len(entries) + 1}", parts[0]
        else:
            mol_id, smiles = parts[0], parts[1]
        entries.append({"id": mol_id, "canonical_smiles": smiles})
    return _unique_ids(entries)


def parse_json_entries(raw) -> list[dict]:
    """Lê o formato do projeto de confôrmeros: lista de {"id", "canonical_smiles"}."""
    data = json.loads(raw)
    if isinstance(data, dict):
        data = [data]
    if not isinstance(data, list):
        raise ValueError("o JSON deve ser uma lista de objetos")

    entries = []
    for i, item in enumerate(data, start=1):
        if not isinstance(item, dict):
            raise ValueError(f"item {i}: esperado um objeto")
        smiles = item.get("canonical_smiles") or item.get("smiles")
        if not smiles:
            raise ValueError(f"item {i}: campo 'canonical_smiles' ausente")
        entries.append({
            "id": item.get("id", f"mol_{i}"),
            "canonical_smiles": str(smiles).strip(),
        })
    return _unique_ids(entries)


def describe_smiles(smiles: str) -> tuple[bool, str]:
    """Valida um SMILES e retorna (válido, fórmula molecular ou erro)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "SMILES inválido"
    return True, rdMolDescriptors.CalcMolFormula(mol)


def _search_molecule(entry: dict, tmp_dir: Path, num_confs: int,
                     prune_rms: float, max_confs: int):
    """Gera e seleciona confôrmeros; tenta MMFF94 e, se falhar, UFF."""
    molecule = Chem.MolFromSmiles(entry["canonical_smiles"])
    if molecule is None:
        raise ValueError("SMILES inválido")

    generator = ConformerGenerator(mol=molecule, num_confs=num_confs, prune_rms=prune_rms)
    last_error = None

    for strategy in (MMFF94Strategy(), UFFStrategy()):
        try:
            mol_3d, res = generator.run_optimization(strategy)
            best = select_low_energy(res, max_confs=max_confs)
            xyz_data = write_xyz_files(mol_3d, best, tmp_dir, entry["id"])
            return strategy, xyz_data
        except Exception as error:
            last_error = error

    raise RuntimeError(f"nenhum campo de força convergiu ({last_error})")


def _clear_previous(output_dir: Path, mol_id: str) -> None:
    """Remove confôrmeros antigos da mesma molécula (apenas {id}_conf<N>.xyz)."""
    pattern = re.compile(rf"^{re.escape(mol_id)}_conf\d+\.xyz$")
    for file in output_dir.iterdir():
        if file.is_file() and pattern.match(file.name):
            file.unlink()


def run_conformer_search(
        entries: list[dict],
        output_dir: Path,
        max_confs: int = 5,
        num_confs: int = 100,
        prune_rms: float = 0.5,
        overwrite: bool = False,
        on_progress: Optional[Callable[[int, int, str], None]] = None,
        logs_dir: Path = Path("./logs"),
    ) -> list[dict]:
    """Executa a busca conformacional e grava os .xyz em output_dir.

    Retomada automática: moléculas cujo `{id}_conf0.xyz` já existe em
    output_dir são ignoradas, salvo se overwrite=True.

    Args:
        entries: Lista de {"id", "canonical_smiles"}.
        output_dir: Destino dos .xyz (normalmente xyz_semi_opt/).
        max_confs: Quantos confôrmeros de menor energia manter por molécula.
        num_confs: Quantos confôrmeros gerar antes da seleção.
        prune_rms: Limiar de RMSD (Å) para descartar confôrmeros redundantes.
        overwrite: Se True, substitui os confôrmeros já existentes.
        on_progress: Callback (concluídas, total, id_atual).
        logs_dir: Diretório dos logs e do backup JSONL de energias.

    Returns:
        Uma lista de dicionários por molécula, com 'id', 'status'
        ('ok' | 'skipped' | 'failed'), 'message' e, se 'ok', 'force_field'
        e 'conformers' (lista de {'rank', 'file', 'energy', 'rel_energy'}
        em kcal/mol).
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    logger = _get_logger(logs_dir)
    results = []

    with tempfile.TemporaryDirectory() as tmp, \
            open(logs_dir / "conformers_backup.jsonl", "a") as backup:
        tmp_dir = Path(tmp)

        for i, entry in enumerate(entries, start=1):
            mol_id = entry["id"]
            try:
                if not overwrite and (output_dir / f"{mol_id}_conf0.xyz").exists():
                    logger.warning(f"{mol_id}: já existe, ignorada")
                    results.append({"id": mol_id, "status": "skipped", "message": "já existe"})
                else:
                    strategy, xyz_data = _search_molecule(
                        entry, tmp_dir, num_confs, prune_rms, max_confs
                    )
                    _clear_previous(output_dir, mol_id)

                    final = []
                    for rank, energy, src in xyz_data:
                        dst = output_dir / src.name
                        shutil.copy2(src, dst)
                        final.append((rank, energy, dst))

                    serializer_json(
                        xyz_data=final, file_output=backup,
                        molecule_id=mol_id, strategy=strategy
                    )
                    e_min = min(energy for _, energy, _ in final)
                    results.append({
                        "id": mol_id, "status": "ok", "message": "ok",
                        "force_field": strategy.name,
                        "conformers": [
                            {"rank": rank, "file": dst.name,
                             "energy": energy, "rel_energy": energy - e_min}
                            for rank, energy, dst in final
                        ],
                    })
                    logger.info(f"{mol_id}: {len(final)} confôrmero(s) com {strategy.name}")

            except Exception as error:
                logger.error(f"{mol_id} falhou: {error}")
                results.append({"id": mol_id, "status": "failed", "message": str(error)})

            if on_progress:
                on_progress(i, len(entries), mol_id)

    return results