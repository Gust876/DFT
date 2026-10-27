from rdkit import Chem
from pathlib import Path
from typing import List, Tuple
import heapq


def select_low_energy(
        res, max_confs: int
        ) -> List[Tuple[int, float]]:

    converged = [
        (conf_id, energy)
        for conf_id, (status, energy) in enumerate(res)
        if status == 0
    ]

    if not converged:
        raise ValueError("Nenhum confôrmero convergiu")

    return heapq.nsmallest(
        max_confs,
        converged,
        key=lambda x: x[1]
    )

def write_xyz_files(
        mol,
        conformers,
        output_dir: Path,
        prefix: str
        ):
    output_dir.mkdir(parents=True, exist_ok=True)

    xyz_paths = []

    for rank, (conf_id, energy) in enumerate(conformers):
        filename = f"{prefix}_conf{rank}.xyz"
        filepath = output_dir / filename

        Chem.MolToXYZFile(mol, str(filepath), confId=conf_id)

        xyz_paths.append((rank, energy, filepath))

    return xyz_paths