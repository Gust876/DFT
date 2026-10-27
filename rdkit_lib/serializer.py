from interfaces.forcefield_interface import ForceField
from pathlib import Path
import json

def serializer_json(
        xyz_data: list,
        file_output: Path,
        molecule_id: int,
        strategy: ForceField
        ):

    conformers_data = {
        "molecule_id": molecule_id,
        "conformers": [
            {
                "rank": rank, 
                "energy": energy,
                "force_field": strategy.name, 
                "xyz_file": str(path)
            }
            for rank, energy, path in xyz_data
        ]
    }

    file_output.write(json.dumps(conformers_data) + "\n")
    file_output.flush()