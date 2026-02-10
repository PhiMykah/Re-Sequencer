from dataclasses import dataclass


@dataclass
class PDBRecord:
    line_idx: int
    record_name: str
    atom_number: int
    atom_name: str
    residue_name: str
    chain_id: str
    residue_number: int
    x_coord: float
    y_coord: float
    z_coord: float
    occupancy: float
    b_factor: float
    segment_id: str
    element_symbol: str
    charge: float
