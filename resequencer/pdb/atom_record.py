from dataclasses import dataclass


@dataclass
class AtomRecord:
    """
    Dataclass Representation of Atoms in PDB file.

    Attributes
    ----------
    line_idx : int
        Index of the line in the PDB file.
    record_name : str
        Record type name (e.g., 'ATOM' or 'HETATM').
    atom_number : int
        Unique identifier for the atom.
    atom_name : str
        Name of the atom (e.g., 'CA', 'CB').
    residue_name : str
        Three-letter code for the residue (e.g., 'ALA', 'GLY').
    chain_id : str
        Single character identifier for the protein chain.
    residue_number : int
        Residue sequence number.
    x_coord : float
        X coordinate of the atom in Angstroms.
    y_coord : float
        Y coordinate of the atom in Angstroms.
    z_coord : float
        Z coordinate of the atom in Angstroms.
    occupancy : float
        Occupancy value (0.0 to 1.0) indicating atom presence.
    b_factor : float
        Temperature factor (B-factor) indicating atom displacement.
    segment_id : str
        Segment identifier for the atom.
    element_symbol : str
        Chemical element symbol (e.g., 'C', 'N', 'O').
    charge : float
        Formal charge on the atom.
    """

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
