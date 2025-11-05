from pathlib import Path

import pandas as pd
from biopandas.pdb import PandasPdb
from pandas import DataFrame

from resequencer.external import run_x3dna, run_pymol
from resequencer.pdb import refactor_column, update_ter

from .add import Addition, load_addition_file


def pdb_addition(
    input_file: str,
    pdb: PandasPdb,
    output: Path | str = Path.cwd().resolve(),
) -> PandasPdb:
    # Collect atoms from pdb
    atoms: DataFrame = pdb.df["ATOM"]

    # Import substitution file if it exists
    addition_file = Path(input_file)
    if not addition_file.is_file():
        raise Exception(f"Addition file '{addition_file}' does not exist!")

    additions: dict[int, Addition] = load_addition_file(addition_file)

    # ---------------------------- Create path objects --------------------------- #

    output_path: Path = (
        output if isinstance(output, Path) else Path(output).parent.resolve()
    )
    Path.mkdir(output_path, exist_ok=True, parents=True)
    new_path: Path = output_path / "new.pdb"
    aligned_path: Path = output_path / "aligned.pdb"

    # Obtain the last n-10 bases and save it as new.pdb
    for idx, addition in additions.items():
        # Determine the target chain type to modify
        target_chain: int = (
            addition.target_chain - 1 if addition.target_chain > 0 else 0
        )
        chain_type: str = addition.chains[target_chain]

        chain: DataFrame = atoms[atoms["chain_id"] == chain_type.upper()]

        # ---------------------------------------------------------------------------- #
        #                               x3DNA Mini Helix                               #
        # ---------------------------------------------------------------------------- #

        is_print_only = run_x3dna(addition, chain, chain_type, new_path)

        # ---------------------------------------------------------------------------- #
        #                                     pymol                                    #
        # ---------------------------------------------------------------------------- #

        run_pymol(
            addition, chain_type, pdb.pdb_path, new_path, aligned_path, is_print_only
        )

        if not aligned_path.exists():
            continue

        # Update working atoms and the pdb object so subsequent iterations see the change
        atoms = append_addition(atoms, addition, aligned_path)

        # change end of chains
        # find last row for chain A and chain B (case-insensitive)
        a_mask = atoms["chain_id"].astype(str).str.upper() == "A"
        b_mask = atoms["chain_id"].astype(str).str.upper() == "B"

        a_idx_list = atoms.index[a_mask].tolist()
        b_idx_list = atoms.index[b_mask].tolist()

        last_a_row = int(a_idx_list[-1]) if a_idx_list else None
        last_b_row = int(b_idx_list[-1]) if b_idx_list else None

        if last_a_row:
            update_ter(pdb, atoms.iloc[last_a_row], "A")
        if last_b_row:
            update_ter(pdb, atoms.iloc[last_b_row], "B")

    pdb.df["ATOM"] = atoms
    return pdb


def append_addition(
    atoms: DataFrame,
    addition: Addition,
    aligned_path: Path,
) -> DataFrame:
    excess_length: int = addition.total_bp - 10
    aligned_pdb = PandasPdb().read_pdb(aligned_path)
    aligned_atoms = aligned_pdb.df["ATOM"]

    # Collect unique residue numbers from the aligned atoms
    if "residue_number" not in aligned_atoms.columns:
        raise KeyError(
            "Aligned PDB DataFrame missing required column: 'residue_number'"
        )

    unique_residues: list[int] = (
        aligned_atoms["residue_number"].drop_duplicates().astype(int).tolist()
    )

    drop_count = excess_length * 2
    if drop_count <= 0:
        # nothing to drop
        target_residues = unique_residues
    else:
        # safely drop the first drop_count entries and redundant entries
        target_residues = unique_residues[drop_count + excess_length : -excess_length]

    # collect matching atom rows (ensure residue_number ints)
    mask = aligned_atoms["residue_number"].astype(int).isin(target_residues)
    collected_atoms = aligned_atoms.loc[mask].copy()

    start_pos: int = int(addition.start_position) + 1
    incremental_residues: list[int] = list(
        range(start_pos, start_pos + len(target_residues))
    )
    # create mapping from original target residues to incremental residues
    res_map = dict(zip(target_residues, incremental_residues))

    # map and replace residue numbers in collected_atoms
    mapped = collected_atoms["residue_number"].astype(int).map(res_map)
    if mapped.isnull().any():
        missing = sorted(
            set(collected_atoms["residue_number"].astype(int).unique())
            - set(res_map.keys())
        )
        raise KeyError(f"Failed to map some residue_number values: {missing}")
    collected_atoms["residue_number"] = mapped.astype(int)

    # Ensure residue_number columns are ints
    atoms["residue_number"] = atoms["residue_number"].astype(int)
    collected_atoms["residue_number"] = collected_atoms["residue_number"].astype(int)

    # Determine insertion point and how many residue indices are being inserted
    insert_after = int(addition.total_bp)
    n_insert = collected_atoms["residue_number"].nunique()

    # Shift residue numbers for atoms that come after the insertion point
    shift_mask_left = atoms["residue_number"] <= insert_after
    shift_mask_right = atoms["residue_number"] > insert_after
    if shift_mask_right.any():
        atoms.loc[shift_mask_right, "residue_number"] = (
            atoms.loc[shift_mask_right, "residue_number"] + n_insert
        )

    first_half = atoms.loc[shift_mask_left]

    # Align indices and assign new atom numbers
    collected_atoms = refactor_column(
        collected_atoms,
        "atom_number",
        int(first_half["atom_number"].iloc[-1]) + 1,
        len(collected_atoms),
    )

    collected_atoms = refactor_column(
        collected_atoms,
        "line_idx",
        int(first_half["line_idx"].iloc[-1]) + 1,
        len(collected_atoms),
    )

    second_half = atoms.loc[shift_mask_right]

    # Align indices and assign new atom numbers
    second_half = refactor_column(
        second_half,
        "atom_number",
        int(collected_atoms["atom_number"].iloc[-1]) + 1,
        len(second_half),
    )

    second_half = refactor_column(
        second_half,
        "line_idx",
        int(collected_atoms["line_idx"].iloc[-1]) + 1,
        len(second_half),
    )

    # Combine atoms
    combined = pd.concat([first_half, collected_atoms, second_half], ignore_index=True)

    combined.to_csv("output/combined.csv")
    return combined.reset_index(drop=True)
