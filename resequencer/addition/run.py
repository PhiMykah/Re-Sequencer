from pathlib import Path

import pandas as pd
from biopandas.pdb import PandasPdb
from pandas import DataFrame

from resequencer.external import run_x3dna, run_pymol
from resequencer.pdb import refactor_column, update_ter, update_hetatm

from .add import Addition, load_addition_file


def pdb_addition(
    input_file: str,
    pdb: PandasPdb,
    output_path: Path | str = Path.cwd().resolve(),
) -> PandasPdb:
    # Collect atoms from pdb
    atoms: DataFrame = pdb.df["ATOM"]

    # Import substitution file if it exists
    addition_file = Path(input_file)
    if not addition_file.is_file():
        raise Exception(f"Addition file '{addition_file}' does not exist!")

    additions: dict[int, Addition] = load_addition_file(addition_file)

    # ---------------------------- Create path objects --------------------------- #

    output_path = Path(output_path)
    output_dir = output_path.parent if output_path.is_file() else output_path
    Path.mkdir(output_dir, exist_ok=True, parents=True)
    new_path: Path = output_dir / "new.pdb"
    aligned_path: Path = output_dir / "aligned.pdb"

    for idx, addition in additions.items():
        # Determine the target chain type to modify
        if addition.target_chain not in addition.chains:
            raise TypeError(
                f"Target chain {addition.target_chain} is not"
                " in list of provided chains {','.join(addition.chains)}"
            )
        chain_type: str = addition.target_chain

        # Set the other_chain type to the chain that isn't target
        other_chain_type = (
            addition.chains[0]
            if addition.chains[1].upper() == addition.target_chain
            else addition.chains[1]
        )
        target_chain: DataFrame = atoms[atoms["chain_id"] == chain_type.upper()]
        other_chain: DataFrame = atoms[atoms["chain_id"] == other_chain_type.upper()]

        if target_chain.empty:
            raise ValueError(f"Target chain '{chain_type}' is not in provided pdb!")

        # ---------------------------------------------------------------------------- #
        #                               x3DNA Mini Helix                               #
        # ---------------------------------------------------------------------------- #

        is_print_only, add_from_start, target_range, other_range = run_x3dna(
            addition, target_chain, other_chain, chain_type, new_path
        )

        # ---------------------------------------------------------------------------- #
        #                                     pymol                                    #
        # ---------------------------------------------------------------------------- #

        run_pymol(
            addition,
            pdb.pdb_path,
            new_path,
            aligned_path,
            target_range,
            other_range,
            is_print_only,
            add_from_start,
        )

        if not aligned_path.exists():
            continue

        # Update working atoms and the pdb object so subsequent iterations see the change
        atoms = append_addition(atoms, addition, aligned_path, add_from_start)

        # change end of chains
        # find last row for chain A and chain B (case-insensitive)
        a_mask = atoms["chain_id"].astype(str).str.upper() == addition.chains[0].upper()
        b_mask = atoms["chain_id"].astype(str).str.upper() == addition.chains[1].upper()

        a_idx_list = atoms.index[a_mask].tolist()
        b_idx_list = atoms.index[b_mask].tolist()

        last_a_row = int(a_idx_list[-1]) if a_idx_list else None
        last_b_row = int(b_idx_list[-1]) if b_idx_list else None

        if last_a_row:
            update_ter(pdb, atoms.iloc[last_a_row], addition.chains[0].upper())
        if last_b_row:
            update_ter(pdb, atoms.iloc[last_b_row], addition.chains[1].upper())

    pdb.df["ATOM"] = atoms
    last_line_idx = update_hetatm(pdb)
    others = pdb.df["OTHERS"]
    if last_line_idx is not None and others is not None and not others.empty:
        keys_order = ["ENDMDL", "MASTER", "END"]
        for key in keys_order:
            # build per-key mask based on available columns
            if "record_name" in others.columns:
                per_mask = others["record_name"].astype(str) == key
            elif "line" in others.columns:
                per_mask = others["line"].astype(str).str.split().str[0] == key
            else:
                per_mask = pd.Series(False, index=others.index)

            if not per_mask.any():
                continue

            new_idx = int(last_line_idx) + 1

            if "line_idx" in others.columns:
                others.loc[per_mask, "line_idx"] = new_idx
            elif "line_number" in others.columns:
                others.loc[per_mask, "line_number"] = new_idx
            elif "line" in others.columns:

                def _replace_idx(s, idx):
                    s = str(s).rstrip()
                    parts = s.split()
                    if parts and parts[-1].isdigit():
                        parts[-1] = str(idx)
                        return " ".join(parts)
                    return s + " " + str(idx)

                others.loc[per_mask, "line"] = others.loc[per_mask, "line"].apply(
                    lambda s: _replace_idx(s, new_idx)
                )

            # increment for the next key
            last_line_idx = new_idx
        pdb.df["OTHERS"] = others
    return pdb


def append_addition(
    atoms: DataFrame,
    addition: Addition,
    aligned_path: Path,
    add_from_start: bool,
) -> DataFrame:
    excess_length: int = Addition.mini_helix_tail()
    aligned_pdb = PandasPdb().read_pdb(aligned_path)
    aligned_atoms = aligned_pdb.df["ATOM"]

    # Collect unique residue numbers from the aligned atoms
    if "residue_number" not in aligned_atoms.columns:
        raise KeyError(
            "Aligned PDB DataFrame missing required column: 'residue_number'"
        )

    res_numbers = aligned_atoms["residue_number"]
    unique_residues: list[int] = res_numbers[
        res_numbers.ne(res_numbers.shift())
    ].tolist()

    # Extract unique residues from tail
    drop_count = excess_length * 2
    if drop_count <= 0:
        # nothing to drop
        target_residues = unique_residues
    else:
        # safely drop the first drop_count entries and redundant entries
        target_residues = unique_residues[drop_count + excess_length : -excess_length]

    # collect matching atom rows (ensure residue_number ints)
    mask = aligned_atoms["residue_number"].astype(int).isin(target_residues)
    collected_atoms: DataFrame = aligned_atoms.loc[mask].copy()

    # Begin tail-end or front-end addition
    if add_from_start:
        start_pos: int = 1
    else:
        start_pos: int = int(addition.start_position) + 1
    incremental_residues: list[int] = list(
        range(start_pos, start_pos + len(target_residues))
    )
    # create mapping from original target residues to incremental residues
    res_map = dict(zip(target_residues, incremental_residues))
    # create mapping from original chain to specific chain
    chain_ids = collected_atoms["chain_id"].unique().tolist()
    chain_pairs = []

    # NOTE: Assumes that the chain ids are in the same order as the addition chains!!!
    for idx, chain_id in enumerate(chain_ids):
        chain_pairs.append((chain_id, addition.chains[idx]))
    chain_map = dict(chain_pairs)

    # ------------ map and replace residue numbers in collected_atoms ------------ #
    residue_mapped = collected_atoms["residue_number"].astype(int).map(res_map)
    if residue_mapped.isnull().any():
        missing = sorted(
            set(collected_atoms["residue_number"].astype(int).unique())
            - set(res_map.keys())
        )
        raise KeyError(f"Failed to map some residue_number values: {missing}")
    collected_atoms["residue_number"] = residue_mapped.astype(int)

    # Ensure residue_number columns are ints
    atoms["residue_number"] = atoms["residue_number"].astype(int)
    collected_atoms["residue_number"] = collected_atoms["residue_number"].astype(int)

    # ---------------- map and replace chain id in collected_atoms --------------- #
    chain_mapped = collected_atoms["chain_id"].astype(str).map(chain_map)
    if chain_mapped.isnull().any():
        missing = sorted(
            set(collected_atoms["chainX_id"].astype(str).unique())
            - set(chain_map.keys())
        )
        raise KeyError(f"Failed to map some chain_id values: {missing}")
    collected_atoms["chain_id"] = chain_mapped.astype(str)

    # Ensure chain_id columns are str
    atoms["chain_id"] = atoms["chain_id"].astype(str)
    collected_atoms["chain_id"] = collected_atoms["chain_id"].astype(str)

    # Determine insertion point and how many residue indices are being inserted
    insert_after = 0 if add_from_start else int(addition.start_position)
    n_insert = collected_atoms["residue_number"].nunique()

    # Shift residue numbers for atoms that come after the insertion point
    shift_mask_left = atoms["residue_number"] <= insert_after
    shift_mask_right = atoms["residue_number"] > insert_after
    if shift_mask_right.any():
        atoms.loc[shift_mask_right, "residue_number"] = (
            atoms.loc[shift_mask_right, "residue_number"] + n_insert
        )

    # ---------------------------------------------------------------------------- #
    #                             Stitch data together                             #
    # ---------------------------------------------------------------------------- #

    first_half = atoms.loc[shift_mask_left]

    # Re-index by aligning indices and assign new atom numbers
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

    # combined.to_csv("output/combined.csv")
    return combined.reset_index(drop=True)
