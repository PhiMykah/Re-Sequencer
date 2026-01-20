import sys
from pathlib import Path

from pandas import DataFrame
from pymol import cmd

from resequencer.addition import Addition
from resequencer.external.x3dna.fiber_bindings import run_fiber


def run_x3dna(
    addition: Addition,
    chain: DataFrame,
    chain_type: str,
    new_path: Path,
):
    # Identify the base count for the original portion of the tail
    base_count: int = Addition.mini_helix_tail()

    res_cols = ["residue_number", "residue_name"]
    if not set(res_cols).issubset(chain.columns):
        raise KeyError(f"Chain DataFrame missing required columns: {res_cols}")

    unique_residues = chain[res_cols].drop_duplicates().reset_index(drop=True)
    # optional: as a list of tuples (residue_number, residue_name)
    unique_residue_list = [tuple(x) for x in unique_residues.to_numpy()]

    if base_count > len(unique_residue_list):
        raise ValueError(
            f"Requested base_count ({base_count}) exceeds available residues ({len(unique_residue_list)})"
        )

    # take the last or first `base_count` residues (preserves their original order)
    if addition.start_position == 1:
        tail_residues = unique_residue_list[:base_count]
    else:
        tail_residues = unique_residue_list[-base_count:]

    # sequence made from residue names (trimmed and uppercased)
    tail_sequence = "".join(
        name.strip().upper().strip("D") for _, name in tail_residues
    )
    new_sequence = "".join(a.strip().upper().strip("D") for a in addition.sequence)

    if addition.start_position == 1:
        mini_helix: str = new_sequence + tail_sequence
    else:
        mini_helix: str = tail_sequence + new_sequence

    is_rna = (
        "-rna" if ("U" in mini_helix.upper() and "T" not in mini_helix.upper()) else ""
    )
    # Build command
    fiber_cmd = ["fiber", f"-{addition.original_geometry.lower()}"]
    if is_rna:
        fiber_cmd.append(is_rna)
    fiber_cmd.append(f"-seq={mini_helix}")
    fiber_cmd.append(f"{str(new_path)}")

    try:
        print("Running fiber...", file=sys.stderr)
        run_fiber(fiber_cmd)
        is_print_only = False
    except Exception as e:
        is_print_only = True
        print(f"Error running fiber: {e}", file=sys.stderr)

    if is_print_only:
        # On Windows or error just print the commands
        print("--- x3DNA Command ---", file=sys.stderr)
        print(" ".join(fiber_cmd), file=sys.stderr)

    return is_print_only


def run_pymol(
    addition: Addition,
    input_file: str,
    new_path: Path,
    aligned_path: Path,
    print_mode: bool = False,
) -> None:
    # ----------------------------- Calculate Extract ---------------------------- #
    excess_length: int = Addition.mini_helix_tail()

    """
    Example: If chain type is A, total_bp is 12, start position is 12
        - excess_length is 12-10 = 2
        - excess_start is 12 - 2 = 10.
        - We then add 1 to not include excess start, or 10 + 1 = 11
        - Lastly, it ends at start position = 12.
    Example: If chain type is B, total_bp is 12, start position is 12,
        - excess length is 12-10 = 2, 
        - excess_start = 12
        - We then add 1 to not enclude start pos, or 12 + 1 = 13
        - excess_end = it ends at start position + excess_length = 12 + 2 = 14
    """

    extract_chain = []

    for chain_val in ("A", "B"):
        excess_start: int = (
            addition.start_position
            if chain_val == "B"
            else addition.start_position - excess_length
        )
        excess_start += 1
        excess_end: int = (
            addition.start_position + excess_length
            if chain_val == "B"
            else addition.start_position
        )
        extract_chain.append(
            " ".join(
                [
                    f"(chain {chain_val}",
                    "and",
                    "resi",
                    f"{excess_start}-{excess_end})",
                ]
            )
        )

    extraction: str = " or ".join(extract_chain)

    # --------------------------------- Run pymol -------------------------------- #
    if not print_mode:
        print("Running pymol...", file=sys.stderr)
        # load /PATH/TO/input_file.pdb
        cmd.load(input_file)
        original_obj: str = cmd.get_names("objects")[0]
        # extract temp, chain A and resi 11-12 or chain B and resi 13-14)
        cmd.extract("temp", extraction)
        # load /PATH/TO/new.pdb
        cmd.load(str(new_path))
        # delete input_file
        cmd.delete(original_obj)
        # super new.pdb, temp
        cmd.super("new", "temp")
        # multisave /PATH/TO/aligned.pdb
        cmd.multisave(str(aligned_path))
    else:
        print_output = []
        original_obj = Path(input_file).stem
        print_output.extend(["pymol", input_file, "-c", "-d"])
        command = []
        command.extend(
            [
                f"extract temp, {extraction}",
                f"load {str(new_path)}",
                f"delete {original_obj}",
                "super new, temp",
                f"multisave {str(aligned_path)}",
            ]
        )
        print_output.append(f"'{';'.join(command)}'")
        print("--- pymol Commands ---", file=sys.stderr)
        print(" ".join(print_output), file=sys.stderr)
