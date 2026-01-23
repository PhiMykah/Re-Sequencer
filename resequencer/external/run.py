import sys
from pathlib import Path

from pandas import DataFrame
from pymol import cmd

from resequencer.addition import Addition
from resequencer.external.x3dna.fiber_bindings import run_fiber


def run_x3dna(
    addition: Addition,
    target_chain: DataFrame,
    other_chain: DataFrame,
    chain_type: str,
    new_path: Path,
) -> tuple[bool, bool, tuple[int, int], tuple[int, int]]:
    # Identify the base count for the original portion of the tail
    base_count: int = Addition.mini_helix_tail()

    res_cols = ["residue_number", "residue_name"]
    if not set(res_cols).issubset(target_chain.columns):
        raise KeyError(f"Target Chain DataFrame missing required columns: {res_cols}")
    if not set(res_cols).issubset(other_chain.columns):
        raise KeyError(f"Pair Chain DataFrame missing required columns: {res_cols}")

    # Represent residues as a list of tuples (residue_number, residue_name)
    target_chain_residues: list[tuple[int, str]] = [
        tuple(x)
        for x in target_chain[res_cols]
        .drop_duplicates()
        .reset_index(drop=True)
        .to_numpy()
    ]

    # Represent pother chain's residues as a list of tuples (residue_number, residue_name)
    other_chain_residues: list[tuple[int, str]] = [
        tuple(x)
        for x in other_chain[res_cols]
        .drop_duplicates()
        .reset_index(drop=True)
        .to_numpy()
    ]

    if base_count > len(target_chain_residues):
        raise ValueError(
            f"Requested base_count ({base_count}) exceeds available residues ({len(target_chain_residues)})"
        )

    # Find smallest and largest residue number
    target_residue_numbers: list[int] = [res[0] for res in target_chain_residues]
    other_residue_numbers: list[int] = [res[0] for res in other_chain_residues]

    min_residue = min(target_residue_numbers)
    max_residue = max(target_residue_numbers)

    # take the last or first `base_count` residues (preserves their original order)
    if abs(addition.start_position - min_residue) < abs(
        addition.start_position - max_residue
    ):
        tail_residues = target_chain_residues[:base_count]
        target_range = (
            target_residue_numbers[0],
            target_residue_numbers[base_count - 1],
        )
        other_range = (
            other_residue_numbers[len(other_residue_numbers) - base_count],
            other_residue_numbers[len(other_residue_numbers) - 1],
        )
        add_from_start = True
    else:
        tail_residues = target_chain_residues[-base_count:]
        add_from_start = False
        target_range: tuple[int, int] = (
            target_residue_numbers[len(other_residue_numbers) - base_count],
            target_residue_numbers[len(other_residue_numbers) - 1],
        )
        other_range: tuple[int, int] = (
            other_residue_numbers[0],
            other_residue_numbers[base_count - 1],
        )

    # sequence made from residue names (trimmed and uppercased)
    tail_sequence = "".join(
        name.strip().upper().strip("D") for _, name in tail_residues
    )
    new_sequence = "".join(a.strip().upper().strip("D") for a in addition.sequence)

    # Check if adding new sequence to start or end of tail
    if add_from_start:
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

    return is_print_only, add_from_start, target_range, other_range


def run_pymol(
    addition: Addition,
    input_file: str,
    new_path: Path,
    aligned_path: Path,
    target_range: tuple[int, int],
    other_range: tuple[int, int],
    print_mode: bool = False,
    add_from_start: bool = False,
) -> None:
    # ----------------------------- Calculate Extract ---------------------------- #

    extract_chain = []

    for chain_val in addition.chains:
        start = (
            target_range[0]
            if chain_val.upper() == addition.target_chain
            else other_range[0]
        )
        end = (
            target_range[1]
            if chain_val.upper() == addition.target_chain
            else other_range[1]
        )
        extract_chain.append(
            " ".join(
                [
                    f"(chain {chain_val}",
                    "and",
                    "resi",
                    f"{start}-{end})",
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
        # extract temp, chain A and resi 11-12 or chain B and resi 13-14
        cmd.extract("temp", extraction)
        # load /PATH/TO/new.pdb
        cmd.load(str(new_path))
        # delete input_file
        cmd.delete(original_obj)
        # super new.pdb, temp
        if add_from_start:
            cmd.super("temp", "new")
        else:
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
                f"{'super temp, new' if add_from_start else 'super new, temp'}",
                f"multisave {str(aligned_path)}",
            ]
        )
        print_output.append(f"'{';'.join(command)}'")
        print("--- pymol Commands ---", file=sys.stderr)
        print(" ".join(print_output), file=sys.stderr)
