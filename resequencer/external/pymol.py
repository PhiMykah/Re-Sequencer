import sys
from resequencer.addition import Addition
from pathlib import Path
# conda install -c conda-forge -c schrodinger pymol-bundle, requires python 3.10
from pymol import cmd

def run_pymol(
    addition: Addition,
    chain_type: str,
    input_file: str,
    new_path: Path,
    aligned_path: Path,
    print_mode: bool = False,
) -> None:
    # ----------------------------- Calculate Extract ---------------------------- #
    excess_length: int = addition.total_bp - 10

    # NOTE: Determine if this is the proper method to calculate the resi range
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
