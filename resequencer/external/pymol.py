import sys
import typing
from pathlib import Path

from pymol import cmd

from .external import mini_helix_tail
from resequencer.pdb import PDB, Chain

if typing.TYPE_CHECKING:
    from resequencer.addition import Addition


def run_pymol(
    addition: "Addition",
    pdb: PDB,
    mini_helix: str,
    helix_orientation: str,
    input_file: Path,
    output_path: Path,
    is_print_only: bool,
) -> list[tuple]:

    # Target_chain old and new
    # Chain class representation of target
    target_chain = pdb[addition.chains[0]]

    # Other chain old and new]
    # Chain class representation of other
    other_chain = pdb[addition.chains[1]]

    assert isinstance(target_chain, Chain)
    assert isinstance(other_chain, Chain)

    target_tail: int = mini_helix_tail()
    other_tail: int = mini_helix_tail()

    # Collect the overlap from the helix
    if helix_orientation.lower() == "start":
        target_start: int = target_chain.starting_residue
        target_end: int = target_start + target_tail - 1

        other_start: int = other_chain[-1].residue_number - other_tail + 1  # type: ignore
        other_end: int = other_chain[-1].residue_number  # type: ignore
    else:
        target_start: int = target_chain[-1].residue_number - target_tail + 1  # type: ignore
        target_end: int = target_chain[-1].residue_number  # type: ignore

        other_start: int = other_chain.starting_residue
        other_end: int = other_start + other_tail - 1

    ranges: list[tuple] = [(target_start, target_end), (other_start, other_end)]
    extract_chain = []
    for idx, chain in enumerate(addition.chains):
        # Assumes addition.chains has 2 elements
        extract_chain.append(
            f"(chain {chain.upper()} and resi {ranges[idx][0]}-{ranges[idx][1]})"
        )

    extraction: str = " or ".join(extract_chain)

    new_path: Path = output_path / "new.pdb"
    aligned_path = output_path / "aligned.pdb"

    # --------------------------------- Run pymol -------------------------------- #
    if not is_print_only:
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
        if helix_orientation.lower() == "start":
            cmd.super("new", "temp")
        else:
            cmd.super("temp", "new")
        # multisave /PATH/TO/aligned.pdb
        cmd.multisave(str(aligned_path))
    else:
        print_output = []
        original_obj = str(Path(input_file).stem)
        print_output.extend(["pymol", str(input_file), "-c", "-d"])
        command = []
        command.extend(
            [
                f"extract temp, {extraction}",
                f"load {str(new_path)}",
                f"delete {original_obj}",
                f"{'super temp, new' if helix_orientation.lower() == 'start' else 'super new, temp'}",
                f"multisave {str(aligned_path)}",
            ]
        )
        print_output.append(f"'{';'.join(command)}'")
        print("--- pymol Commands ---", file=sys.stderr)
        print(" ".join(print_output), file=sys.stderr)

    if not is_print_only:
        print_output = []
        original_obj = str(Path(input_file).stem)
        print_output.extend(["pymol", str(input_file), "-c", "-d"])
        command = []
        command.extend(
            [
                f"extract temp, {extraction}",
                f"load {str(new_path)}",
                f"delete {original_obj}",
                f"{'super temp, new' if helix_orientation.lower() == 'start' else 'super new, temp'}",
                f"multisave {str(aligned_path)}",
            ]
        )
        print_output.append(f"'{';'.join(command)}'")
        print("Ran pymol commands:", " ".join(print_output), file=sys.stderr)

    return ranges
