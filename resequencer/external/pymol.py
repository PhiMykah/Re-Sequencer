import logging
import sys
import typing
from pathlib import Path

from pymol import cmd

from resequencer.pdb import PDB, Chain

from .external import MINI_HELIX_TAIL

if typing.TYPE_CHECKING:
    from resequencer.addition import Addition


def run_pymol(
    addition: "Addition",
    pdb: PDB,
    helix_orientation: str,
    input_file: Path,
    output_path: Path,
    is_print_only: bool,
    mini_helix_length: int = MINI_HELIX_TAIL,
    alignment_type: str = "align",
) -> list[tuple]:
    """
    Execute PyMOL commands to extract and align protein chain segments based on helix orientation.
    This function extracts overlapping regions from two protein chains and performs structural
    alignment using PyMOL. The extraction is based on the helix orientation (start or end),
    and can either execute PyMOL commands directly or print them for inspection.

    Parameters
    ----------
    addition : Addition
        Addition object containing information about the chains to process.
    pdb : PDB
        PDB object representing the protein structure.
    helix_orientation : str
        Orientation of the helix ("start" or "end") determining which
        residues are extracted for alignment.
    input_file : Path
        Path to the input PDB file to be processed.
    output_path : Path
        Directory path where output PDB files (new.pdb and aligned.pdb) are saved.
    is_print_only: bool
        If True, prints PyMOL commands without executing them.
        If False, executes PyMOL commands and logs the output.
    mini_helix_length: int
        Length to overlap mini-helix, by default MINI_HELIX_TAIL
    alignment_type: str
        Type of alignment to use to align new and selection, by default 'align'

    Returns
    -------
    list[tuple,tuple]
        A list of two tuples containing the residue number ranges for target and other chains
        in the format [(target_start, target_end), (other_start, other_end)].

    Raises
    ------
    AssertionError
        If target_chain or other_chain are not Chain instances.
    """
    # Target_chain old and new
    # Chain class representation of target
    target_chain = pdb[addition.chains[0]]

    # Other chain old and new]
    # Chain class representation of other
    other_chain = pdb[addition.chains[1]]

    assert isinstance(target_chain, Chain)
    assert isinstance(other_chain, Chain)

    target_tail: int = mini_helix_length
    other_tail: int = mini_helix_length

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
        # select selection, (chain A and resi 11-12 or chain B and resi 13-14)
        cmd.select("selection", extraction)
        # # create temp, selection
        # cmd.create("temp", "selection")
        # load /PATH/TO/new.pdb
        cmd.load(str(new_path))
        # super new.pdb, selection
        if alignment_type.lower() == "super":
            cmd.super("new", "selection")
        elif alignment_type.lower() == "cealign":
            if mini_helix_length < 6 or (target_end - target_start + 1) < 6:
                cmd.align("new", "selection")
            else:
                cmd.cealign("selection", "new", window=3)
        elif alignment_type.lower() == "align":
            cmd.align("new", "selection")
        else:
            cmd.align("new", "selection")
        # deselect
        cmd.deselect()
        # delete input_file
        cmd.delete(original_obj)
        # delete temp
        # cmd.delete("temp")
        # multisave /PATH/TO/aligned.pdb
        cmd.multisave(str(aligned_path))
        cmd.delete("*")

    if alignment_type.lower() == "super":
        alignment_method = "super new, selection"
    elif alignment_type.lower() == "cealign":
        alignment_method = "cealign selection, new"
    elif alignment_type.lower() == "align":
        alignment_method = "align new, selection"
    else:
        alignment_method = "align new, selection"
    print_output = []
    original_obj = str(Path(input_file).stem)
    print_output.extend(["pymol", str(input_file), "-c", "-d"])
    command = []
    command.extend(
        [
            f"select selection, {extraction}",
            # "create temp, selection",
            f"{alignment_method}",
            "deselect",
            f"delete {original_obj}",
            # "delete temp",
            f"multisave {str(aligned_path)}",
        ]
    )
    print_output.append(f"'{';'.join(command)}'")
    logging.info("Ran pymol commands: " + " ".join(print_output))
    if is_print_only:
        print("--- pymol Commands ---", file=sys.stderr)
        print(" ".join(print_output), file=sys.stderr)
    return ranges
