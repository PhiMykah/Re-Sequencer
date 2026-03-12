import logging
import sys
import typing
from pathlib import Path

from resequencer.external.x3dna.fiber_bindings import run_fiber
from resequencer.pdb import PDB

from .external import MINI_HELIX_TAIL

if typing.TYPE_CHECKING:
    from resequencer.addition import Addition


def run_x3dna(addition: "Addition", pdb: PDB, output_path: Path):
    """
    Generate a mini-helix structure using x3DNA's fiber command based on sequence changes.
    This function analyzes the differences between old and new DNA/RNA sequences,
    extracts a relevant portion to create a mini-helix, and uses the fiber tool to
    generate a 3D structure. If the sequences match completely in the overlap region,
    the mini-helix is extracted from the tail; otherwise, it's extracted from the start.

    Parameters
    ----------
    addition : Addition
        An Addition object containing old and new sequences for both
        target and other chains, as well as the original geometry.
    pdb : PDB
        A PDB object used to determine if the target chain is DNA or RNA.
    output_path : Path
        The directory path where the generated PDB file will be saved.

    Returns
    -------
    tuple
        A tuple containing:
        - mini_helix (str): The sequence string used to generate the fiber structure.
        - helix_orientation (str): Either "start" or "end" indicating where the differences
                                    in the sequence begin relative to the mini-helix.
        - is_print_only (bool): True if fiber command failed and was only printed (e.g., on Windows),
                                False if the command executed successfully.

    Raises
    ------
    ValueError
        If the new target chain is smaller than or equal to the original target chain,
        or if the new other chain is smaller than or equal to the original other chain.
    """

    # Identify the base count for the original portion of the tail
    base_count: int = MINI_HELIX_TAIL

    # Target_chain old and new
    target_chain: list[str] = addition.old_seq[0]
    new_target_chain: list[str] = addition.new_seq[0]

    # Other chain old and new
    other_chain: list[str] = addition.old_seq[1]
    new_other_chain: list[str] = addition.new_seq[1]

    # Find where sequences differ (first mismatch in the overlap)
    if len(new_target_chain) <= len(target_chain):
        raise ValueError(
            f"New chain ({len(new_target_chain)}) is smaller than or equal to original chain ({len(target_chain)})"
        )
    if len(new_other_chain) <= len(other_chain):
        raise ValueError(
            f"New chain ({len(new_other_chain)}) is smaller than or equal to original chain ({len(other_chain)})"
        )
    chain_len: int = len(target_chain)
    diff_index: int = next(
        (i for i in range(chain_len) if target_chain[i] != new_target_chain[i]),
        -1,
    )

    # Build the mini-helix from the new target sequence.
    # If the overlap matches entirely, the mini-helix starts at the tail.
    if diff_index == -1:
        helix_orientation = "end"
        mini_helix_start: int = max(chain_len - base_count, 0)
        mini_helix_list: list[str] = new_target_chain[mini_helix_start:]
    else:
        helix_orientation = "start"
        mini_helix_end: int = len(new_target_chain) - len(target_chain) + base_count
        mini_helix_list: list[str] = new_target_chain[:mini_helix_end]

    mini_helix: str = "".join(
        atom.strip().upper().strip("D") for atom in mini_helix_list
    )

    # Determine if it is RNA or DNA
    is_rna: str = "" if pdb.is_dna(addition.target_chain) else "-rna"

    fiber_cmd = ["fiber", f"-{addition.original_geometry.lower()}"]
    if is_rna:
        fiber_cmd.append(is_rna)
    fiber_cmd.append(f"-seq={mini_helix}")
    new_path: Path = output_path / "new.pdb"

    fiber_cmd.append(f"{str(new_path)}")

    try:
        print("Running fiber...", file=sys.stderr)
        run_fiber(fiber_cmd)
        logging.info(f"Ran fiber command: {' '.join(fiber_cmd)}")
        is_print_only = False
    except Exception as e:
        is_print_only = True
        logging.error(f"Error running fiber: {e}")

    if is_print_only:
        # On Windows or error just print the commands
        print("--- x3DNA Command ---", file=sys.stderr)
        print(" ".join(fiber_cmd), file=sys.stderr)

    return mini_helix, helix_orientation, is_print_only
