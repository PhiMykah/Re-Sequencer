import logging
import re
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from resequencer.external import run_pymol, run_x3dna_fiber
from resequencer.pdb.pdb import PDB

# Precompiled regex for DNA/RNA bases
# Matches either: "DA", "DG", "DC", "DT" (DNA) OR "A", "G", "C", "U" (RNA)
BASE_PATTERN = re.compile(r"(?:D[AGCT]|[AGCTU])", re.IGNORECASE)


class ChainRole(Enum):
    """Enum to represent the role of a chain in an addition."""

    TARGET = 0
    OTHER = 1


@dataclass
class Addition:
    """
    Represents a single sequence addition in a biomolecular sequence.

    Attributes
    ----------
    chains: tuple[str, ...]
        The chain type of each nucleic acid strand, expects 2.
    old_seq: tuple[list[str], ...]
        Current sequences for both chains.
    new_seq: tuple[list[str], ...]
        New sequences for both chains.
    original_geometry: str
        Type of nucleic acid geometry to add.
    target_chain: str
        Name of target strain.
    """

    chains: tuple[str, ...]
    old_seq: tuple[list[str], ...]
    new_seq: tuple[list[str], ...]
    original_geometry: str
    target_chain: str

    @classmethod
    def load_addition_file(cls, file: Path, form: str):
        """
        Load and parse an addition file to create Addition objects.
        This class method reads a file containing chain addition specifications,
        where each pair of consecutive lines defines a single Addition object.
        Each line must contain exactly 3 whitespace-separated values: chain identifier,
        old sequence, and new sequence. Empty lines and comments (lines starting with '#')
        are ignored.

        Parameters
        ----------
        file: Path
            Path object pointing to the addition file to be read.
        form: str
            The original geometry form to assign to all created Addition objects.

        Returns
        -------
        list[Addition]
            A list of Addition objects created from pairs of lines in the file.
            Each Addition object contains:
            - chains: tuple of two chain identifiers
            - old_seq: tuple of two lists of parsed base sequences
            - new_seq: tuple of two lists of parsed base sequences
            - original_geometry: the form parameter
            - target_chain: the chain identifier from the first line of each pair

        Raises
        ------
        ValueError: If any line does not contain exactly 3 whitespace-separated values.

        Note
        ----
        Lines are processed in pairs. The chain from the first line of each pair
        is stored as the target_chain for the resulting Addition object.
        """
        with open(file, "r") as f:
            lines: list[str] = f.readlines()

        additions = []
        chains: list[str] = []
        old_seqs: list[list[str]] = []
        new_seqs: list[list[str]] = []
        target_chain: str = ""

        for idx, line in enumerate(lines):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            if len(line.split()) != 3:
                raise ValueError(
                    f"Error reading Line #{idx + 1}: Expected 3 values separated by whitespace, got {len(line.split())}"
                )

            chain, old_seq, new_seq = line.split(None, 2)

            chains.append(chain.lower())
            old_seqs.append(BASE_PATTERN.findall(old_seq))
            new_seqs.append(BASE_PATTERN.findall(new_seq))

            # Create Addition object every two lines
            if len(chains) == 1:
                target_chain = chain.lower()
            elif len(chains) == 2:
                additions.append(
                    cls(
                        chains=tuple(chains),
                        old_seq=tuple(old_seqs),
                        new_seq=tuple(new_seqs),
                        original_geometry=form,
                        target_chain=target_chain,
                    )
                )
                chains = []
                old_seqs = []
                new_seqs = []

        return additions


def pdb_addition(
    pdb: PDB, input_path: Path, add_input: Path, output_path: Path, form: str
):
    """
    Perform PDB addition operations by aligning and appending modification sequences.
    This function processes a list of additions from an input file, performs structural
    alignment using x3DNA and PyMOL, and appends the aligned structures to the original
    PDB file. For each addition, it generates a mini helix structure, determines helix
    orientation, performs alignment based on specified ranges, and integrates the result
    into the original PDB with reindexed atom numbering.

    Parameters
    ----------
    pdb : PDB
        The original PDB structure object to which additions will be appended.
    input_path : Path
        Path to the input PDB file.
    add_input : Path
        Path to the addition input file containing modification specifications.
    output_path : Path
        Directory Path where output files will be written.
    form : str
        Type of nucleic acid geometry to add

    Raises
    ------
    FileNotFoundError
        If the addition input file at add_input does not exist.
    """

    if not add_input.is_file():
        raise FileNotFoundError(f"Addition file '{add_input}' does not exist!")

    additions: list[Addition] = Addition.load_addition_file(add_input, form)

    _run_addition(pdb, additions, input_path, output_path)


def _run_addition(
    pdb: PDB, additions: list[Addition], input_path: Path, output_path: Path
):
    from .append import append_addition
    # ---------------------------- Create path objects --------------------------- #

    output_path = Path(output_path)
    output_dir = output_path.parent if output_path.is_file() else output_path
    Path.mkdir(output_dir, exist_ok=True, parents=True)

    # ------------------------- Iterate through Additions ------------------------ #

    for idx, addition in enumerate(additions):
        logging.info(f"Performing Addition ({idx + 1} of {len(additions)})")

        # ---------------------------------------------------------------------------- #
        #                               x3DNA Mini Helix                               #
        # ---------------------------------------------------------------------------- #

        mini_helix, helix_orientation, is_print_only = run_x3dna_fiber(
            addition, pdb, output_dir
        )

        # ---------------------------------------------------------------------------- #
        #                                     pymol                                    #
        # ---------------------------------------------------------------------------- #

        aligned_ranges = run_pymol(
            addition,
            pdb,
            helix_orientation,
            input_path,
            output_dir,
            is_print_only,
        )

        logging.info("Adding Aligned PDB to original PDB...")
        append_addition(pdb, addition, aligned_ranges, output_dir, helix_orientation)
        # Reorder all atom numbers for updated pdb
        pdb.reindex_atom_num(1)
