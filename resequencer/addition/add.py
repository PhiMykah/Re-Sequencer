import re
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from sys import stderr

from resequencer.external import run_pymol, run_x3dna
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
        The chain type of each nucleic acid strand, expects 2
    old_seq: tuple[list[str], ...]
        Current sequences for both chains
    new_seq: tuple[list[str], ...]
        New sequences for both chains
    original_geometry: str
        Type of nucleic acid geometry to add
    target_chain: str
        Name of target strain
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
                    f"Error reading Line #{idx}: Expected 3 values separated by whitespace, got {len(line.split())}"
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
    from .append import append_addition

    if not add_input.is_file():
        raise Exception(f"Addition file '{add_input}' does not exist!")

    additions: list[Addition] = Addition.load_addition_file(add_input, form)

    # ---------------------------- Create path objects --------------------------- #

    output_path = Path(output_path)
    output_dir = output_path.parent if output_path.is_file() else output_path
    Path.mkdir(output_dir, exist_ok=True, parents=True)

    # ------------------------- Iterate through Additions ------------------------ #

    for idx, addition in enumerate(additions):
        print(f"Performing Addition ({idx + 1} of {len(additions)})", file=stderr)

        # ---------------------------------------------------------------------------- #
        #                               x3DNA Mini Helix                               #
        # ---------------------------------------------------------------------------- #

        mini_helix, helix_orientation, is_print_only = run_x3dna(
            addition, pdb, output_dir
        )

        # ---------------------------------------------------------------------------- #
        #                                     pymol                                    #
        # ---------------------------------------------------------------------------- #

        aligned_ranges = run_pymol(
            addition,
            pdb,
            mini_helix,
            helix_orientation,
            input_path,
            output_dir,
            is_print_only,
        )

        print("Adding Aligned PDB to original PDB...", file=stderr)
        append_addition(pdb, addition, aligned_ranges, output_dir, helix_orientation)
        pdb.reindex_atom_num(1)
