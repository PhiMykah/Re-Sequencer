from dataclasses import dataclass
from pathlib import Path
from sys import stderr

from resequencer.pdb import PDB


@dataclass
class Substitution:
    """
    Represents a single amino acid or nucleotide substitution in a biomolecular sequence.

    Attributes
    ----------
    chain (str)
        The chain identifier where the substitution occurs.
    residue (int)
        The residue number in the sequence.
    base (str)
        The original nucleotide base  at the specified position.
    new_base (str)
        The new nucleotide base replacing the original.
    """

    chain: str
    residue: int
    base: str
    new_base: str


def load_substitution_file(file: Path) -> dict[int, Substitution]:
    """
    Loads a substitution file and returns a dictionary mapping residue numbers to Substitution objects.
    The substitution file should contain lines with four fields: chain, residue, base, and new_base,
    separated by either commas or whitespace.

    Lines starting with '#' or empty lines are ignored.
    Lines with more or fewer than four fields are also ignored.

    Parameters
    ----------
        file (Path)
            The path to the substitution file.
    Returns
    -------
        dict[int, Substitution]
            A dictionary where keys are residue numbers (int) and values are
            Substitution objects representing the substitutions to be made.
    """

    substitutions: dict[int, Substitution] = {}
    with file.open("r") as f:
        for line in f:
            line: str = line.strip()
            # Iterate through the loop skipping comments
            if not line or line.startswith("#"):
                continue

            # Split by comma separated values or by whitespace
            parts = (
                [p.strip() for p in line.split(",")] if "," in line else line.split()
            )

            # Ignore lines with items greater than 4
            if len(parts) != 4:
                continue
            chain, residue, base, new_base = parts

            substitutions[int(residue)] = Substitution(
                chain, int(residue), base, new_base
            )
    return substitutions


def pdb_substitution(pdb: PDB, sub_input: Path) -> None:

    if not sub_input.is_file():
        raise Exception(f"Substitution file '{sub_input}' does not exist!")

    substitutions: dict[int, Substitution] = load_substitution_file(sub_input)

    for res_num, sub in substitutions.items():
        print(
            f"Chain {sub.chain}: Substituting {sub.base} with {sub.new_base} on residue {res_num}",
            file=stderr,
        )
        pdb[sub.chain][res_num - 1].substitute(sub.base, sub.new_base)  # type: ignore
