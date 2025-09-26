from dataclasses import dataclass
from pathlib import Path


@dataclass
class Substitution:
    """
    Represents a single amino acid or nucleotide substitution in a biomolecular sequence.

    Attributes
    ----------
        chain (str)
            The chain identifier where the substitution occurs.
        residue (int)
            he residue number in the sequence.
        base (str)
            The original base or amino acid at the specified position.
        new_base (str)
            The new base or amino acid replacing the original.
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
            if "," in line:
                parts = [p.strip() for p in line.split(",")]
            else:
                parts = line.split()

            # Ignore lines with items greater than 4
            if len(parts) != 4:
                continue
            chain, residue, base, new_base = parts

            substitutions[int(residue)] = Substitution(
                chain, int(residue), base, new_base
            )
    return substitutions
