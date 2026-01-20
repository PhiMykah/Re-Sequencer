import re
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Addition:
    """
    Represents a single sequence addition in a biomolecular sequence.

    Attributes
    ----------
    chains (tuple[str, ...])
        The chain type of each nucleic acid strand, expects 2
    target_chain (int)
        The Chain of where base pairs are added
    start_position (int)
        Position on chain after which to add bps

        (should be the end of chain to add)
    sequence (list[str])
        Sequence to add to original sequence,

        DA,DG,DC,DT: DNA bases

        A, G, C, U: RNA bases
    original_geometry (str)
        Geometry of bases to be added – A, B or Z

        (same as geometry of bases in the original pdb)
    total_bp (int)
        Total number of base pairs in original biomolecular sequence.
    """

    chains: tuple[str, ...]
    target_chain: str
    start_position: int
    sequence: list[str]
    original_geometry: str

    @classmethod
    def mini_helix_tail(cls):
        return 2


# Precompiled regex for DNA/RNA bases
# Matches either: "DA", "DG", "DC", "DT" (DNA) OR "A", "G", "C", "U" (RNA)
BASE_PATTERN = re.compile(r"(?:D[AGCT]|[AGCU])", re.IGNORECASE)


def load_addition_file(file: Path) -> dict[int, Addition]:
    additions: dict[int, Addition] = {}
    # TODO - Combine chain to add and start position in chain to add into residue 3'
    with file.open("r") as f:
        for line_number, line in enumerate(f, start=1):
            line: str = line.strip()
            # Iterate through the loop skipping comments
            if not line or line.startswith("#"):
                continue

            # Split by comma separated values or by whitespace
            parts = (
                [p.strip() for p in line.split(",")] if "," in line else line.split()
            )

            # Ignore lines with items greater than 6
            if len(parts) < 6:
                continue
            else:
                parts = parts[:6]

            *chains, target_chain, start_pos, seq, original_geometry = parts

            seq_list = BASE_PATTERN.findall(seq)

            reconstructed = "".join(seq_list)
            if reconstructed.upper() != seq.upper():
                invalid_chars = re.sub(BASE_PATTERN, "", seq, flags=re.IGNORECASE)
                raise ValueError(
                    f"Invalid base(s) '{invalid_chars}' found in sequence '{seq}' "
                    f"on line {line_number}"
                )

            additions[line_number] = Addition(
                tuple(chains),
                target_chain,
                int(start_pos),
                seq_list,
                original_geometry,
            )

    return additions
