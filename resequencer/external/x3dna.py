import subprocess
import sys
import tempfile
from pathlib import Path

from pandas import DataFrame

from resequencer.addition import Addition


def run_x3dna(
    addition: Addition,
    chain: DataFrame,
    chain_type: str,
    new_path: Path,
):
    # Calculate the length of the first mini helix part
    base_count: int = addition.total_bp - 10  # + len(addition.sequence)

    res_cols = ["residue_number", "residue_name"]
    if not set(res_cols).issubset(chain.columns):
        raise KeyError(f"Chain DataFrame missing required columns: {res_cols}")

    unique_residues = chain[res_cols].drop_duplicates().reset_index(drop=True)
    # optional: as a list of tuples (residue_number, residue_name)
    unique_residue_list = [tuple(x) for x in unique_residues.to_numpy()]

    if base_count <= 0:
        raise ValueError("base_count must be > 0")
    if base_count > len(unique_residue_list):
        raise ValueError(
            f"Requested base_count ({base_count}) exceeds available residues ({len(unique_residue_list)})"
        )

    # take the last `base_count` residues (preserves their original order)
    last_residues = unique_residue_list[-base_count:]

    # sequence made from residue names (trimmed and uppercased)
    tail_sequence = "".join(
        name.strip().upper().strip("D") for _, name in last_residues
    )
    new_sequence = "".join(a.strip().upper().strip("D") for a in addition.sequence)
    mini_helix: str = tail_sequence + new_sequence

    is_rna = (
        "-rna" if ("U" in mini_helix.upper() and "T" not in mini_helix.upper()) else ""
    )
    # Build command
    fiber_cmd = ["fiber", f"-{chain_type.lower()}"]
    if is_rna:
        fiber_cmd.append(is_rna)
    fiber_cmd.append(f"-seq={mini_helix}")
    fiber_cmd.append(str(new_path))

    is_print_only = False
    # ------------------------------ x3DNA Commands ------------------------------ #
    # Run on Linux/macOS, print on Windows
    if sys.platform.startswith("linux") or sys.platform == "darwin":
        with tempfile.TemporaryDirectory() as temp:
            try:
                print("Running fiber...", file=sys.stderr)
                subprocess.run(fiber_cmd, cwd=temp, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Command failed ({e.returncode}): {e}", file=sys.stderr)
                is_print_only = True
    else:
        is_print_only = True

    if is_print_only:
        # On Windows or error just print the commands
        print("--- x3DNA Command ---", file=sys.stderr)
        print(" ".join(fiber_cmd), file=sys.stderr)

    return is_print_only
