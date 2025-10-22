import sys
from pathlib import Path
import subprocess
import tempfile

from biopandas.pdb import PandasPdb
from pandas import DataFrame

from .add import Addition, load_addition_file


def pdb_addition(
    input_file: str, pdb: PandasPdb, output: Path | str = Path.cwd().resolve()
) -> PandasPdb:
    # Collect atoms from pdb
    atoms: DataFrame = pdb.df["ATOM"]

    # Import substitution file if it exists
    addition_file = Path(input_file)
    if not addition_file.is_file():
        raise Exception(f"Addition file '{addition_file}' does not exist!")

    additions: dict[int, Addition] = load_addition_file(addition_file)

    # Obtain the last n-10 bases and save it as new.pdb
    for idx, addition in additions.items():
        # Determine the target chain type to modify
        target_chain: int = (
            addition.target_chain - 1 if addition.target_chain > 0 else 0
        )
        chain_type: str = addition.chains[target_chain]

        chain: DataFrame = atoms[atoms["chain_id"] == chain_type.upper()]

        # ---------------------------------------------------------------------------- #
        #                               x3DNA Mini Helix                               #
        # ---------------------------------------------------------------------------- #
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

        # last_residue_numbers = [num for num, _ in last_residues]
        # last_residue_names = [
        #     name.strip().upper().strip("D") for _, name in last_residues
        # ]

        # if chain_type.upper() in ("A", "B"):
        #     print(f"x3dna_utils cp_std {chain_type.upper()}DNA", file=sys.stderr)

        is_rna = (
            "-rna"
            if ("U" in mini_helix.upper() and "T" not in mini_helix.upper())
            else ""
        )
        # Build commands
        fiber_cmd = ["fiber", f"-{chain_type.lower()}"]
        if is_rna:
            fiber_cmd.append(is_rna)
        fiber_cmd.append(f"-seq={mini_helix}")
        fiber_cmd.append("fiber.pdb")

        find_pair_cmd = ["find_pair", "fiber.pdb"]
        analyze_cmd = ["analyze"]

        mini_helix_path: Path = (
            output if isinstance(output, Path) else Path(output).parent.resolve()
        )
        Path.mkdir(mini_helix_path, exist_ok=True, parents=True)

        rebuild_cmd = None
        if chain_type.upper() == "A" or chain_type.upper() == "B":
            rebuild_cmd = [
                "rebuild",
                "-atomic",
                "bp_step.par",
                str(mini_helix_path / "new.pdb"),
            ]

        # ------------------------------ x3DNA Commands ------------------------------ #
        # Run on Linux/macOS, print on Windows
        if sys.platform.startswith("linux") or sys.platform == "darwin":
            with tempfile.TemporaryDirectory() as temp:
                try:
                    print("Running fiber...", file=sys.stderr)
                    subprocess.run(fiber_cmd, cwd=temp, check=True)
                    print("Running find_pair...", file=sys.stderr)
                    p = subprocess.run(
                        find_pair_cmd, cwd=temp, check=True, capture_output=True
                    )
                    print("Analyzing and rebuilding...", file=sys.stderr)
                    subprocess.run(analyze_cmd, input=p.stdout, cwd=temp, check=True)
                    if rebuild_cmd:
                        subprocess.run(rebuild_cmd, cwd=temp, check=True)
                    else:
                        print(f"Chain type: {chain_type.upper()}", file=sys.stderr)
                except subprocess.CalledProcessError as e:
                    print(f"Command failed ({e.returncode}): {e}", file=sys.stderr)
        else:
            # On Windows just print the commands
            print(" ".join(fiber_cmd), file=sys.stderr)
            print("find_pair new.pdb | analyze", file=sys.stderr)
            if rebuild_cmd:
                print(" ".join(rebuild_cmd), file=sys.stderr)
            else:
                print(f"Chain type: {chain_type.upper()}", file=sys.stderr)

        # ---------------------------------------------------------------------------- #
        #                                     pymol                                    #
        # ---------------------------------------------------------------------------- #

        
    return pdb
