import sys
from pathlib import Path
import subprocess
import tempfile

from biopandas.pdb import PandasPdb
from pandas import DataFrame

# conda install -c conda-forge -c schrodinger pymol-bundle, requires python 3.10
from pymol import cmd

from .add import Addition, load_addition_file


def pdb_addition(
    input_file: str,
    pdb: PandasPdb,
    output: Path | str = Path.cwd().resolve(),
) -> PandasPdb:
    # Collect atoms from pdb
    atoms: DataFrame = pdb.df["ATOM"]

    # Import substitution file if it exists
    addition_file = Path(input_file)
    if not addition_file.is_file():
        raise Exception(f"Addition file '{addition_file}' does not exist!")

    additions: dict[int, Addition] = load_addition_file(addition_file)

    # ---------------------------- Create path objects --------------------------- #

    output_path: Path = (
        output if isinstance(output, Path) else Path(output).parent.resolve()
    )
    Path.mkdir(output_path, exist_ok=True, parents=True)
    new_path: Path = output_path / "new.pdb"
    aligned_path: Path = output_path / "aligned.pdb"

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

        is_print_only = run_x3dna(addition, chain, chain_type, new_path)

        # ---------------------------------------------------------------------------- #
        #                                     pymol                                    #
        # ---------------------------------------------------------------------------- #

        # TODO:
        # Questions:
        # 1. Best way to calculate resi range?
        run_pymol(
            addition, chain_type, pdb.pdb_path, new_path, aligned_path, is_print_only
        )
    return pdb


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


def run_pymol(
    addition: Addition,
    chain_type: str,
    input_file: str,
    new_path: Path,
    aligned_path: Path,
    print_mode: bool = False,
) -> None:
    # ----------------------------- Calculate Extract ---------------------------- #
    excess_length: int = addition.total_bp - 10

    # NOTE: Determine if this is the proper method to calculate the resi range
    """
    Example: If chain type is A, total_bp is 12, start position is 12
        - excess_length is 12-10 = 2
        - excess_start is 12 - 2 = 10.
        - We then add 1 to not include excess start, or 10 + 1 = 11
        - Lastly, it ends at start position = 12.
    Example: If chain type is B, total_bp is 12, start position is 12,
        - excess length is 12-10 = 2, 
        - excess_start = 12
        - We then add 1 to not enclude start pos, or 12 + 1 = 13
        - excess_end = it ends at start position + excess_length = 12 + 2 = 14
    """

    excess_start: int = (
        addition.start_position
        if chain_type.upper() == "B"
        else addition.start_position - excess_length
    )
    excess_start += 1
    excess_end: int = (
        addition.start_position + excess_length
        if chain_type.upper() == "B"
        else addition.start_position
    )
    extraction: str = " ".join(
        [f"chain {chain_type.upper()}", "and", "resi", f"{excess_start}-{excess_end}"]
    )

    # --------------------------------- Run pymol -------------------------------- #
    if not print_mode:
        print("Running pymol...", file=sys.stderr)
        # load /PATH/TO/input_file.pdb
        cmd.load(input_file)
        original_obj: str = cmd.get_names("objects")[0]
        # extract temp, chain A and resi 11-12 or chain B and resi 13-14)
        cmd.extract("temp", extraction)
        # load /PATH/TO/new.pdb
        cmd.load(str(new_path))
        # delete input_file
        cmd.delete(original_obj)
        # super new.pdb, temp
        cmd.super("new", "temp")
        # multisave /PATH/TO/aligned.pdb
        cmd.multisave(str(aligned_path))
    else:
        print_output = []
        original_obj = Path(input_file).stem
        print_output.extend(["pymol", input_file, "-c", "-d"])
        command = []
        command.extend(
            [
                f"extract temp, {extraction}",
                f"load {str(new_path)}",
                f"delete {original_obj}",
                "super new, temp",
                f"multisave {str(aligned_path)}",
            ]
        )
        print_output.append(f"'{';'.join(command)}'")
        print("--- pymol Commands ---", file=sys.stderr)
        print(" ".join(print_output), file=sys.stderr)
