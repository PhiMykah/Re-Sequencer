import logging
import sys
from pathlib import Path
import os
from resequencer.external.x3dna.find_pair_bindings import run_find_pair


def run_x3dna_pair(input_path: Path, output_path: Path):
    is_print_only = False

    # Prepare Find Pair
    curr_dir: str = os.getcwd()
    find_pair_path = output_path / "find_pair"
    Path.mkdir(find_pair_path, exist_ok=True, parents=True)
    full_input_path: Path = Path.absolute(input_path)
    full_output_path: Path = Path.absolute(output_path)

    output_file: Path = full_output_path / f"{input_path.name}.par"
    pair_cmd: list[str] = [
        "find_pair",
        "-a",
        str(full_input_path),
        str(output_file),
    ]
    try:
        print("Running find_pair...", file=sys.stderr)
        os.chdir(find_pair_path)
        run_find_pair()
        os.chdir(curr_dir)
        logging.info(f"Ran find_pair command: {' '.join(pair_cmd)}")
    except Exception as e:
        is_print_only = True
        logging.error(f"Error running fiber: {e}")

    if is_print_only:
        # On Windows or error just print the commands
        print("--- x3DNA Command ---", file=sys.stderr)
        print(" ".join(pair_cmd), file=sys.stderr)

    return output_file, is_print_only
