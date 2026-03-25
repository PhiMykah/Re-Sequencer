import logging
import sys
from pathlib import Path

from Bio.SeqRecord import SeqRecord
from biopandas.pdb.pandas_pdb import PandasPdb

from resequencer.addition import pdb_addition
from resequencer.io import get_input_fasta, get_input_pdb
from resequencer.parse import parse_args
from resequencer.pdb import PDB
from resequencer.substitute import pdb_substitution
from resequencer.walk import run_protein_walk


def entry(argv: list[str] = sys.argv[1:]) -> None:
    """entry-point function for command-line and script.

    Parameters
    ----------
    argv : list[str], optional
        Arguments to pass to Re-Sequencer, by default sys.argv[1:].
    """

    # --------------------------------- Load PDB --------------------------------- #
    # Parse command-line arguments
    cl_args = parse_args(argv)  # Re-Sequencer command-line arguments

    # Set logs level and output to stderr
    logging.basicConfig(
        level=cl_args.verbose_flag, stream=sys.stderr, format="%(message)s"
    )

    input: PandasPdb = get_input_pdb(cl_args.pdb_input)
    fasta: dict[str, SeqRecord] = get_input_fasta(cl_args.fasta_input)

    pdb = PDB(input.df, fasta)
    # ----------------------- Calculate Tasks for Printing ----------------------- #

    # Increment operation count based on number of tasks
    operation_count: int = (
        int(bool(cl_args.sub_input))
        + int(bool(cl_args.add_input))
        + int(bool(cl_args.walk_input))
    )
    operation_idx: int = 1

    # ------------------------------- Substitution ------------------------------- #
    if cl_args.sub_input:
        print("Performing Substitution...", file=sys.stderr)
        pdb_substitution(pdb, cl_args.sub_input)
        print(
            f"Substitution Complete! Operation ({operation_idx} of {operation_count}) ",
            file=sys.stderr,
        )
        operation_idx += 1

    # --------------------------------- Addition --------------------------------- #
    if cl_args.add_input:
        print("Performing Addition...", file=sys.stderr)
        pdb_addition(
            pdb,
            Path(cl_args.pdb_input),
            cl_args.add_input,
            cl_args.output_path,
            cl_args.nucleic_structure,
        )
        print(
            f"Sequence Addition Complete! Operation ({operation_idx} of {operation_count}) ",
            file=sys.stderr,
        )
        operation_idx += 1

    # ------------------------------- Protein Walk ------------------------------- #
    if cl_args.walk_input:
        print("Performing Protein Walk...", file=sys.stderr)
        run_protein_walk(
            pdb,
            Path(cl_args.pdb_input),
            Path(cl_args.walk_input),
            cl_args.output_path,
            cl_args.nucleic_structure,
        )
        print(
            f"Protein Walk Complete! Operation ({operation_idx} of {operation_count}) ",
            file=sys.stderr,
        )
        operation_idx += 1

    # ------------------------------- Output Result ------------------------------ #

    if operation_count != 0:
        output_path = cl_args.output_path
        output_file = cl_args.output_file
        print(
            f'Outputting new pdb to path: "{output_path}"',
            file=sys.stderr,
        )
        output_path.mkdir(parents=True, exist_ok=True)
        output: PandasPdb = pdb.output_pdb(input)
        output.to_pdb(path=(output_path / output_file), records=None, gz=False)
        print("Done!", file=sys.stderr)
    else:
        print("No operations performed. Closing...", file=sys.stderr)


if __name__ == "__main__":
    entry()
