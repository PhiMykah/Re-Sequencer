import sys
from pathlib import Path

from biopandas.pdb import PandasPdb

from resequencer.parse import parse_args
from resequencer.substitute import pdb_substitution
from resequencer.addition import pdb_addition


def main(argv: list = sys.argv[1:]) -> None:
    # --------------------------------- Load PDB --------------------------------- #
    # Parse command-line arguments
    args = parse_args(argv)

    # Collect pdb from file or from
    pdb_path = Path(args.input)
    if pdb_path.is_file():
        pdb: PandasPdb = PandasPdb().read_pdb(pdb_path)
    else:
        pdb: PandasPdb = PandasPdb().fetch_pdb(args.input)

    # ----------------------- Calculate Tasks for Printing ----------------------- #

    if args.sub and args.add:
        operation_count = 2
    elif args.sub or args.add:
        operation_count = 1
    else:
        operation_count = 0

    operation_idx = 1

    # ------------------------------- Substitution ------------------------------- #

    if args.sub:
        pdb = pdb_substitution(args.sub, pdb)
        print(
            f"Substitution Complete! Operation ({operation_idx} of {operation_count}) ",
            file=sys.stderr,
        )
        operation_idx += 1

    # --------------------------------- Addition --------------------------------- #

    if args.add:
        pdb = pdb_addition(args.add, pdb, args.output_path)
        print(
            f"Sequence Addition Complete! Operation ({operation_idx} of {operation_count}) ",
            file=sys.stderr,
        )

    # ------------------------------- Output Result ------------------------------ #

    if operation_count != 0:
        output_path = Path(args.output_path)
        output_dir = output_path.parent if output_path.is_file() else output_path
        print(
            f'Outputting new pdb to path: "{output_dir}"',
            file=sys.stderr,
        )
        output_dir.mkdir(parents=True, exist_ok=True)
        pdb.to_pdb(
            path=str(output_dir.absolute() / "output.pdb"), records=None, gz=False
        )
    else:
        print("No operations performed. Closing...", file=sys.stderr)


if __name__ == "__main__":
    main()
