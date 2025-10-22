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
    pdb_path = Path(args.pdb)
    if pdb_path.is_file():
        pdb: PandasPdb = PandasPdb().read_pdb(pdb_path)
    else:
        pdb: PandasPdb = PandasPdb().fetch_pdb(args.pdb)

    # ----------------------- Calculate Tasks for Printing ----------------------- #

    if args.input and args.add:
        operation_count = 2
    elif args.input or args.add:
        operation_count = 1
    else:
        operation_count = 0

    operation_idx = 1

    # ------------------------------- Substitution ------------------------------- #

    if args.input:
        pdb = pdb_substitution(args.input, pdb)
        print(
            f"Substitution Complete! Operation ({operation_idx} of {operation_count}) ",
            file=sys.stderr,
        )
        operation_idx += 1

    # --------------------------------- Addition --------------------------------- #

    if args.add:
        pdb = pdb_addition(args.add, pdb, args.output)
        print(
            f"Sequence Addition Complete! Operation ({operation_idx} of {operation_count}) ",
            file=sys.stderr,
        )

    # ------------------------------- Output Result ------------------------------ #

    if operation_count != 0:
        print(
            f"Outputting new pdb to {args.output}...",
            file=sys.stderr,
        )
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        pdb.to_pdb(path=args.output, records=None, gz=False)
    else:
        print("No operations performed. Closing...", file=sys.stderr)


if __name__ == "__main__":
    main()
