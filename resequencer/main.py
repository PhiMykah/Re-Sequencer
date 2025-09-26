import sys
from pathlib import Path

from biopandas.pdb import PandasPdb
from pandas import DataFrame, concat

from resequencer.parse import parse_args
from resequencer.residue import Residue
from resequencer.substitute import Substitution, load_substitution_file


def main(argv: list = sys.argv[1:]) -> None:
    # Parse command-line arguments
    args = parse_args(argv)

    # Collect pdb from file or from
    pdb_path = Path(args.pdb)
    if pdb_path.is_file():
        pdb: PandasPdb = PandasPdb().read_pdb(pdb_path)
    else:
        pdb: PandasPdb = PandasPdb().fetch_pdb(args.pdb)

    # Collect atoms from pdb
    atoms: DataFrame = pdb.df["ATOM"]

    # Import substitution file if it exists
    substitution_file = Path(args.input)
    if not substitution_file.is_file():
        raise Exception(f"Substitution file '{substitution_file}' does not exist!")

    substitutions: dict[int, Substitution] = load_substitution_file(substitution_file)

    residues: dict[int, Residue] = {}

    # Convert atoms DataFrame into a list of Residue objects, one per unique residue_number
    unique_res_nums = atoms["residue_number"].unique()
    for res_num in unique_res_nums:
        residues[res_num] = Residue(atoms, res_num)

    for res_num, sub in substitutions.items():
        print(
            f"Chain {sub.chain}: Substituting {sub.base} with {sub.new_base} on residue {res_num}",
            file=sys.stderr,
        )
        residues[res_num].substitute(sub.base, sub.new_base)

    # Combine all residues into a single dataframe
    print("Combining Residues...", file=sys.stderr)
    residue_dfs = [residue.to_dataframe() for residue in residues.values()]
    if residue_dfs:
        combined_residues = DataFrame(concat(residue_dfs, ignore_index=True))
    else:
        raise Exception("Unable to combine residues into single ATOM record!")

    pdb.df["ATOM"] = combined_residues
    print(
        f"Re-sequencing Complete! Outputting new pdb to {args.output}...",
        file=sys.stderr,
    )
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    pdb.to_pdb(path=args.output, records=None, gz=False)


if __name__ == "__main__":
    main()
