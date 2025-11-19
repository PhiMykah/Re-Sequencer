import argparse


def build_parser() -> argparse.ArgumentParser:
    """
    Creates and configures an ArgumentParser for the Re-Sequencer command line interface.
    An ArgumentParser object configured with the following arguments:

     -help: Shows the help message and exits.
     -pdb, --pdb: Path to a PDB file or a PDB ID (string).
     --input, -in: Path to the substitution input file (.in) (string).
     --output, -out: Path to the output PDB file after substitution (string, default: 'output.pdb').

    Returns
    -------
    argparse.ArgumentParser
        ArgumentParser object to parse incoming arguments
    """
    parser = argparse.ArgumentParser(description="Re-Sequencer Command Line Arguments")

    # Help
    parser.add_argument(
        "-help",
        action="help",
        default=argparse.SUPPRESS,
        help="",
    )
    # PDB File or PDB ID
    parser.add_argument(
        "--input",
        "-input",
        "-in",
        type=str,
        metavar="'PDB File or PDB ID'",
        required=True,
        dest="input",
        help="Path to PDB file or PDB ID",
    )
    # Input File (.in)
    parser.add_argument(
        "--sub",
        "-substitute",
        "-sub",
        type=str,
        metavar="'Substitution Input File (.in)'",
        dest="sub",
        help="Input file for substitution",
        default="",
    )
    # Chain Addition (.chain)
    parser.add_argument(
        "--add",
        "-add",
        "-a",
        "--add-chain",
        type=str,
        dest="add",
        help="Input file for chain additions",
        default="",
    )
    # Output File (.pdb)
    parser.add_argument(
        "--output",
        "-output",
        "-out",
        type=str,
        metavar="'Output File (.pdb)",
        dest="output",
        help="Output PDB File after substitution",
        default="output.pdb",
    )
    return parser


class CLargs(argparse.Namespace):
    """
    A namespace class for storing command-line arguments.

    Attributes
    ----------
    pdb (str)
        Path to the PDB file.
    input (str)
        Path to the input file (provided as 'sub_input' in the constructor).
    output (str)
        Path to the output file.
    """

    def __init__(
        self,
        pdb_input: str,
        sub: str,
        output: str,
        add: str,
    ) -> None:
        self.input: str = pdb_input
        self.sub: str = sub
        self.output: str = output
        self.add: str = add


def parse_args(argv: list[str] | None = None) -> CLargs:
    """
    Parses command-line arguments and returns them as a CLargs object.

    Parameters
    ----------
        argv (list[str] | None, optional)
            List of command-line arguments to parse.
            If None, parses arguments from sys.argv.
    Returns
    -------
        CLargs
            An object containing the parsed command-line arguments (pdb, input, output).
    """

    parser: argparse.ArgumentParser = build_parser()
    args: argparse.Namespace = parser.parse_args(argv)
    return CLargs(args.input, args.sub, args.output, args.add)


if __name__ == "__main__":
    args: CLargs | None = parse_args()
    print(args)
