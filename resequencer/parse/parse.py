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
        help="Show this help message and exit",
    )
    # PDB File or PDB ID
    parser.add_argument(
        "-pdb",
        "--pdb",
        type=str,
        metavar="'PDB File or PDB ID'",
        dest="pdb",
        help="Path to PDB file or PDB ID",
    )
    # Input File (.in)
    parser.add_argument(
        "--input",
        "-in",
        type=str,
        metavar="'Substitution Input File (.in)'",
        dest="input",
        help="Input file for substitution",
    )
    # Output File (.pdb)
    parser.add_argument(
        "--output",
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
        pdb: str,
        sub_input: str,
        output: str,
    ) -> None:
        self.pdb: str = pdb
        self.input: str = sub_input
        self.output: str = output


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
    return CLargs(args.pdb, args.input, args.output)


if __name__ == "__main__":
    args: CLargs | None = parse_args()
    print(args)
