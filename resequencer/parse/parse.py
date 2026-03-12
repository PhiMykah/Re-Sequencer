import argparse
import logging
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    """
    Creates and configures an ArgumentParser for the Re-Sequencer command line interface.
    An ArgumentParser object configured with the following arguments:

     --help: Shows the help message and exits.
     --input: Path to PDB file or PDB ID.
     --fasta: Path to the fasta file.
     --form: Desired Structure of DNA/RNA, eg. b-form, a-form, z-form.
     --sub: Input file for substitution.
     --add: Path to chain addition input file.
     --walk: Path to protein walk input file.
     --output: Output Path for Finished PDB and Intermediate Steps.
     --file: File name for final output.
     --verbose: Verbose output printing.
     
    Returns
    -------
    argparse.ArgumentParser
        ArgumentParser object to containing incoming arguments.
    """
    parser = argparse.ArgumentParser(
        prog="Re-Sequencer", description="Re-Sequencer Command Line Arguments"
    )

    # Help Message
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
        type=str,
        metavar="'File or PDB ID'",
        required=True,
        dest="pdb_input",
        help="Path to PDB file or PDB ID.",
    )
    # Fasta File
    parser.add_argument(
        "--fasta",
        "-fasta",
        "-f",
        type=Path,
        metavar="'file.fasta'",
        required=True,
        dest="fasta_input",
        help="Path to the fasta file.",
    )
    # Nucleic Acid Form
    parser.add_argument(
        "--form",
        "-form",
        type=str,
        choices=["a", "b", "z"],
        dest="nucleic_structure",
        help="Desired Structure of DNA/RNA, eg. b-form, a-form, z-form.",
        default="",
    )
    # Substitution File (.in)
    parser.add_argument(
        "--sub",
        "-sub",
        type=str,
        metavar="'substitute.in'",
        dest="sub_input",
        help="Input file for substitution.",
        default="",
    )
    # Chain Addition (.chain)
    parser.add_argument(
        "--add",
        "-add",
        type=str,
        metavar="'add.chain'",
        dest="add_input",
        help="Path to chain addition input file.",
        default="",
    )
    # Protein Walk (.walk)
    parser.add_argument(
        "--walk",
        "-walk",
        type=str,
        metavar="'file.walk'",
        dest="walk_input",
        help="Path to protein walk input file.",
        default="",
    )
    # Output Path
    parser.add_argument(
        "--output",
        "-output",
        type=Path,
        metavar="'Path'",
        dest="output_path",
        help="Output Path for Finished PDB and Intermediate Steps.",
        default=Path("output"),
    )
    # Output File (.pdb)
    parser.add_argument(
        "--file",
        "-file",
        type=str,
        metavar="'Output.pdb'",
        dest="output_file",
        help="File name for final output.",
        default="output.pdb",
    )
    # Verbose Flag
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_const",
        dest="verbose_flag",
        help="Verbose output printing.",
        const=logging.INFO,
    )
    return parser


class CLargs(argparse.Namespace):
    """
    A namespace class for storing command-line arguments.

    Attributes
    ----------
    pdb_input : Path | str
        Path to the PDB file.
    fasta_input : Path
        Path to the fasta file.
    nucleic_structure: str
        Desired Structure of DNA/RNA, eg. b-form, a-form, z-form;
        Empty string if not adding.
    sub_input : Path | None
        Path to substitution input file.
    add_input : Path | None
        Path to chain addition input file.
    walk_input : Path | None
        Path to protein walk input file.
    output_path : Path
        Output directory path.
    output_file : str
        File name for final output.
    verbose_flag : bool
        Whether or not to print verbose data.
    """

    def __init__(
        self,
        pdb_input: str,
        fasta_input: Path,
        nucleic_structure: str,
        sub_input: str,
        add_input: str,
        walk_input: str,
        output_path: Path,
        output_file: str,
        verbose_flag: bool,
    ) -> None:
        self.pdb_input: Path | str = (
            Path(pdb_input) if pdb_input.endswith("pdb") else pdb_input
        )
        self.fasta_input: Path = fasta_input
        self.nucleic_structure: str = nucleic_structure
        self.sub_input: Path | None = Path(sub_input) if sub_input != "" else None
        self.add_input: Path | None = Path(add_input) if add_input != "" else None
        self.walk_input: Path | None = Path(walk_input) if walk_input != "" else None
        self.output_path: Path = output_path
        self.output_file: str = output_file
        self.verbose_flag: bool = verbose_flag


def parse_args(argv: list[str]) -> CLargs:
    """Parse arguments from provided list and return CLargs object.

    Parameters
    ----------
    argv : list[str]
        List of command-line arguments.

    Returns
    -------
    CLargs
        Namespace Object with Resequencer attributes.
    """
    parser: argparse.ArgumentParser = build_parser()
    args: argparse.Namespace = parser.parse_args(argv)

    if args.add_input and args.nucleic_structure == "":
        parser.error("--form is required when --add is specified")

    if args.walk_input and args.nucleic_structure == "":
        parser.error("--form is required when --walk is specified")

    if (
        isinstance(args.pdb_input, str)
        and not (args.pdb_input.endswith("pdb"))
        and len(args.pdb_input) != 4
    ):
        parser.error("--input pdb code from Protein Databank must be 4 characters")

    return CLargs(
        args.pdb_input,
        args.fasta_input,
        args.nucleic_structure,
        args.sub_input,
        args.add_input,
        args.walk_input,
        args.output_path,
        args.output_file,
        args.verbose_flag,
    )


if __name__ == "__main__":
    import sys

    # Print arguments from command-line for testing
    print(parse_args(sys.argv[1:]))
