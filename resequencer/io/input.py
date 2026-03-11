from pathlib import Path
from sys import stderr

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from biopandas.pdb import PandasPdb


def get_input_pdb(target_pdb: Path | str) -> PandasPdb:
    """
    Create a PandasPdb Object based on given path or PDB code from Protein Databank.

     Assumes PDB codes are 4-characters in length.

    Parameters
    ----------
    target_pdb : Path | str
        Input pdb as file path or 4-char code

    Returns
    -------
    PandasPdb
        PandasPdb object representing inputted PDB.

    Raises
    ------
    FileExistsError
        Occurs when attempting to convert input file into pdb suffix
    FileNotFoundError
        Occurs when input path does not exist
    ValueError
        Occurs when failing to fetch pdb code
    """
    # Cast string to path if it does not fit code structure
    if isinstance(target_pdb, str) and len(target_pdb) != 4:
        target_pdb = Path(target_pdb)

    # Test if pdb input is path or code
    if isinstance(target_pdb, Path):
        print("Attempting to read input PDB...", file=stderr)
        if target_pdb.is_file():
            # Try to read file as a pdb file
            try:
                input: PandasPdb = PandasPdb().read_pdb(
                    target_pdb.rename(target_pdb.with_suffix(".pdb"))
                )
            except FileExistsError:
                # Raise error if renaming with pdb suffix causes file clash
                raise FileExistsError(
                    "Unable to read input file as pdb, a pdb file already exists!"
                )
        else:
            raise FileNotFoundError(f"Unable to find input PDB file: {target_pdb}")
    else:
        print(f"Attempting to fetch pdb with code {target_pdb}", file=stderr)
        try:
            input: PandasPdb = PandasPdb().fetch_pdb(target_pdb, source="pdb")
        except Exception:
            raise ValueError(
                f"Could not find pdb with code: {target_pdb}. "
                f"Either pdb does not exist or could not reach website."
            )

    print("Successfully loaded PDB file!", file=stderr)
    return input


def get_input_fasta(target_fasta: Path) -> dict[str, SeqRecord]:
    """
    Create a fasta SeqRecord dictionary based on given path.

    Parameters
    ----------
    target_fasta : Path
        Input fasta file

    Returns
    -------
    dict[str, SeqRecord]
        Dictionary matching record name to each sequence

    Raises
    ------
    FileNotFoundError
        Occurs when input path does not exist
    ValueError
        Occurs when failing to read fasta file
    """
    print("Attempting to read fasta file...", file=stderr)
    if target_fasta.is_file():
        try:
            input: dict[str, SeqRecord] = SeqIO.to_dict(
                SeqIO.parse(target_fasta, "fasta")
            )
        except Exception:
            raise ValueError(f"Unable to read fasta pdb: {target_fasta}")
    else:
        raise FileNotFoundError(f"Unable to find input fasta file: {target_fasta}")

    print("Successfully loaded fasta file!", file=stderr)
    return input
