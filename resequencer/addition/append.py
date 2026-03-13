import typing
from pathlib import Path

from biopandas.pdb.pandas_pdb import PandasPdb

from resequencer.external import MINI_HELIX_TAIL
from resequencer.pdb import PDB, Chain

from .add import ChainRole

if typing.TYPE_CHECKING:
    from .add import Addition


def append_addition(
    pdb: PDB,
    addition: "Addition",
    aligned_ranges: list[tuple],
    output_path: Path,
    helix_orientation: str,
) -> None:
    """
    Append aligned chain segments to PDB chains based on helix orientation.
    This function modifies a PDB structure by removing specified residue ranges
    from aligned chains and concatenating them with the original PDB chains in
    an order determined by the helix orientation. Residue numbers are then
    renumbered sequentially for both chains.

    Parameters
    ----------
    pdb : PDB
        The PDB structure to be modified with appended chain segments.
    addition : Addition
        An Addition object containing information about the chains to be modified,
        including target and other chain identifiers.
    aligned_ranges : list[tuple]
        A list of two tuples, each containing a pair of integers representing
        the start and end residue positions to remove from aligned chains.
    output_path : Path
        The file path where the aligned PDB file is located.
    helix_orientation : str
        The orientation of the helix ("start" or "end"). Determines the order
        in which chains are concatenated and which residues are removed.

    Raises
    ------
    ValueError
        If aligned_ranges does not contain exactly 2 tuples with 2 values each.
    """
    TARGET: int = ChainRole.TARGET.value
    OTHER: int = ChainRole.OTHER.value

    # Obtain aligned pdb
    aligned_path = output_path / "aligned.pdb"
    aligned_pandas = PandasPdb().read_pdb(aligned_path)
    aligned = PDB(aligned_pandas.df, pdb.fasta)

    if len(aligned_ranges) != 2 or all(len(x) != 2 for x in aligned_ranges):
        raise ValueError(
            "Expected 2 tuples of 2 values for each chain in append_addition!"
        )

    # Remove mini_helix additions from aligned pdb
    for idx, aligned_range in enumerate(aligned_ranges):
        chain: Chain = aligned[aligned.get_chain_idx(addition.chains[idx])]  # type: ignore
        chain.remove(aligned_range[0], False)
        chain.remove(aligned_range[1], False)

    # collect target chain and its aligned addition
    target_pdb_chain: Chain = pdb[pdb.get_chain_idx(addition.chains[TARGET])]  # type: ignore
    target_aligned_chain: Chain = aligned[
        aligned.get_chain_idx(addition.chains[TARGET])
    ]  # type: ignore
    target_len: int = MINI_HELIX_TAIL

    # collect other chain and its aligned addition
    other_pdb_chain: Chain = pdb[pdb.get_chain_idx(addition.chains[OTHER])]  # type: ignore
    other_aligned_chain: Chain = aligned[aligned.get_chain_idx(addition.chains[OTHER])]  # type: ignore
    other_len: int = MINI_HELIX_TAIL

    if helix_orientation.lower() == "start":
        # Remove the additional overlap from aligned chain
        if target_len < len(target_aligned_chain):
            for _ in range(target_len):
                pdb[addition.chains[TARGET]].remove_at(-1)  # type: ignore
        if other_len < len(other_aligned_chain):
            for _ in range(other_len):
                pdb[addition.chains[OTHER]].remove_at(0)  # type: ignore

        # Add updated chains to original pdb
        pdb[addition.chains[TARGET]] = target_aligned_chain + target_pdb_chain
        pdb[addition.chains[OTHER]] = other_pdb_chain + other_aligned_chain
    else:
        # Remove the additional overlap from aligned chain
        if target_len < len(target_aligned_chain):
            for _ in range(target_len):
                pdb[addition.chains[TARGET]].remove_at(0)  # type: ignore
        if other_len < len(other_aligned_chain):
            for _ in range(other_len):
                pdb[addition.chains[OTHER]].remove_at(-1)  # type: ignore

        # Add updated chains to original pdb
        pdb[addition.chains[TARGET]] = target_pdb_chain + target_aligned_chain
        pdb[addition.chains[OTHER]] = other_aligned_chain + other_pdb_chain

    chain_len: int = len(pdb[addition.chains[TARGET]])

    # Reorder residue numbers to original
    pdb[addition.chains[TARGET]]._reorder_residue_numbers(1)  # type: ignore
    pdb[addition.chains[OTHER]]._reorder_residue_numbers(chain_len + 1)  # type: ignore
