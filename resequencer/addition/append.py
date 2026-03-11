import typing
from pathlib import Path

from biopandas.pdb.pandas_pdb import PandasPdb

from resequencer.external import mini_helix_tail
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
):
    TARGET = ChainRole.TARGET.value
    OTHER = ChainRole.OTHER.value
    aligned_path = output_path / "aligned.pdb"
    aligned_pandas = PandasPdb().read_pdb(aligned_path)
    aligned = PDB(aligned_pandas.df, pdb.fasta)

    if len(aligned_ranges) != 2 or all(len(x) != 2 for x in aligned_ranges):
        raise ValueError(
            "Expected 2 tuples of 2 values for each chain in append_addition!"
        )

    for idx, aligned_range in enumerate(aligned_ranges):
        chain: Chain = aligned[aligned.get_chain_idx(addition.chains[idx])]  # type: ignore
        chain.remove(aligned_range[0], False)
        chain.remove(aligned_range[1], False)

    target_pdb_chain: Chain = pdb[pdb.get_chain_idx(addition.chains[TARGET])]  # type: ignore
    target_aligned_chain: Chain = aligned[
        aligned.get_chain_idx(addition.chains[TARGET])
    ]  # type: ignore
    target_len: int = mini_helix_tail()

    other_pdb_chain: Chain = pdb[pdb.get_chain_idx(addition.chains[OTHER])]  # type: ignore
    other_aligned_chain: Chain = aligned[aligned.get_chain_idx(addition.chains[OTHER])]  # type: ignore
    other_len: int = mini_helix_tail()

    if helix_orientation.lower() == "start":
        for i in range(target_len):
            target_aligned_chain.remove_at(-1)
        for i in range(other_len):
            other_aligned_chain.remove_at(0)
        pdb[addition.chains[TARGET]] = target_aligned_chain + target_pdb_chain
        pdb[addition.chains[OTHER]] = other_pdb_chain + other_aligned_chain
    else:
        if target_len > mini_helix_tail():
            for i in range(target_len):
                target_aligned_chain.remove_at(0)
        if other_len > mini_helix_tail():
            for i in range(other_len):
                other_aligned_chain.remove_at(-1)
        pdb[addition.chains[TARGET]] = target_pdb_chain + target_aligned_chain
        pdb[addition.chains[OTHER]] = other_aligned_chain + other_pdb_chain

    chain_len: int = len(pdb[addition.chains[TARGET]])
    pdb[addition.chains[TARGET]]._reorder_residue_numbers(1)  # type: ignore
    pdb[addition.chains[OTHER]]._reorder_residue_numbers(chain_len + 1)  # type: ignore
