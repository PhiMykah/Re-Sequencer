from copy import deepcopy

import pytest

from resequencer.pdb import AtomRecord, Chain, Residue


@pytest.fixture
def atom_records():
    return [
        AtomRecord(
            record_name="test",
            atom_number=1,
            atom_name="N1",
            line_idx=1,
            residue_name="DA",
            residue_number=1,
            chain_id="A",
            x_coord=0.0,
            y_coord=1.0,
            z_coord=-1.0,
            occupancy=1.0,
            b_factor=1.0,
            segment_id="",
            element_symbol="A",
            charge=1.0,
        ),
        AtomRecord(
            record_name="test",
            atom_number=2,
            atom_name="C2",
            line_idx=2,
            residue_name="DA",
            residue_number=1,
            chain_id="A",
            x_coord=0.0,
            y_coord=1.0,
            z_coord=-1.0,
            occupancy=1.0,
            b_factor=1.0,
            segment_id="",
            element_symbol="A",
            charge=1.0,
        ),
        AtomRecord(
            record_name="test",
            atom_number=3,
            atom_name="N8",
            line_idx=3,
            residue_name="DA",
            residue_number=1,
            chain_id="A",
            x_coord=0.0,
            y_coord=1.0,
            z_coord=-1.0,
            occupancy=1.0,
            b_factor=1.0,
            segment_id="",
            element_symbol="A",
            charge=1.0,
        ),
    ]


@pytest.fixture
def sample_residues(atom_records):
    """Create sample residues for testing."""
    residues = []
    for i in range(3):
        res = Residue(
            atom_records,
            residue_number=i + 1,
            residue_name="DT",
            chain_id="A",
            starting_idx=i + 1,
        )
        residues.append(res)
    return residues


@pytest.fixture
def sample_chain(sample_residues):
    """Create a sample chain for testing."""
    return Chain(sample_residues, chain_id="A", starting_idx=1)


class TestChainInitialization:
    def test_chain_init_with_residues(self, sample_residues):
        chain = Chain(sample_residues, "A", 1)
        assert chain.chain_id == "A"
        assert len(chain) == 3
        assert chain.starting_idx == 1

    def test_chain_init_with_default_starting_idx(self, sample_residues):
        chain = Chain(sample_residues, "A")
        assert (
            chain.starting_idx == -1
            or chain.starting_idx == sample_residues[0].starting_idx
        )

    def test_chain_init_empty_residues(self):
        chain = Chain([], "B", 0)
        assert chain.chain_id == "B"
        assert len(chain) == 0


class TestChainProperties:
    def test_residues_getter(self, sample_chain):
        assert isinstance(sample_chain.residues, list)
        assert len(sample_chain.residues) == 3

    def test_chain_id_getter(self, sample_chain):
        assert sample_chain.chain_id == "A"

    def test_chain_id_setter(self, sample_chain):
        sample_chain.chain_id = "B"
        assert sample_chain.chain_id == "B"
        assert all(res.chain_id == "B" for res in sample_chain.residues)

    def test_starting_idx_setter(self, sample_chain):
        sample_chain.starting_idx = 10
        assert sample_chain.starting_idx == 10

    def test_starting_residue_getter(self, sample_chain):
        assert isinstance(sample_chain.starting_residue, int)

    def test_starting_residue_setter(self, sample_chain):
        sample_chain.starting_residue = 100
        assert sample_chain.starting_residue == 100


class TestChainModification:
    def test_append_residue(self, sample_chain):
        new_res = Residue(
            [], residue_number=4, residue_name="GLY", chain_id="A", starting_idx=4
        )
        initial_len = len(sample_chain)
        sample_chain.append(new_res)
        assert len(sample_chain) == initial_len + 1
        assert new_res.chain_id == "A"

    def test_insert_residue(self, sample_chain):
        new_res = Residue(
            [], residue_number=2, residue_name="GLY", chain_id="A", starting_idx=2
        )
        sample_chain.insert(new_res, 1)
        assert len(sample_chain) == 4
        assert sample_chain[1] == new_res

    def test_remove_by_number(self, sample_chain):
        initial_len = len(sample_chain)
        sample_chain.remove(2)
        assert len(sample_chain) == initial_len - 1

    def test_remove_invalid_number(self, sample_chain):
        with pytest.raises(ValueError):
            sample_chain.remove(999)

    def test_remove_at_index(self, sample_chain):
        initial_len = len(sample_chain)
        sample_chain.remove_at(0)
        assert len(sample_chain) == initial_len - 1

    def test_remove_at_no_reorder(self, sample_chain):
        sample_chain.remove_at(0, reorder=False)
        assert len(sample_chain) == 2


class TestChainCharacteristics:
    def test_total_length(self, sample_chain):
        total = sample_chain.total_length()
        assert isinstance(total, int)
        assert total > 0

    def test_chain_length(self, sample_chain):
        assert len(sample_chain) == 3


class TestChainMagicMethods:
    def test_deepcopy(self, sample_chain):
        chain_copy = deepcopy(sample_chain)
        assert len(chain_copy) == len(sample_chain)
        assert chain_copy is not sample_chain
        assert chain_copy.chain_id == sample_chain.chain_id

    def test_add_chains(self, sample_residues):
        chain1 = Chain(sample_residues[:2], "A", 1)
        chain2 = Chain(sample_residues[2:], "A", 3)
        combined = chain1 + chain2
        assert len(combined) == 3
        assert combined.chain_id == "A"

    def test_getitem_single_index(self, sample_chain):
        res = sample_chain[0]
        assert isinstance(res, Residue)

    def test_getitem_slice(self, sample_chain):
        residues = sample_chain[0:2]
        assert isinstance(residues, list)
        assert len(residues) == 2

    def test_iterator(self, sample_chain):
        count = 0
        for residue in sample_chain:
            assert isinstance(residue, Residue)
            count += 1
        assert count == 3
