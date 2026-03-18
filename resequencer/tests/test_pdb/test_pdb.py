from unittest.mock import Mock

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

from resequencer.pdb import PDB, Chain


class TestPDBInitialization:
    def test_init_with_valid_atom_dataframe(self):
        atom_df = DataFrame(
            {
                "chain_id": ["A", "A"],
                "residue_number": [1, 1],
                "record_name": ["ATOM", "ATOM"],
                "atom_number": [1, 2],
                "atom_name": ["N", "CA"],
                "residue_name": ["ALA", "ALA"],
                "x_coord": [0.0, 1.0],
                "y_coord": [0.0, 1.0],
                "z_coord": [0.0, 1.0],
                "occupancy": [1.0, 1.0],
                "b_factor": [0.0, 0.0],
                "segment_id": ["", ""],
                "element_symbol": ["N", "C"],
                "charge": [0, 0],
            }
        )
        pdb_data = {
            "ATOM": atom_df,
            "HETATM": DataFrame(),
            "ANISOU": DataFrame(),
            "OTHERS": DataFrame(),
        }
        pdb = PDB(pdb_data)
        assert len(pdb.chains) == 1
        assert pdb.chains[0].chain_id == "A"

    def test_init_with_empty_atom_dataframe_raises_error(self):
        pdb_data = {
            "ATOM": DataFrame(),
            "HETATM": DataFrame(),
            "ANISOU": DataFrame(),
            "OTHERS": DataFrame(),
        }
        with pytest.raises(ValueError, match="Atom Dataframe for PDB is Empty!"):
            PDB(pdb_data)

    def test_init_with_fasta(self):
        atom_df = DataFrame(
            {
                "chain_id": ["A"],
                "residue_number": [1],
                "record_name": ["ATOM"],
                "atom_number": [1],
                "atom_name": ["N"],
                "residue_name": ["ALA"],
                "x_coord": [0.0],
                "y_coord": [0.0],
                "z_coord": [0.0],
                "occupancy": [1.0],
                "b_factor": [0.0],
                "segment_id": [""],
                "element_symbol": ["N"],
                "charge": [0],
            }
        )
        fasta = {"seq1": SeqRecord(seq=Seq("MVLWAALLVTFLAGCQAKVEQAESMGQVD"), id="seq1")}
        pdb_data = {"ATOM": atom_df}
        pdb = PDB(pdb_data, fasta)
        assert pdb.fasta == fasta


class TestPDBProperties:
    @pytest.fixture
    def pdb_instance(self):
        atom_df = DataFrame(
            {
                "chain_id": ["A"],
                "residue_number": [1],
                "record_name": ["ATOM"],
                "atom_number": [1],
                "atom_name": ["N"],
                "residue_name": ["ALA"],
                "x_coord": [0.0],
                "y_coord": [0.0],
                "z_coord": [0.0],
                "occupancy": [1.0],
                "b_factor": [0.0],
                "segment_id": [""],
                "element_symbol": ["N"],
                "charge": [0],
            }
        )
        return PDB({"ATOM": atom_df})

    def test_chains_getter(self, pdb_instance):
        assert isinstance(pdb_instance.chains, list)
        assert len(pdb_instance.chains) > 0

    def test_chains_setter(self, pdb_instance):
        mock_chains = [Mock(spec=Chain)]
        pdb_instance.chains = mock_chains
        assert pdb_instance.chains == mock_chains

    def test_fasta_getter(self, pdb_instance):
        assert isinstance(pdb_instance.fasta, dict)

    def test_fasta_setter(self, pdb_instance):
        fasta = {"seq1": SeqRecord(seq=Seq("ACGT"), id="seq1")}
        pdb_instance.fasta = fasta
        assert pdb_instance.fasta == fasta


class TestPDBDataAccess:
    @pytest.fixture
    def pdb_instance(self):
        atom_df = DataFrame(
            {
                "chain_id": ["A", "B"],
                "residue_number": [1, 1],
                "record_name": ["ATOM", "ATOM"],
                "atom_number": [1, 2],
                "atom_name": ["N", "N"],
                "residue_name": ["ALA", "GLY"],
                "x_coord": [0.0, 1.0],
                "y_coord": [0.0, 1.0],
                "z_coord": [0.0, 1.0],
                "occupancy": [1.0, 1.0],
                "b_factor": [0.0, 0.0],
                "segment_id": ["", ""],
                "element_symbol": ["N", "N"],
                "charge": [0, 0],
            }
        )
        return PDB({"ATOM": atom_df})

    def test_get_chain_idx_valid(self, pdb_instance):
        idx = pdb_instance.get_chain_idx("A")
        assert idx == 0

    def test_get_chain_idx_case_insensitive(self, pdb_instance):
        idx_lower = pdb_instance.get_chain_idx("a")
        idx_upper = pdb_instance.get_chain_idx("A")
        assert idx_lower == idx_upper

    def test_get_chain_idx_invalid(self, pdb_instance):
        with pytest.raises(ValueError, match="Chain with id Z does not exist in PDB!"):
            pdb_instance.get_chain_idx("Z")


class TestPDBMagicMethods:
    @pytest.fixture
    def pdb_instance(self):
        atom_df = DataFrame(
            {
                "chain_id": ["A", "B"],
                "residue_number": [1, 1],
                "record_name": ["ATOM", "ATOM"],
                "atom_number": [1, 2],
                "atom_name": ["N", "N"],
                "residue_name": ["ALA", "GLY"],
                "x_coord": [0.0, 1.0],
                "y_coord": [0.0, 1.0],
                "z_coord": [0.0, 1.0],
                "occupancy": [1.0, 1.0],
                "b_factor": [0.0, 0.0],
                "segment_id": ["", ""],
                "element_symbol": ["N", "N"],
                "charge": [0, 0],
            }
        )
        return PDB({"ATOM": atom_df})

    def test_getitem_by_index(self, pdb_instance):
        chain = pdb_instance[0]
        assert chain.chain_id == "A"

    def test_getitem_by_chain_id(self, pdb_instance):
        chain = pdb_instance["A"]
        assert chain.chain_id == "A"

    def test_getitem_slice(self, pdb_instance):
        chains = pdb_instance[0:2]
        assert len(chains) == 2

    def test_setitem_by_index(self, pdb_instance):
        mock_chain = Mock(spec=Chain)
        pdb_instance[0] = mock_chain
        assert pdb_instance[0] == mock_chain

    def test_setitem_by_chain_id(self, pdb_instance):
        mock_chain = Mock(spec=Chain)
        pdb_instance["A"] = mock_chain
        assert pdb_instance["A"] == mock_chain

    def test_len(self, pdb_instance):
        assert len(pdb_instance) == 2
