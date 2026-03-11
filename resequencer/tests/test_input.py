from pathlib import Path

import pytest
from Bio.SeqRecord import SeqRecord
from biopandas.pdb import PandasPdb

from resequencer.io.input import get_input_fasta, get_input_pdb

# Assuming data folder exists in the same directory
DATA_DIR = Path(__file__).parent / "test_data"
VALID_PDB_FILE = DATA_DIR / "1BNA.pdb"
VALID_FASTA_FILE = DATA_DIR / "rcsb_pdb_1BNA.fasta"
ALT_FASTA_FILE = DATA_DIR / "rcsb_pdb_1LMB.fasta"
INVALID_PDB_FILE = DATA_DIR / "nonexistent.pdb"
INVALID_FASTA_FILE = DATA_DIR / "nonexistent.fasta"


class TestGetInputPdb:
    def test_get_input_pdb_from_file(self):
        """Test loading PDB from valid file path"""
        result = get_input_pdb(VALID_PDB_FILE)
        assert isinstance(result, PandasPdb)

    def test_get_input_pdb_from_valid_code(self):
        """Test fetching PDB from valid PDB code"""
        result = get_input_pdb("1MBN")
        assert isinstance(result, PandasPdb)

    def test_get_input_pdb_file_not_found(self):
        """Test error handling for missing file"""
        with pytest.raises(FileNotFoundError):
            get_input_pdb(INVALID_PDB_FILE)

    def test_get_input_pdb_invalid_code(self):
        """Test error handling for invalid PDB code"""
        with pytest.raises(ValueError):
            get_input_pdb("XXXX")

    def test_get_input_pdb_string_path(self):
        """Test loading PDB from string path"""
        result = get_input_pdb(str(VALID_PDB_FILE))
        assert isinstance(result, PandasPdb)


class TestGetInputFasta:
    def test_get_input_fasta_valid_file(self):
        """Test loading FASTA from valid file"""
        result = get_input_fasta(VALID_FASTA_FILE)
        assert isinstance(result, dict)
        assert all(isinstance(v, SeqRecord) for v in result.values())

    def test_get_input_fasta_file_not_found(self):
        """Test error handling for missing FASTA file"""
        with pytest.raises(FileNotFoundError):
            get_input_fasta(INVALID_FASTA_FILE)

    def test_get_input_fasta_returns_dict(self):
        """Test that FASTA returns properly formatted dictionary"""
        result = get_input_fasta(VALID_FASTA_FILE)
        assert len(result) > 0
        assert all(isinstance(k, str) for k in result.keys())

    def test_get_input_fasta_multiline_file(self):
        """Tests that FASTA returns single and multi-line files"""
        result = get_input_fasta(VALID_FASTA_FILE)
        assert len(result) == 1
        multi_result = get_input_fasta(ALT_FASTA_FILE)
        assert len(multi_result) == 3
