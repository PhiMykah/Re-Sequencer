import pytest
from pathlib import Path
from unittest.mock import Mock, patch
from resequencer.addition.add import Addition, ChainRole, BASE_PATTERN, pdb_addition


class TestBasePattern:
    """Test BASE_PATTERN regex for DNA/RNA bases."""

    def test_dna_bases(self):
        assert BASE_PATTERN.findall("DAGCTDAGCT") == [
            "DA",
            "G",
            "C",
            "T",
            "DA",
            "G",
            "C",
            "T",
        ]

    def test_rna_bases(self):
        assert BASE_PATTERN.findall("AGCUAGCU") == [
            "A",
            "G",
            "C",
            "U",
            "A",
            "G",
            "C",
            "U",
        ]

    def test_mixed_case(self):
        assert BASE_PATTERN.findall("DaGcDaGc") == ["Da", "G", "c", "Da", "G", "c"]


class TestChainRole:
    """Test ChainRole enum."""

    def test_target_chain_role(self):
        assert ChainRole.TARGET.value == 0

    def test_other_chain_role(self):
        assert ChainRole.OTHER.value == 1


class TestAddition:
    """Test Addition dataclass."""

    def test_addition_initialization(self):
        addition = Addition(
            chains=("A", "B"),
            old_seq=(["DA", "G"], ["C", "T"]),
            new_seq=(["DA", "G", "C"], ["C", "T", "U"]),
            original_geometry="B-form",
            target_chain="A",
        )
        assert addition.chains == ("A", "B")
        assert addition.target_chain == "A"
        assert addition.original_geometry == "B-form"

    def test_load_addition_file_valid(self, tmp_path):
        file = tmp_path / "additions.txt"
        file.write_text("A DAG CT\nB C TU\n")

        additions = Addition.load_addition_file(file, "B-form")

        assert len(additions) == 1
        assert additions[0].chains == ("a", "b")
        assert additions[0].target_chain == "a"
        assert additions[0].original_geometry == "B-form"

    def test_load_addition_file_with_comments(self, tmp_path):
        file = tmp_path / "additions.txt"
        file.write_text("# Comment\nA DAG CT\nB C TU\n# Another comment\n")

        additions = Addition.load_addition_file(file, "A-form")

        assert len(additions) == 1

    def test_load_addition_file_empty_lines(self, tmp_path):
        file = tmp_path / "additions.txt"
        file.write_text("A DAG CT\n\nB C TU\n\n")

        additions = Addition.load_addition_file(file, "Z-form")

        assert len(additions) == 1

    def test_load_addition_file_invalid_values(self, tmp_path):
        file = tmp_path / "additions.txt"
        file.write_text("A DAG\n")

        with pytest.raises(ValueError):
            Addition.load_addition_file(file, "B-form")

    def test_load_addition_file_multiple(self, tmp_path):
        file = tmp_path / "additions.txt"
        file.write_text("A DAG CT\nB C TU\nC DAGC AGCU\nD U AG\n")

        additions = Addition.load_addition_file(file, "B-form")

        assert len(additions) == 2
        assert additions[0].target_chain == "a"
        assert additions[1].target_chain == "c"


class TestPdbAddition:
    """Test pdb_addition function."""

    @patch("resequencer.addition.add.run_pymol")
    @patch("resequencer.addition.add.run_x3dna")
    @patch("resequencer.addition.add.Addition.load_addition_file")
    def test_pdb_addition_file_not_found(self, mock_load, mock_x3dna, mock_pymol):
        pdb = Mock()
        input_path = Path("input.pdb")
        add_input = Path("nonexistent.txt")
        output_path = Path("output")

        with pytest.raises(FileNotFoundError):
            pdb_addition(pdb, input_path, add_input, output_path, "B-form")

    @patch("resequencer.addition.append.append_addition")
    @patch("resequencer.addition.add.run_pymol")
    @patch("resequencer.addition.add.run_x3dna")
    @patch("resequencer.addition.add.Addition.load_addition_file")
    def test_pdb_addition_success(
        self, mock_load, mock_x3dna, mock_pymol, mock_append, tmp_path
    ):
        pdb = Mock()
        pdb.reindex_atom_num = Mock()
        input_path = Path("input.pdb")
        add_input = tmp_path / "additions.txt"
        add_input.write_text("A DAG CT\nB C TU\n")
        output_path = tmp_path / "output"

        mock_addition = Mock(spec=Addition)
        mock_load.return_value = [mock_addition]
        mock_x3dna.return_value = (Mock(), "forward", False)
        mock_pymol.return_value = {"range1": (1, 5)}

        pdb_addition(pdb, input_path, add_input, output_path, "B-form")

        assert output_path.exists()
        mock_x3dna.assert_called_once()
        mock_pymol.assert_called_once()
        mock_append.assert_called_once()
        pdb.reindex_atom_num.assert_called_once_with(1)
