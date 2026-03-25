import tempfile
from pathlib import Path
from unittest.mock import Mock, patch
import pytest
from resequencer.substitute.sub import (
    Substitution,
    load_substitution_file,
    pdb_substitution,
)


class TestSubstitution:
    def test_substitution_creation(self):
        sub = Substitution(chain="A", residue=42, base="A", new_base="G")
        assert sub.chain == "A"
        assert sub.residue == 42
        assert sub.base == "A"
        assert sub.new_base == "G"


class TestLoadSubstitutionFile:
    def test_load_comma_separated_values(self):
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
            f.write("A,10,A,G\n")
            f.write("B,20,T,C\n")
            f.flush()
            result = load_substitution_file(Path(f.name))
        assert len(result) == 2
        assert result[10].chain == "A"
        assert result[10].new_base == "G"
        assert result[20].chain == "B"
        Path(f.name).unlink()

    def test_load_whitespace_separated_values(self):
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
            f.write("A 10 A G\n")
            f.write("B 20 T C\n")
            f.flush()
            result = load_substitution_file(Path(f.name))
        assert len(result) == 2
        assert result[10].residue == 10
        assert result[20].residue == 20
        Path(f.name).unlink()

    def test_skip_comments_and_empty_lines(self):
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
            f.write("# This is a comment\n")
            f.write("\n")
            f.write("A,10,A,G\n")
            f.write("  \n")
            f.write("# Another comment\n")
            f.write("B,20,T,C\n")
            f.flush()
            result = load_substitution_file(Path(f.name))
        assert len(result) == 2
        Path(f.name).unlink()

    def test_ignore_malformed_lines(self):
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
            f.write("A,10,A\n")  # 3 fields
            f.write("A,10,A,G,extra\n")  # 5 fields
            f.write("A,10,A,G\n")  # valid
            f.flush()
            result = load_substitution_file(Path(f.name))
        assert len(result) == 1
        assert 10 in result
        Path(f.name).unlink()

    def test_empty_file(self):
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
            f.flush()
            result = load_substitution_file(Path(f.name))
        assert result == {}
        Path(f.name).unlink()


class TestPdbSubstitution:
    def test_file_not_found(self):
        pdb = Mock()
        with pytest.raises(FileNotFoundError):
            pdb_substitution(pdb, Path("/nonexistent/file.txt"))

    def test_substitution_applied(self):
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
            f.write("A,10,A,G\n")
            f.write("B,20,T,C\n")
            f.flush()

            mock_residue_a = Mock()
            mock_residue_b = Mock()
            mock_pdb = {"A": {9: mock_residue_a}, "B": {19: mock_residue_b}}

            with patch("resequencer.substitute.sub.logging.info"):
                pdb_substitution(mock_pdb, Path(f.name))  # type: ignore

            mock_residue_a.substitute.assert_called_once_with("A", "G")
            mock_residue_b.substitute.assert_called_once_with("T", "C")
        Path(f.name).unlink()

    def test_logging_output(self):
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
            f.write("A,10,A,G\n")
            f.flush()

            mock_residue = Mock()
            mock_pdb = {"A": {9: mock_residue}}

            with patch("resequencer.substitute.sub.logging.info") as mock_log:
                pdb_substitution(mock_pdb, Path(f.name))  # type: ignore
                mock_log.assert_called_once()
                assert "Chain A" in mock_log.call_args[0][0]
                assert "A" in mock_log.call_args[0][0]
                assert "G" in mock_log.call_args[0][0]
        Path(f.name).unlink()
