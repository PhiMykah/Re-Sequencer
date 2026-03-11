import pytest
from pathlib import Path
from resequencer.parse.parse import parse_args, build_parser, CLargs


class TestBuildParser:
    """Tests for build_parser function"""

    def test_parser_creation(self):
        """Test that parser is created successfully"""
        parser = build_parser()
        assert parser is not None
        assert parser.prog == "Re-Sequencer"


class TestParseArgsBasic:
    """Tests for basic parse_args functionality"""

    def test_required_arguments_only(self):
        """Test parsing with only required arguments"""
        argv = ["--input", "1ABC", "--fasta", "test.fasta"]
        result = parse_args(argv)
        assert result.pdb_input == "1ABC"
        assert result.fasta_input == Path("test.fasta")
        assert result.nucleic_structure == ""
        assert result.sub_input is None
        assert result.add_input is None
        assert result.walk_input is None

    def test_missing_input_argument(self):
        """Test that missing --input raises error"""
        argv = ["--fasta", "test.fasta"]
        with pytest.raises(SystemExit):
            parse_args(argv)

    def test_missing_fasta_argument(self):
        """Test that missing --fasta raises error"""
        argv = ["--input", "1ABC"]
        with pytest.raises(SystemExit):
            parse_args(argv)


class TestParseArgsNucleicStructure:
    """Tests for nucleic structure argument"""

    def test_form_a(self):
        """Test a-form nucleic structure"""
        argv = ["--input", "1ABC", "--fasta", "test.fasta", "--form", "a"]
        result = parse_args(argv)
        assert result.nucleic_structure == "a"

    def test_form_b(self):
        """Test b-form nucleic structure"""
        argv = ["--input", "1ABC", "--fasta", "test.fasta", "--form", "b"]
        result = parse_args(argv)
        assert result.nucleic_structure == "b"

    def test_form_z(self):
        """Test z-form nucleic structure"""
        argv = ["--input", "1ABC", "--fasta", "test.fasta", "--form", "z"]
        result = parse_args(argv)
        assert result.nucleic_structure == "z"

    def test_invalid_form(self):
        """Test invalid nucleic structure choice"""
        argv = ["--input", "1ABC", "--fasta", "test.fasta", "--form", "x"]
        with pytest.raises(SystemExit):
            parse_args(argv)


class TestParseArgsOptionalFiles:
    """Tests for optional file arguments"""

    def test_substitution_input(self):
        """Test substitution file argument"""
        argv = ["--input", "1ABC", "--fasta", "test.fasta", "--sub", "sub.in"]
        result = parse_args(argv)
        assert result.sub_input == Path("sub.in")

    def test_chain_addition_with_form(self):
        """Test chain addition with required form argument"""
        argv = [
            "--input",
            "1ABC",
            "--fasta",
            "test.fasta",
            "--add",
            "add.chain",
            "--form",
            "b",
        ]
        result = parse_args(argv)
        assert result.add_input == Path("add.chain")
        assert result.nucleic_structure == "b"

    def test_chain_addition_without_form_error(self):
        """Test that --add without --form raises error"""
        argv = ["--input", "1ABC", "--fasta", "test.fasta", "--add", "add.chain"]
        with pytest.raises(SystemExit):
            parse_args(argv)

    def test_walk_input_with_form(self):
        """Test walk input with required form argument"""
        argv = [
            "--input",
            "1ABC",
            "--fasta",
            "test.fasta",
            "--walk",
            "file.walk",
            "--form",
            "a",
        ]
        result = parse_args(argv)
        assert result.walk_input == Path("file.walk")
        assert result.nucleic_structure == "a"

    def test_walk_input_without_form_error(self):
        """Test that --walk without --form raises error"""
        argv = ["--input", "1ABC", "--fasta", "test.fasta", "--walk", "file.walk"]
        with pytest.raises(SystemExit):
            parse_args(argv)


class TestParseArgsOutput:
    """Tests for output arguments"""

    def test_default_output_path(self):
        """Test default output path"""
        argv = ["--input", "1ABC", "--fasta", "test.fasta"]
        result = parse_args(argv)
        assert result.output_path == Path("output")

    def test_custom_output_path(self):
        """Test custom output path"""
        argv = ["--input", "1ABC", "--fasta", "test.fasta", "--output", "/custom/path"]
        result = parse_args(argv)
        assert result.output_path == Path("/custom/path")

    def test_default_output_file(self):
        """Test default output file name"""
        argv = ["--input", "1ABC", "--fasta", "test.fasta"]
        result = parse_args(argv)
        assert result.output_file == "output.pdb"

    def test_custom_output_file(self):
        """Test custom output file name"""
        argv = ["--input", "1ABC", "--fasta", "test.fasta", "--file", "custom.pdb"]
        result = parse_args(argv)
        assert result.output_file == "custom.pdb"


class TestParseArgsAlternativeFlags:
    """Tests for alternative argument flags"""

    def test_input_short_flag(self):
        """Test -input short flag"""
        argv = ["-input", "1ABC", "--fasta", "test.fasta"]
        result = parse_args(argv)
        assert result.pdb_input == "1ABC"

    def test_fasta_short_flag_f(self):
        """Test -f short flag for fasta"""
        argv = ["--input", "1ABC", "-f", "test.fasta"]
        result = parse_args(argv)
        assert result.fasta_input == Path("test.fasta")

    def test_form_short_flag(self):
        """Test -form short flag"""
        argv = ["--input", "1ABC", "--fasta", "test.fasta", "-form", "b"]
        result = parse_args(argv)
        assert result.nucleic_structure == "b"

    def test_sub_short_flag(self):
        """Test -sub short flag"""
        argv = ["--input", "1ABC", "--fasta", "test.fasta", "-sub", "sub.in"]
        result = parse_args(argv)
        assert result.sub_input == Path("sub.in")


class TestParseArgsComplex:
    """Tests for complex argument combinations"""

    def test_all_arguments(self):
        """Test with all arguments specified"""
        argv = [
            "--input",
            "1ABC",
            "--fasta",
            "test.fasta",
            "--form",
            "b",
            "--sub",
            "sub.in",
            "--add",
            "add.chain",
            "--walk",
            "file.walk",
            "--output",
            "/out/path",
            "--file",
            "final.pdb",
        ]
        result = parse_args(argv)
        assert result.pdb_input == "1ABC"
        assert result.fasta_input == Path("test.fasta")
        assert result.nucleic_structure == "b"
        assert result.sub_input == Path("sub.in")
        assert result.add_input == Path("add.chain")
        assert result.walk_input == Path("file.walk")
        assert result.output_path == Path("/out/path")
        assert result.output_file == "final.pdb"

    def test_pdb_id_vs_file_path(self):
        """Test with both PDB ID and file path as input"""
        argv1 = ["--input", "1ABC", "--fasta", "test.fasta"]
        argv2 = ["--input", "/path/to/protein.pdb", "--fasta", "test.fasta"]
        result1 = parse_args(argv1)
        result2 = parse_args(argv2)
        assert result1.pdb_input == "1ABC"
        assert result2.pdb_input == Path("/path/to/protein.pdb")

    def test_incorrect_pdb_id(self):
        """Test for Error when inputting Incorrect Length PDB ID"""
        argv1 = ["--input", "1ABCD", "--fasta", "test.fasta"]
        argv2 = ["--input", "1AB", "--fasta", "test.fasta"]
        with pytest.raises(SystemExit):
            parse_args(argv1)
        with pytest.raises(SystemExit):
            parse_args(argv2)


class TestCLargsClass:
    """Tests for CLargs namespace class"""

    def test_clargs_initialization(self):
        """Test CLargs object creation"""
        clargs = CLargs(
            pdb_input="1ABC",
            fasta_input=Path("test.fasta"),
            nucleic_structure="b",
            sub_input="",
            add_input="add.chain",
            walk_input="",
            output_path=Path("output"),
            output_file="output.pdb",
        )
        assert clargs.pdb_input == "1ABC"
        assert clargs.nucleic_structure == "b"
        assert clargs.add_input == Path("add.chain")
        assert clargs.sub_input is None
        assert clargs.walk_input is None

    def test_clargs_empty_string_conversion(self):
        """Test that empty strings convert to None for optional paths"""
        clargs = CLargs(
            pdb_input="1ABC",
            fasta_input=Path("test.fasta"),
            nucleic_structure="",
            sub_input="",
            add_input="",
            walk_input="",
            output_path=Path("output"),
            output_file="output.pdb",
        )
        assert clargs.sub_input is None
        assert clargs.add_input is None
        assert clargs.walk_input is None
