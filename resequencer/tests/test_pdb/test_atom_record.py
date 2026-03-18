import pytest

from resequencer.pdb import AtomRecord


class TestAtomRecord:
    """Test suite for AtomRecord dataclass"""

    @pytest.fixture
    def sample_atom(self):
        """Create a sample AtomRecord for testing"""
        return AtomRecord(
            line_idx=0,
            record_name="ATOM",
            atom_number=1,
            atom_name="N",
            residue_name="ALA",
            chain_id="A",
            residue_number=1,
            x_coord=1.0,
            y_coord=2.0,
            z_coord=3.0,
            occupancy=1.0,
            b_factor=20.0,
            segment_id="SEG1",
            element_symbol="N",
            charge=0.0,
        )

    def test_atom_record_creation(self, sample_atom):
        """Test basic AtomRecord instantiation"""
        assert sample_atom.atom_number == 1
        assert sample_atom.atom_name == "N"
        assert sample_atom.residue_name == "ALA"

    def test_atom_coordinates(self, sample_atom):
        """Test coordinate storage"""
        assert sample_atom.x_coord == 1.0
        assert sample_atom.y_coord == 2.0
        assert sample_atom.z_coord == 3.0

    def test_atom_record_fields(self, sample_atom):
        """Test all fields are properly assigned"""
        assert sample_atom.line_idx == 0
        assert sample_atom.record_name == "ATOM"
        assert sample_atom.chain_id == "A"
        assert sample_atom.residue_number == 1

    def test_hetatm_record(self, sample_atom):
        """Test HETATM record type"""
        hetatm = AtomRecord(
            line_idx=100,
            record_name="HETATM",
            atom_number=500,
            atom_name="O",
            residue_name="HOH",
            chain_id="B",
            residue_number=200,
            x_coord=5.0,
            y_coord=6.0,
            z_coord=7.0,
            occupancy=0.5,
            b_factor=30.0,
            segment_id="SEG2",
            element_symbol="O",
            charge=-1.0,
        )
        assert hetatm.record_name == "HETATM"
        assert hetatm.occupancy == 0.5
        assert hetatm.charge == -1.0

    def test_occupancy_range(self):
        """Test occupancy values within valid range"""
        for occupancy in [0.0, 0.5, 1.0]:
            atom = AtomRecord(
                line_idx=0,
                record_name="ATOM",
                atom_number=1,
                atom_name="CA",
                residue_name="GLY",
                chain_id="A",
                residue_number=1,
                x_coord=0.0,
                y_coord=0.0,
                z_coord=0.0,
                occupancy=occupancy,
                b_factor=0.0,
                segment_id="",
                element_symbol="C",
                charge=0.0,
            )
            assert 0.0 <= atom.occupancy <= 1.0
