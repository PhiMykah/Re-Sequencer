import pytest
from resequencer.pdb.residue import Residue
from resequencer.pdb.atom_record import AtomRecord


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


def test_residue_initialization(atom_records):
    residue = Residue(atom_records, 1, "DA", 1, "A")
    assert len(residue.records) == 3
    assert residue.residue_name == "DA"
    assert residue.residue_number == 1
    assert residue.chain_id == "A"
    assert residue.starting_idx == 1


def test_append(atom_records):
    residue = Residue(atom_records, 1, "DA", 1, "A")
    new_atom = AtomRecord(
        record_name="test",
        atom_number=2,
        atom_name="C4",
        line_idx=4,
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
    )
    residue.append(new_atom)
    assert len(residue.records) == 4
    assert residue.records[-1].atom_name == "C4"


def test_insert(atom_records):
    residue = Residue(atom_records, 1, "DA", 1, "A")
    new_atom = AtomRecord(
        record_name="test",
        atom_number=2,
        atom_name="C4",
        line_idx=4,
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
    )
    residue.insert(new_atom, 1)
    assert len(residue.records) == 4
    assert residue.records[1].atom_name == "C4"


def test_remove(atom_records):
    residue = Residue(atom_records, 1, "DA", 1, "A")
    residue.remove(1)
    assert len(residue.records) == 2
    assert residue.records[0].atom_name == "N1"
    assert residue.records[1].atom_name == "N8"


def test_substitute(atom_records):
    residue = Residue(atom_records, 1, "DA", 1, "A")
    residue.substitute("DA", "DT")
    assert residue.residue_name == "DT"


def test_print(atom_records):
    residue = Residue(atom_records, 1, "DA", 1, "A")
    assert residue.print() == "N1 C2 N8"


def test_getitem(atom_records):
    residue = Residue(atom_records, 1, "DA", 1, "A")
    assert residue[0].atom_name == "N1"  # type: ignore
    assert isinstance(residue[1:], list)
    assert residue[1].atom_name == "C2"  # type: ignore


def test_len(atom_records):
    residue = Residue(atom_records, 1, "DA", 1, "A")
    assert len(residue) == 3


def test_setters(atom_records):
    residue = Residue(atom_records, 1, "DA", 1, "A")
    residue.residue_name = "DT"
    assert residue.residue_name == "DT"
    residue.residue_number = 2
    assert residue.residue_number == 2
    residue.chain_id = "B"
    assert residue.chain_id == "B"
    residue.starting_idx = 5
    assert residue.starting_idx == 5
