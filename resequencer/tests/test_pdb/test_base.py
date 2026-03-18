import pytest

from resequencer.pdb import base


class DummyRecord:
    def __init__(self, atom_name):
        self.atom_name = atom_name


def test_base_to_int_enum():
    assert base._base_to_int(base.NucleotideBase.A) == 0
    assert base._base_to_int(base.NucleotideBase.G) == 1
    assert base._base_to_int(base.NucleotideBase.C) == 2
    assert base._base_to_int(base.NucleotideBase.U) == 3
    assert base._base_to_int(base.NucleotideBase.T) == 4


def test_base_to_int_str():
    assert base._base_to_int("A") == 0
    assert base._base_to_int("G") == 1
    assert base._base_to_int("C") == 2
    assert base._base_to_int("U") == 3
    assert base._base_to_int("T") == 4
    assert base._base_to_int("adenine") == 0
    assert base._base_to_int("guanine") == 1
    assert base._base_to_int("cytosine") == 2
    assert base._base_to_int("uracil") == 3
    assert base._base_to_int("thymine") == 4
    assert base._base_to_int("DA") == 0
    assert base._base_to_int("DG") == 1
    assert base._base_to_int("DC") == 2
    assert base._base_to_int("DT") == 4


def test_base_to_int_int():
    assert base._base_to_int(0) == 0
    assert base._base_to_int(1) == 1
    assert base._base_to_int(2) == 2
    assert base._base_to_int(3) == 3
    assert base._base_to_int(4) == 4


def test_base_to_int_invalid():
    with pytest.raises(ValueError):
        base._base_to_int("X")
    with pytest.raises(ValueError):
        base._base_to_int(-1)
    with pytest.raises(ValueError):
        base._base_to_int(99)
    # Future proofing for in 3.12 when __contains__ will not raise TypeError
    try:
        assert not base._base_to_int(None)  # type: ignore
    except TypeError:
        with pytest.deprecated_call():
            with pytest.raises(TypeError):
                base._base_to_int(None)  # type: ignore
    # Future proofing for in 3.12 when __contains__ will not raise TypeError
    try:
        assert not base._base_to_int(3.14)  # type: ignore
    except TypeError:
        with pytest.deprecated_call():
            with pytest.raises(TypeError):
                base._base_to_int(3.14)  # type: ignore


def test_swap_base_atoms_basic():
    records = [
        DummyRecord("N9"),
        DummyRecord("C4"),
        DummyRecord("C8"),
        DummyRecord("N1"),
    ]
    mapping = base.PURINE_TO_PYRIMIDINE
    swapped = base._swap_base_atoms(records, mapping)
    assert [r.atom_name for r in swapped] == ["N1", "C2", "C6", "N1"]


def test_swap_base_atoms_no_change():
    records = [DummyRecord("X1"), DummyRecord("Y2")]
    mapping = base.PURINE_TO_PYRIMIDINE
    swapped = base._swap_base_atoms(records, mapping)
    assert [r.atom_name for r in swapped] == ["X1", "Y2"]


def test_swap_base_atoms_partial():
    records = [DummyRecord("N9"), DummyRecord("Y2")]
    mapping = base.PURINE_TO_PYRIMIDINE
    swapped = base._swap_base_atoms(records, mapping)
    assert [r.atom_name for r in swapped] == ["N1", "Y2"]


def test_deleted_atoms_keys():
    # Ensure all keys in DELETED_ATOMS are valid base indices
    for k in base.DELETED_ATOMS.keys():
        assert k in [e.value for e in base.NucleotideBase]


def test_purine_and_pyrimidine_lists():
    assert set(base.PURINE) == {
        base.NucleotideBase.A.value,
        base.NucleotideBase.G.value,
    }
    assert set(base.PYRIMIDINE) == {
        base.NucleotideBase.C.value,
        base.NucleotideBase.U.value,
        base.NucleotideBase.T.value,
    }
