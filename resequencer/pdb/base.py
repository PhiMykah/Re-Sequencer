from enum import Enum


class NucleotideBase(Enum):
    """
    Enumeration of nucleotide bases with multiple aliases for DNA and RNA bases.
    """

    A = 0
    G = 1
    C = 2
    U = 3
    T = 4
    adenine = A
    guanine = G
    cytosine = C
    uracil = U
    thymine = 4
    DA = A
    DG = G
    DC = C
    DT = T


# Maps nucleotide base indices to lists of atom names that are deleted for each base
DELETED_ATOMS: "dict[int, list[str]]" = {
    0: ["N1", "C2", "N3", "C5", "C6", "N6", "N7"],  # A
    1: ["N1", "C2", "N2", "N3", "C5", "C6", "O6", "N7"],  # G
    2: ["O2", "N3", "C4", "N4", "C5"],  # C
    3: ["O2", "N3", "C4", "O4", "C5"],  # U
    4: ["O2", "N3", "C4", "O4", "C5", "C7"],  # T
}

# Maps purine atom names to their corresponding pyrimidine atom names.
PURINE_TO_PYRIMIDINE: "dict[str, str]" = {
    "N9": "N1",
    "C4": "C2",
    "C8": "C6",
}

# Maps pyrimidine atom names to their corresponding purine atom names.
PYRIMIDINE_TO_PURINE: "dict[str, str]" = {
    "N1": "N9",
    "C6": "C8",
    "C2": "C4",
}

# List of indices corresponding to purine bases (Adenine and Guanine).
PURINE: list[int] = [0, 1]

# List of indices corresponding to pyrimidine bases (Cytosine, Uracil, Thymine).
PYRIMIDINE: list[int] = [2, 3, 4]


def _base_to_int(base: "NucleotideBase | str | int") -> int:
    """
    Converts a nucleotide base representation to its corresponding integer value.

    Parameters
    ----------
        base : NucleotideBase | str | int
        The nucleotide base to convert. This can be:

            - An instance of the NucleotideBase enum,
            - A string representing the base (case-insensitive, must match a valid NucleotideBase member),
            - An integer corresponding to a valid NucleotideBase value.
    Returns
    -------
        int
            The integer value corresponding to the provided nucleotide base.
    Raises
    ------
        ValueError
            If the input is not a valid NucleotideBase, string, or integer representation,
            or if the integer is negative.
    """

    # Test for string instance
    if isinstance(base, str):
        base = base.lower()
        if base in [val.lower() for val in NucleotideBase.__members__.keys()]:
            nucleotide_base: int = NucleotideBase[base.upper()].value
        else:
            valid_bases: str = ",".join(
                [val.lower() for val in NucleotideBase.__members__.keys()]
                + list({str(val.value) for val in NucleotideBase.__members__.values()})
            )
            raise ValueError(f"Invalid base type {base}! Must be {valid_bases}")
    # Test for int instance
    elif isinstance(base, int):
        if base >= 0:
            nucleotide_base: int = NucleotideBase(base).value
        else:
            raise ValueError(f"base must be a positive int; got {base}")
    # Test for NucleotideBase enum
    elif base in NucleotideBase:
        nucleotide_base: int = base.value
    # Raise error if base does not match any valid instances
    else:
        raise ValueError(
            f"base must be an instance of {type(NucleotideBase)}, str, or int; got {type(base)}"
        )

    return nucleotide_base


def _swap_base_atoms(records: list, mapping: dict[str, str]) -> list:
    """
    Swap atom names in an object list according to a provided mapping.
    This function replaces values in the 'atom_name' attribute of the given object list
    based on the specified mapping dictionary.

    Parameters
    ----------
        records : list
            The input list containing an 'atom_name' column.
        mapping : dict[str, str]
            A dictionary mapping original atom names to new atom names.

    Returns
    -------
        list
            Updated list with new atom names.
    """
    for rec in records:
        if rec.atom_name in mapping:
            rec.atom_name = mapping[rec.atom_name]
    return records
