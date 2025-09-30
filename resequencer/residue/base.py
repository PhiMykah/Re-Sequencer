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
deleted_atoms: "dict[int, list[str]]" = {
    0: ["N1", "C2", "N3", "C5", "C6", "N6", "N7"], # A 
    1: ["N1", "C2", "N2", "N3", "C5", "C6", "O6", "N7"], # G
    2: ["O2", "N3", "C4", "N4", "C5"], # C
    3: ["O2", "N3", "C4", "O4", "C5"], # U
    4: ["O2", "N3", "C4", "O4", "C5", "C7"], # T
}

# Maps purine atom names to their corresponding pyrimidine atom names.
purine_to_pyrimidine: "dict[str, str]" = {
    "N9": "N1",
    "C4": "C6",
    "C8": "C2",
}

# Maps pyrimidine atom names to their corresponding purine atom names.
pyrimidine_to_purine: "dict[str, str]" = {
    "N1": "N9",
    "C6": "C4",
    "C2": "C8",
}

# List of indices corresponding to purine bases (Adenine and Guanine).
PURINE: list[int] = [0, 1]

# List of indices corresponding to pyrimidine bases (Cytosine, Uracil, Thymine).
PYRIMIDINE: list[int] = [2, 3, 4]
