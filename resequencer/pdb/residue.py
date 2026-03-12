from .atom_record import AtomRecord
from .base import (
    DELETED_ATOMS,
    PURINE,
    PURINE_TO_PYRIMIDINE,
    PYRIMIDINE,
    PYRIMIDINE_TO_PURINE,
    NucleotideBase,
    _base_to_int,
    _swap_base_atoms,
)


class Residue:
    """
    Class representation of Residue in PDB and PandasPDB.

    Attributes
    ----------
    records : list[AtomRecord]
        List of Atoms in Residue.
    starting_idx : int
        Starting index of first atom in residue.
    residue_name : str
        Name of given residue.
    residue_number : int
        Current index of residue based on chain.
    chain_id : str
        Chain that residue is contained within.
    """

    # ---------------------------------------------------------------------------- #
    #                        Initialization and Constructor                        #
    # ---------------------------------------------------------------------------- #
    def __init__(
        self,
        pdb_records: list[AtomRecord],
        starting_idx: int,
        residue_name: str,
        residue_number: int,
        chain_id: str,
    ) -> None:
        self._records = []
        self._residue_name = residue_name
        self._residue_number = residue_number
        self._chain_id = chain_id
        self.starting_idx = starting_idx
        self._set_records(
            pdb_records, starting_idx, residue_name, residue_number, chain_id
        )

    # ---------------------------------------------------------------------------- #
    #                              Getters and Setters                             #
    # ---------------------------------------------------------------------------- #

    # ---------------------------------- records --------------------------------- #
    @property
    def records(self) -> list[AtomRecord]:
        return self._records

    @records.setter
    def records(self, value: list[AtomRecord]) -> None:
        self._set_records(
            value,
            self._starting_idx,
            self._residue_name,
            self._residue_number,
            self._chain_id,
        )

    # ------------------------------- starting_idx ------------------------------- #
    @property
    def starting_idx(self) -> int:
        return self._starting_idx

    @starting_idx.setter
    def starting_idx(self, value: int) -> None:
        if self._records:
            self._set_records(
                self._records,
                value,
                self._residue_name,
                self._residue_number,
                self._chain_id,
            )
        self._starting_idx = value

    # ------------------------------- residue_name ------------------------------- #
    @property
    def residue_name(self) -> str:
        return self._residue_name

    @residue_name.setter
    def residue_name(self, value: str) -> None:
        if self._records:
            self._set_records(
                self._records,
                self._starting_idx,
                value,
                self._residue_number,
                self._chain_id,
            )

    # ------------------------------ residue_number ------------------------------ #
    @property
    def residue_number(self) -> int:
        return self._residue_number

    @residue_number.setter
    def residue_number(self, value: int) -> None:
        if self._records:
            self._set_records(
                self._records,
                self._starting_idx,
                self._residue_name,
                value,
                self.chain_id,
            )
        self._residue_number: int = value

    # --------------------------------- chain_id --------------------------------- #
    @property
    def chain_id(self) -> str:
        return self._chain_id

    @chain_id.setter
    def chain_id(self, value: str) -> None:
        if self._records:
            self._set_records(
                self._records,
                self._starting_idx,
                self._residue_name,
                self._residue_number,
                value,
            )
        self._chain_id: str = value

    # ---------------------------------------------------------------------------- #
    #                               Modifying Record                               #
    # ---------------------------------------------------------------------------- #

    def append(self, pdb_record: AtomRecord) -> None:
        """
        Append a new atom record to the residue and re-index.

        Parameters
        ----------
        pdb_record : AtomRecord
            The atom record to append to this residue.
        """

        if not self._records:
            return

        new_records = self._records + [pdb_record]
        self._records = new_records
        starting_idx = self._records[0].line_idx
        # re-index after appending new record
        self._set_records(
            self._records,
            starting_idx,
            self._residue_name,
            self._residue_number,
            self.chain_id,
        )

    def insert(self, pdb_record: AtomRecord, idx: int) -> None:
        """
        Insert a new atom record to the residue at the target index and re-index.

        Parameters
        ----------
        pdb_record : AtomRecord
            The atom record to append to this residue.
        idx : int
            Position to add new atom record.
        """

        if not self._records:
            return
        insert_position = idx % len(self._records)
        new_records = (
            self._records[0:insert_position]
            + [pdb_record]
            + self._records[insert_position:]
        )
        self._records = new_records
        starting_idx = self._records[0].line_idx
        # re-index after inserting new record
        self._set_records(
            self._records,
            starting_idx,
            self._residue_name,
            self._residue_number,
            self.chain_id,
        )

    def remove(self, idx: int) -> None:
        """
        Remove an atom record from the residue and re-index.

        Parameters
        ----------
        idx : int
            Position to remove atom record.
        """

        if not self._records:
            return

        removal_position = idx % len(self._records)
        self._records.pop(removal_position)
        starting_idx = self._records[0].line_idx
        # re-index after removal of record
        self._set_records(
            self._records,
            starting_idx,
            self._residue_name,
            self._residue_number,
            self.chain_id,
        )

    # ---------------------------------------------------------------------------- #
    #                                Modify Residue                                #
    # ---------------------------------------------------------------------------- #

    def substitute(
        self,
        old_base: str,
        new_base: str,
    ) -> None:
        """
        Substitutes the base in the residue with a new base. If update_residue is True, modifies the current residue in place.

        Otherwise, returns a new Residue with the substitution.

        Parameters
        ----------
        old_base : str
            Current base to resequence from.
        new_base : str
            New base for current Residue object.
        """
        # Remove the base atoms before conversion
        self._delete_base_atoms(old_base)
        # Change remaining base atoms to new base
        self._change_base_atoms(old_base, new_base)
        # Update the residue name
        self.residue_name = new_base

    # ---------------------------------------------------------------------------- #
    #                                Display Residue                               #
    # ---------------------------------------------------------------------------- #

    def print(self) -> str:
        """
        Generate a string representation of atom names in the residue.

        Returns
        -------
        str
            A space-separated string of atom names from all records in the residue.

            Example:
            residue.print()
            'DA DT DG DC'
        """

        print_list = []
        for rec in self._records:
            print_list.append(rec.atom_name)

        return " ".join(print_list)

    # ---------------------------------------------------------------------------- #
    #                               Helper Functions                               #
    # ---------------------------------------------------------------------------- #

    def _set_records(
        self,
        pdb_records: list[AtomRecord],
        starting_idx: int,
        residue_name: str,
        residue_number: int,
        chain_id: str,
    ) -> None:
        """
        Set atom record list for Residue and ensure proper indexing for line_idx.

        Parameters
        ----------
        pdb_records : list[AtomRecord]
            List of Atoms in Residue.
        starting_idx : int
            Starting index of first atom in residue.
        residue_name : str
            Name of given residue.
        residue_number : int
            Current index of residue based on chain.
        chain_id : str
            Chain that residue is contained within.
        """
        if not pdb_records:
            return
        self._records: list[AtomRecord] = []
        for idx, pdb_record in enumerate(pdb_records):
            pdb_record.line_idx = starting_idx + idx
            pdb_record.residue_name = residue_name
            pdb_record.residue_number = residue_number
            pdb_record.chain_id = chain_id

            self.records.append(pdb_record)
        self._starting_idx = starting_idx
        self._residue_name = residue_name
        self._residue_number = residue_number
        self._chain_id = chain_id

    def _delete_base_atoms(self, base: "NucleotideBase | str | int") -> None:
        """
        Removes atoms corresponding to a specified nucleotide base from a DataFrame.

        Parameters
        ----------
        base : "NucleotideBase | str | int"
            The nucleotide base to remove atoms for. Can be a NucleotideBase enum, string, or integer.

        Raises
        ------
        ValueError
            If the provided nucleotide base is not recognized or out of range.
        """
        nucleotide_base: int = _base_to_int(base)

        if nucleotide_base not in DELETED_ATOMS.keys():
            raise ValueError(
                f"Nucleotide base out of range; got {nucleotide_base}",
            )

        new_records = list(
            filter(
                lambda rec: rec.atom_name not in DELETED_ATOMS[nucleotide_base],
                self._records,
            )
        )

        self._set_records(
            new_records,
            new_records[0].line_idx,
            self.residue_name,
            self.residue_number,
            self.chain_id,
        )

    def _change_base_atoms(
        self,
        old_base: "NucleotideBase | str | int",
        new_base: "NucleotideBase | str | int",
    ) -> None:
        """
        Change the base atoms in a nucleotide DataFrame from one base type to another.
        This function handles conversions between purines and pyrimidines.
        If the conversion is not between purine and pyrimidine (or vice versa),
        the residue is unchanged.

        Parameters
        ----------
        old_base : "NucleotideBase | str | int"
            The original nucleotide base to be changed.
        new_base : "NucleotideBase | str | int"
            The target nucleotide base to change to.

        Raises
        ------
        ValueError
            If either the input or output base is not a recognized nucleotide base.
        """
        input_base: int = _base_to_int(old_base)
        output_base: int = _base_to_int(new_base)

        # Ensure input and output base are either purine or pyrimidine
        if input_base not in PURINE + PYRIMIDINE:
            raise ValueError(f"Input Nucleotide base out of range; got {input_base}")
        if output_base not in PURINE + PYRIMIDINE:
            raise ValueError(f"Input Nucleotide base out of range; got {output_base}")

        if input_base in PURINE and output_base in PYRIMIDINE:
            new_records: list[AtomRecord] = _swap_base_atoms(
                self._records, PURINE_TO_PYRIMIDINE
            )
        elif input_base in PYRIMIDINE and output_base in PURINE:
            new_records: list[AtomRecord] = _swap_base_atoms(
                self._records, PYRIMIDINE_TO_PURINE
            )
        else:
            return

        self._set_records(
            new_records,
            new_records[0].line_idx,
            self.residue_name,
            self.residue_number,
            self.chain_id,
        )

    # ---------------------------------------------------------------------------- #
    #                                 Magic Methods                                #
    # ---------------------------------------------------------------------------- #

    def __len__(self) -> int:
        """
        Get length of residue based on number of atom records.

        Returns
        -------
        int
            Number of atom records.
        """
        return len(self._records)

    def __getitem__(self, indices) -> AtomRecord | list[AtomRecord]:
        """
        Residue indexer.

        Parameters
        ----------
        indices : int | slice
            Index or indices of records.

        Returns
        -------
        AtomRecord | list[AtomRecord]
            Current record(s) at provided index/indices.
        """
        return self._records[indices]
