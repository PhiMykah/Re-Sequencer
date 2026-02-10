from .record import PDBRecord


class PDBResidue:
    # ---------------------------------------------------------------------------- #
    #                        Initialization and Constructor                        #
    # ---------------------------------------------------------------------------- #
    def __init__(
        self,
        pdb_records: list[PDBRecord],
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
    def records(self) -> list[PDBRecord]:
        return self._records

    @records.setter
    def records(self, value: list[PDBRecord]) -> None:
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

    def append(self, pdb_record: PDBRecord) -> None:
        if not self._records:
            return

        new_records = self._records + [pdb_record]
        self._records = new_records
        starting_idx = self._records[0].line_idx
        self._set_records(
            self._records,
            starting_idx,
            self._residue_name,
            self._residue_number,
            self.chain_id,
        )

    def insert(self, pdb_record: PDBRecord, idx: int) -> None:
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
        self._set_records(
            self._records,
            starting_idx,
            self._residue_name,
            self._residue_number,
            self.chain_id,
        )

    def remove(self, idx: int) -> None:
        if not self._records:
            return

        removal_position = idx % len(self._records)
        self._records.pop(removal_position)
        starting_idx = self._records[0].line_idx
        self._set_records(
            self._records,
            starting_idx,
            self._residue_name,
            self._residue_number,
            self.chain_id,
        )

    # ---------------------------------------------------------------------------- #
    #                               Helper Functions                               #
    # ---------------------------------------------------------------------------- #

    def _set_records(
        self,
        pdb_records: list[PDBRecord],
        starting_idx: int,
        residue_name: str,
        residue_number: int,
        chain_id: str,
    ):
        if not pdb_records:
            return
        self._records: list[PDBRecord] = []
        for idx, pdb_record in enumerate(pdb_records):
            pdb_record.line_idx = starting_idx + idx
            pdb_record.residue_name = residue_name
            pdb_record.residue_number = residue_number
            pdb_record.chain_id = chain_id

            self.records.append(pdb_record)

    # ---------------------------------------------------------------------------- #
    #                                 Magic Methods                                #
    # ---------------------------------------------------------------------------- #

    def __len__(self) -> int:
        return len(self._records)
