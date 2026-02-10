from pandas import DataFrame, to_numeric

from resequencer.pdb.record import PDBRecord
from .residue import PDBResidue


class PDBChain:
    # ---------------------------------------------------------------------------- #
    #                        Initialization and Constructor                        #
    # ---------------------------------------------------------------------------- #
    def __init__(
        self, residues: list[PDBResidue], chain_id: str, starting_idx: int = -1
    ) -> None:
        self._residues = []
        self._chain_id = chain_id
        self._starting_idx = starting_idx
        self._set_residues(residues, chain_id, starting_idx)
        pass

    # ---------------------------------------------------------------------------- #
    #                              Getters and Setters                             #
    # ---------------------------------------------------------------------------- #

    # --------------------------------- residues --------------------------------- #
    @property
    def residues(self) -> list[PDBResidue]:
        return self._residues

    @residues.setter
    def residues(self, value: list[PDBResidue]) -> None:
        self._set_residues(value, self._chain_id, self._starting_idx)

    # --------------------------------- chain_id --------------------------------- #
    @property
    def chain_id(self) -> str:
        return self._chain_id

    @chain_id.setter
    def chain_id(self, value) -> None:
        if self._residues:
            self._set_residues(self._residues, value, self._starting_idx)
        self._chain_id = value

    # ------------------------------- starting_idx ------------------------------- #
    @property
    def starting_idx(self) -> int:
        return self._starting_idx

    @starting_idx.setter
    def starting_idx(self, value: int) -> None:
        if self._residues:
            self._set_residues(self._residues, self._chain_id, value)
        self._starting_idx = value

    # ---------------------------------------------------------------------------- #
    #                                Modifying Chain                               #
    # ---------------------------------------------------------------------------- #

    def append(self, residue: PDBResidue) -> None:
        pass

    def insert(self, residue: PDBResidue, idx: int) -> None:
        pass

    def remove(self, idx: int) -> None:
        pass

    # ---------------------------------------------------------------------------- #
    #                            Characteristic Methods                            #
    # ---------------------------------------------------------------------------- #

    def total_length(self) -> int:
        return sum(len(res) for res in self._residues)

    # ---------------------------------------------------------------------------- #
    #                               Helper Functions                               #
    # ---------------------------------------------------------------------------- #

    def _set_residues(
        self, residues: list[PDBResidue], chain_id: str, starting_idx: int
    ) -> None:
        if not residues:
            return
        self._residues: list[PDBResidue] = []
        for idx, residue in enumerate(residues):
            if starting_idx == -1 and idx == 0:
                self._starting_idx = residue.starting_idx
                starting_idx = residue.starting_idx
            residue.chain_id = chain_id
            # residue.residue_number = starting_idx + idx
            self._residues.append(residue)

    # ---------------------------------------------------------------------------- #
    #                                 Magic Methods                                #
    # ---------------------------------------------------------------------------- #

    def __len__(self) -> int:
        return len(self._residues)


def dataframe_to_chains(atom: DataFrame):
    chains: list[PDBChain] = []
    residues: list[PDBResidue] = []
    records: list[PDBRecord] = []

    chain_groups: dict = {
        chain_id: chain_data for chain_id, chain_data in atom.groupby("chain_id")
    }
    for chain_id, chain in chain_groups.items():
        residue_groups: dict = {
            residue_number: residue_data
            for residue_number, residue_data in chain.groupby("residue_number")
        }
        for res_id, res in residue_groups.items():
            starting_idx: int = -1
            residue_name: str = ""
            first = True
            for row in res.itertuples():
                if first:
                    starting_idx = int(to_numeric(row.Index))
                    residue_name = str(row.residue_name)
                    first = False
                records.append(
                    PDBRecord(
                        int(to_numeric(row.Index)),
                        str(row.record_name),
                        int(to_numeric(row.atom_number)),
                        str(row.atom_name),
                        str(row.residue_name),
                        str(row.chain_id),
                        int(to_numeric(row.residue_number)),
                        to_numeric(row.x_coord),
                        to_numeric(row.y_coord),
                        to_numeric(row.z_coord),
                        to_numeric(row.occupancy),
                        to_numeric(row.b_factor),
                        str(row.segment_id),
                        str(row.element_symbol),
                        to_numeric(row.charge),
                    )
                )
            residues.append(
                PDBResidue(
                    records, starting_idx, residue_name, int(res_id), str(chain_id)
                )
            )
            records = []

        chains.append(PDBChain(residues, str(chain_id)))
        residues = []

    return chains
    # chain_id = ""
    # curr_residue_number: int = -1
    # curr_chain_id: str = ""

    # for row in atom.itertuples():
    #     idx: int = int(to_numeric(row.Index))
    #     record_name: str = str(row.record_name)
    #     atom_number: int = int(to_numeric(row.atom_number))
    #     atom_name: str = str(row.atom_name)
    #     residue_name: str = str(row.residue_name)
    #     chain_id = str(row.chain_id)
    #     residue_number = int(to_numeric(row.residue_number))
    #     x_coord: float = to_numeric(row.x_coord)
    #     y_coord: float = to_numeric(row.y_coord)
    #     z_coord: float = to_numeric(row.z_coord)
    #     occupancy: float = to_numeric(row.occupancy)
    #     b_factor: float = to_numeric(row.b_factor)
    #     segment_id: str = str(row.segment_id)
    #     element_symbol: str = str(row.element_symbol)
    #     charge: float = to_numeric(row.charge)

    #     if idx == 0:
    #         curr_residue_number = residue_number
    #         curr_chain_id = chain_id

    #     records.append(
    #         PDBRecord(
    #             idx,
    #             record_name,
    #             atom_number,
    #             atom_name,
    #             residue_name,
    #             chain_id,
    #             residue_number,
    #             x_coord,
    #             y_coord,
    #             z_coord,
    #             occupancy,
    #             b_factor,
    #             segment_id,
    #             element_symbol,
    #             charge,
    #         )
    #     )

    #     if chain_id != curr_chain_id:
    #         curr_chain_id = chain_id
    #         if residues:
    #             chains.append(
    #                 PDBChain(residues, residues[0].chain_id, residues[0].starting_idx)
    #             )
    #         residues = []

    #     if residue_number != curr_residue_number:
    #         curr_residue_number = residue_number
    #         if records:
    #             residues.append(
    #                 PDBResidue(
    #                     records,
    #                     records[0].residue_number,
    #                     residue_name,
    #                     records[0].line_idx,
    #                     chain_id,
    #                 )
    #             )
    #         records = []

    # if residues:
    #     chains.append(
    #         PDBChain(residues, residues[0].chain_id, residues[0].starting_idx)
    #     )
    # return chains


def chains_to_dataframe(chains: list[PDBChain]) -> DataFrame:
    records = []
    for chain in chains:
        for residue in chain.residues:
            for record in residue.records:
                records.append(
                    {
                        "record_name": record.record_name,
                        "atom_number": record.atom_number,
                        "blank_1": "",
                        "atom_name": record.atom_name,
                        "alt_loc": "",
                        "residue_name": record.residue_name,
                        "blank_2": "",
                        "chain_id": record.chain_id,
                        "residue_number": record.residue_number,
                        "insertion": "",
                        "blank_3": "",
                        "x_coord": record.x_coord,
                        "y_coord": record.y_coord,
                        "z_coord": record.z_coord,
                        "occupancy": record.occupancy,
                        "b_factor": record.b_factor,
                        "blank_4": "",
                        "segment_id": record.segment_id,
                        "element_symbol": record.element_symbol,
                        "charge": record.charge,
                        "line_idx": record.line_idx,
                    }
                )
    return DataFrame(records)
