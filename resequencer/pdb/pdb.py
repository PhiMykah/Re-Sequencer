from Bio.SeqRecord import SeqRecord
from biopandas.pdb.pandas_pdb import PandasPdb
from pandas import DataFrame, Series, to_numeric

from .atom_record import AtomRecord
from .chain import Chain
from .residue import Residue


class PDB:
    """
    Class representation of PDB from PandasPDB.

    Attributes
    ----------
    chains : list[Chain]
        List of Chains in PDB.
    fasta : dict[str, SeqRecord]
        Fasta dictionary representation of PDB.
    """

    # ---------------------------------------------------------------------------- #
    #                        Initialization and Constructor                        #
    # ---------------------------------------------------------------------------- #
    def __init__(
        self, dict_dataframe: dict[str, DataFrame], fasta: dict[str, SeqRecord] = {}
    ) -> None:
        atom = dict_dataframe.get("ATOM", DataFrame())
        chains: list[Chain] | None = PDB._to_chains(atom)
        if chains is not None:
            self.chains = chains
            self._chain_mapping: dict[str, int] = {}
            for idx, chain in enumerate(chains):
                self._chain_mapping[chain.chain_id.lower()] = idx
        else:
            raise ValueError("Atom Dataframe for PDB is Empty!")

        self._hetatm: list[Chain] | None = PDB._to_chains(
            dict_dataframe.get("HETATM", DataFrame())
        )
        self._anisou: DataFrame = dict_dataframe.get("ANISOU", DataFrame())
        self._others: DataFrame = dict_dataframe.get("OTHERS", DataFrame())

        self.fasta = fasta

    # ---------------------------------------------------------------------------- #
    #                              Getters and Setters                             #
    # ---------------------------------------------------------------------------- #

    # ---------------------------------- chains ---------------------------------- #
    @property
    def chains(self) -> list[Chain]:
        return self._chains

    @chains.setter
    def chains(self, value: list[Chain]) -> None:
        self._chains: list[Chain] = value

    # ----------------------------------- fasta ---------------------------------- #
    @property
    def fasta(self) -> dict[str, SeqRecord]:
        return self._fasta

    @fasta.setter
    def fasta(self, value: dict[str, SeqRecord]) -> None:
        self._fasta: dict = value

    # ---------------------------------------------------------------------------- #
    #                                  Data Return                                 #
    # ---------------------------------------------------------------------------- #

    def output_pdb(self, input: PandasPdb) -> PandasPdb:
        """
        Convert internal chain data to a PandasPdb object.
        Takes the processed chains stored in this object and constructs a new
        PandasPdb instance with the chains data, preserving metadata from the
        input PandasPdb object.

        Parameters
        ----------
        input : PandasPdb
            Source PandasPdb object containing original metadata
            (header, code, pdb_text, pdb_path).

        Returns
        -------
        PandasPdb
            A new PandasPdb object with the processed chains data and
            original metadata from the input object.
        """

        output = PandasPdb()
        output._df = self._chains_to_full_dataframe()
        output.pdb_text = input.pdb_text
        output.pdb_path = input.pdb_path
        output.header = input.header
        output.code = input.code

        return output

    # ---------------------------------------------------------------------------- #
    #                            Characteristic Methods                            #
    # ---------------------------------------------------------------------------- #

    def total_length(self) -> int:
        """
        Collect the lengh of the chain based on the length of each residue.

        Returns
        -------
        int
            Total length of chain by residue length.
        """
        return sum(chain.total_length() for chain in self.chains)

    def is_dna(self, chain_id: str) -> bool:
        """
        Determine if a specified chain in the PDB structure is DNA.

        Parameters
        ----------
        chain_id: str
            The identifier of the chain to check.

        Returns
        -------
        bool
            True if the chain is DNA, False otherwise.

        Raises
        ------
        ValueError
            If the chain_id is not found in the PDB structure's chain mapping.
        ValueError
            If the chain cannot be found in the PDB's FASTA file.
        """

        if chain_id.lower() not in self._chain_mapping.keys():
            raise ValueError(
                "Inputted Chain ID '{chain_id.upper()}' is not in PDB! (expected {",
                ".join(self._chain_mapping.keys())})",
            )

        target_key: str = ""
        for key, value in self.fasta.items():
            possible_strings = [
                f"chains {chain_id}",
                f"chain {chain_id}",
                f", {chain_id}",
            ]
            if any(s in value.description.lower() for s in possible_strings):
                target_key = key

        if target_key == "":
            raise ValueError(
                f"Unable to find Chain in PDB's fasta File! Looking for Chain {chain_id.upper()}"
            )

        if "-d(" in self.fasta[target_key].description.lower():
            return True
        return False

    # ---------------------------------------------------------------------------- #
    #                               Modify Whole PDB                               #
    # ---------------------------------------------------------------------------- #

    def reindex_atom_num(self, starting_idx: int) -> None:
        """
        Reindex all atom numbers in the PDB starting from the specified index.

        Parameters
        ----------
        starting_idx : int
            The starting atom number for reindexing.
        """
        atom_num: int = starting_idx
        for chain in self.chains:
            for residue in chain.residues:
                for record in residue.records:
                    record.atom_number = atom_num
                    atom_num += 1

    # ---------------------------------------------------------------------------- #
    #                               Helper Functions                               #
    # ---------------------------------------------------------------------------- #

    @staticmethod
    def _to_chains(df: DataFrame) -> list[Chain] | None:
        """
        Simple static method to convert dataframe to chain object list.

        Parameters
        ----------
        df : DataFrame
            Target Dataframe to access chains.

        Returns
        -------
        list[Chain] | None
            Returns a list of chains if the dataframe is not empty, otherwise None.
        """
        if not df.empty:
            return PDB.dataframe_to_chains(df)
        else:
            return None

    # ---------------------------------------------------------------------------- #
    #                                  Data Access                                 #
    # ---------------------------------------------------------------------------- #

    def get_chain_idx(self, chain_id: str) -> int:
        """
        Collect Chain index for chain if it exists in pdb.

        Parameters
        ----------
        chain_id : str
            Target chain ID to look for in PDB.

        Returns
        -------
        int
            index of target chain ID.

        Raises
        ------
        ValueError
            If the chain ID is not found in the pdb.
        """
        idx: int | None = self._chain_mapping.get(chain_id.lower(), None)
        if idx is None:
            raise ValueError(f"Chain with id {chain_id} does not exist in PDB!")
        return idx

    # ---------------------------------------------------------------------------- #
    #                             DataFrame Conversions                            #
    # ---------------------------------------------------------------------------- #

    @staticmethod
    def dataframe_to_chains(atom: DataFrame) -> list[Chain]:
        """
        Convert a DataFrame of atomic data into a hierarchical structure of Chain objects.
        This function parses atomic records from a DataFrame and organizes them into
        a nested hierarchy: chains contain residues, which contain atomic records.

        Parameters
        ----------
        atom : DataFrame
            A pandas DataFrame where each row represents an atomic record with columns for:
            - record_name: Type of PDB record
            - atom_number: Unique atom identifier
            - atom_name: Name of the atom
            - residue_name: Name of the parent residue
            - chain_id: Chain identifier
            - residue_number: Sequence number of the residue
            - x_coord, y_coord, z_coord: Cartesian coordinates
            - occupancy: Occupancy value
            - b_factor: Temperature factor
            - segment_id: Segment identifier
            - element_symbol: Chemical element symbol
            - charge: Atomic charge
            - line_idx: Original line index in PDB file
            - Additional blank columns for PDB format compliance

        Returns
        -------
        list[Chain]
            A list of Chain objects, each containing organized residues
            and their associated atomic records in hierarchical order.
        """

        chains: list[Chain] = []
        residues: list[Residue] = []
        records: list[AtomRecord] = []

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
                        AtomRecord(
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
                    Residue(
                        records, starting_idx, residue_name, int(res_id), str(chain_id)
                    )
                )
                records = []

            chains.append(Chain(residues, str(chain_id)))
            residues = []

        return chains

    @staticmethod
    def chains_to_dataframe(chains: list[Chain], starting_idx: int = 0):
        """
        Convert a list of Chain objects to a pandas DataFrame.
        Iterates through all chains, residues, and records to extract atomic information
        and construct a DataFrame with columns representing PDB ATOM record fields.

        Parameters
        ----------
        chains : list[Chain]
            A list of Chain objects containing residues and their atomic records.
        starting_idx : int
            starting index of chain list, zero if unchanged, by default 0.

        Returns
        -------
        DataFrame
            A pandas DataFrame where each row represents an atomic record with columns for:
            - record_name: Type of PDB record
            - atom_number: Unique atom identifier
            - atom_name: Name of the atom
            - residue_name: Name of the parent residue
            - chain_id: Chain identifier
            - residue_number: Sequence number of the residue
            - x_coord, y_coord, z_coord: Cartesian coordinates
            - occupancy: Occupancy value
            - b_factor: Temperature factor
            - segment_id: Segment identifier
            - element_symbol: Chemical element symbol
            - charge: Atomic charge
            - line_idx: Original line index in PDB file
            - Additional blank columns for PDB format compliance
        """
        idx = starting_idx
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
                            "line_idx": idx,
                        }
                    )
                    idx += 1

        return DataFrame(records)

    def _chains_to_full_dataframe(
        self,
    ) -> dict[str, DataFrame]:
        """
        Convert this PDB into a PandasPDB dictionary of dataframes.
        Iterates through all chains, residues, and records to extract atomic information
        and construct an updated dictionary.

        Returns
        -------
        dict
            Dictionary of every PDB component, "ATOM", "HETATM", "ANISOU", and Others after modification.

            Each dataframe is a structure where each row represents an atomic record with columns for:
            - record_name: Type of PDB record
            - atom_number: Unique atom identifier
            - atom_name: Name of the atom
            - residue_name: Name of the parent residue
            - chain_id: Chain identifier
            - residue_number: Sequence number of the residue
            - x_coord, y_coord, z_coord: Cartesian coordinates
            - occupancy: Occupancy value
            - b_factor: Temperature factor
            - segment_id: Segment identifier
            - element_symbol: Chemical element symbol
            - charge: Atomic charge
            - line_idx: Original line index in PDB file
            - Additional blank columns for PDB format compliance
        """

        others: DataFrame = self._others.copy()

        # Extract "TER" records from others dataframe
        ter_records_mask: Series[bool] = others["record_name"].isin(["TER"])
        ter_records = others[ter_records_mask].iterrows()

        # Extract "ENDMDL", "MASTER", and "END" records from others dataframe
        end_records_mask = others["record_name"].isin(["ENDMDL", "MASTER", "END"])
        end_records: Series = others[end_records_mask].iterrows()

        # Collect the Rows at the start of the pdb by ignoring the special records
        pdb_start = others[~(ter_records_mask + end_records_mask)]
        # Collect the final line idx of the start of the pdb
        last_line_idx = pdb_start["line_idx"].max()
        idx: int = last_line_idx + 1

        records: list = []
        for chain in self.chains:
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
                            "line_idx": idx,
                        }
                    )
                    # Increment index as records go up
                    idx += 1
                # end of record
            # end of residue
            # ------------------------------ Update TER Card ----------------------------- #
            last_record = chain[-1][-1]
            row_idx, ter = next(ter_records)

            if not isinstance(last_record, AtomRecord):
                raise ValueError(
                    f"Did not Receive Atom Record when Accessing Last Record of Chain: {chain.chain_id}"
                )
            ter_format = "  {}       {} {}  {}"
            new_ter_entry: str = ter_format.format(
                last_record.atom_number,
                last_record.residue_name,
                last_record.chain_id,
                last_record.residue_number,
            )

            others.loc[row_idx, "line_idx"] = idx
            others.loc[row_idx, "entry"] = new_ter_entry
            # Skip index to account for ter card
            idx += 1
        # end of chain

        # update hetatm indices if it exists
        idx += 1
        if self._hetatm is not None:
            hetatm = PDB.chains_to_dataframe(self._hetatm, idx)
            hetatm_length = sum([chain.total_length() for chain in self._hetatm])
            end_idx = idx + hetatm_length
        else:
            hetatm = DataFrame()
            end_idx = idx

        for row_idx, row in end_records:
            others.loc[row_idx, "line_idx"] = end_idx
            end_idx += 1

        return {
            "ATOM": DataFrame(records),
            "HETATM": hetatm,
            "ANISOU": self._anisou,
            "OTHERS": others,
        }

    # ---------------------------------------------------------------------------- #
    #                                 Magic Methods                                #
    # ---------------------------------------------------------------------------- #

    def __getitem__(self, indices) -> Chain | list[Chain]:
        """
        Chain indexer.

        Parameters
        ----------
        indices : int | slice
            Index or indices of chain.

        Returns
        -------
        Residue | list[Residue]
            Current residue(s) at provided index/indices.
        """
        if isinstance(indices, str):
            return self._chains[self.get_chain_idx(indices)]
        return self._chains[indices]

    def __setitem__(self, index: int | str, value: Chain) -> None:
        """
        Chain setter.

        Parameters
        ----------
        index : int | str
            Index or chain id of chain.
        value : Chain
            Chain object to set.

        Raises
        ------
        ValueError
            If chain id does not exist when using string index.
        """
        if isinstance(index, str):
            idx = self.get_chain_idx(index)
            self._chains[idx] = value
        else:
            self._chains[index] = value

    def __len__(self) -> int:
        """
        Get length of pdb based on number of chains.

        Returns
        -------
        int
            Number of Chains.
        """
        return len(self._chains)
