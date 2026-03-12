from copy import deepcopy

from .residue import Residue


class Chain:
    """
    Class representation of Chain in PDB and PandasPDB.

    Attributes
    ----------
    residues : list[Residue]
        List of Residues in Chain.
    chain_id : str
        Name or id of current chain.
    starting_idx : int
        Starting index of first residue.
    """

    # ---------------------------------------------------------------------------- #
    #                        Initialization and Constructor                        #
    # ---------------------------------------------------------------------------- #
    def __init__(
        self, residues: list[Residue], chain_id: str, starting_idx: int = -1
    ) -> None:
        self._residues = []
        self._chain_id = chain_id
        self._starting_idx = starting_idx
        self._set_residues(residues, chain_id, starting_idx)

    # ---------------------------------------------------------------------------- #
    #                              Getters and Setters                             #
    # ---------------------------------------------------------------------------- #

    # --------------------------------- residues --------------------------------- #
    @property
    def residues(self) -> list[Residue]:
        return self._residues

    @residues.setter
    def residues(self, value: list[Residue]) -> None:
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

    # ------------------------------- starting_res ------------------------------- #
    @property
    def starting_residue(self) -> int:
        return self._starting_residue

    @starting_residue.setter
    def starting_residue(self, value: int) -> None:
        self._starting_residue = value
        self._reorder_residue_numbers(value)

    # ---------------------------------------------------------------------------- #
    #                                Modifying Chain                               #
    # ---------------------------------------------------------------------------- #

    def append(self, residue: Residue) -> None:
        """
        Append a new residue to the chain and re-index residues.

        Parameters
        ----------
        residue : Residue
            The residue to append to this chain.
        """
        residue.chain_id = self._chain_id
        self._residues.append(residue)
        self._reorder_residue_numbers(self._residues[0].residue_number)

    def insert(self, residue: Residue, idx: int) -> None:
        """
        Insert a new residue to the chain at the target index and re-index residues.

        Parameters
        ----------
        residue : Residue
            The residue to append to this chain.
        idx : int
            Position to add new residue.
        """
        residue.chain_id = self._chain_id
        self._residues.insert(idx, residue)
        self._reorder_residue_numbers(self._residues[0].residue_number)

    def remove(self, residue_number: int, reorder: bool = True) -> None:
        """Remove a residue from the chain by residue number and re-index if necessary.

        Parameters
        ----------
        residue_number : int
            Residue with given residue number to remove.
        reorder : bool, optional
            Whether or not to reorder residue numbers, by default True.
        """
        target_idx: int | None = None
        for idx, res in enumerate(self.residues):
            if res.residue_number == residue_number:
                target_idx = idx

        if target_idx is not None:
            self.remove_at(target_idx, reorder)
        else:
            raise ValueError("Residue Number '{residue_number}', not in chain!")

    def remove_at(self, idx: int, reorder: bool = True) -> None:
        """
        Remove a residue from the chain at index and re-index if necessary.

        Parameters
        ----------
        idx : int
            Position to remove residue.
        reorder: bool
            Whether or not to reorder residue numbers, by default True.
        """

        self._residues.pop(idx)
        if reorder:
            self._reorder_residue_numbers(self._residues[0].residue_number)

    # ---------------------------------------------------------------------------- #
    #                              Modifying Residues                              #
    # ---------------------------------------------------------------------------- #

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
        return sum(len(res) for res in self._residues)

    # ---------------------------------------------------------------------------- #
    #                               Helper Functions                               #
    # ---------------------------------------------------------------------------- #

    def _reorder_residue_numbers(self, starting_number: int) -> None:
        """
        Re-index residues in chain by starting number.

        Parameters
        ----------
        starting_number : int
            Residue number of first residue in chain.
        """
        for idx, residue in enumerate(self._residues):
            residue.residue_number = starting_number + idx
        self._starting_residue = starting_number

    def _set_residues(
        self, residues: list[Residue], chain_id: str, starting_idx: int
    ) -> None:
        """
        Set residue list for Chain and ensure proper indexing for each line_idx.

        Parameters
        ----------
        residues : list[Residue]
            List of Residues in Chain.
        chain_id : str
            Name or id of current chain.
        starting_idx : int
            Starting index of first residue.
        """
        if not residues:
            return
        self._residues: list[Residue] = []
        for idx, residue in enumerate(residues):
            if starting_idx == -1 and idx == 0:
                self._starting_idx = residue.starting_idx
                starting_idx = residue.starting_idx
            residue.chain_id = chain_id
            # residue.residue_number = starting_idx + idx
            self._residues.append(residue)
        self._starting_residue: int = self._residues[0].residue_number

    # ---------------------------------------------------------------------------- #
    #                                 Magic Methods                                #
    # ---------------------------------------------------------------------------- #
    def __deepcopy__(self, memo):
        """
        Create a deep copy of the chain.

        Returns
        -------
        Chain
            Deep copy of current chain.
        """
        new_residues = [deepcopy(res, memo) for res in self._residues]
        return Chain(new_residues, self._chain_id, self._starting_idx)

    def __add__(self, other: "Chain") -> "Chain":
        """
        Add two chains together to create a new chain.

        Parameters
        ----------
        other : Chain
            Chain to add to this chain.

        Returns
        -------
        Chain
            New chain containing residues from both chains.
        """
        combined_residues = self._residues + other._residues
        new_chain = Chain(combined_residues, self._chain_id, self._starting_idx)
        new_chain._reorder_residue_numbers(self._starting_residue)
        return new_chain

    def __len__(self) -> int:
        """
        Get length of chain based on number of residues.

        Returns
        -------
        int
            Number of residues.
        """
        return len(self._residues)

    def __iter__(self):
        """Iterator for Chain.

        Returns
        -------
        self
        """
        self.iteraton: int = 0
        return self

    def __next__(self) -> Residue:
        """Next iteration for Chain iterator.

        Returns
        -------
        Residue
            Next residue in iterator.

        Raises
        ------
        StopIteration
            When iteration has completed.
        """
        if self.iteraton < len(self._residues):
            res = self._residues[self.iteraton]
            self.iteraton += 1
            return res
        else:
            raise StopIteration

    def __getitem__(self, indices) -> Residue | list[Residue]:
        """
        Chain indexer.

        Parameters
        ----------
        indices : int | slice
            Index or indices of residue.

        Returns
        -------
        Residue | list[Residue]
            Current residue(s) at provided index/indices.
        """
        return self._residues[indices]
