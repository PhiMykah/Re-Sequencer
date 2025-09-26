from pandas import DataFrame

from .base import (
    NucleotideBase,
    deleted_atoms,
    purine_to_pyrimidine,
    pyrimidine_to_purine,
    PURINE,
    PYRIMIDINE,
)


class Residue:
    """
    A class representing a residue extracted from an atom PDB DataFrame, providing methods to access and modify its properties.

    Attributes
    ----------
        _residue_number (int)
            The residue number in the atom Dataframe.
        _residue (DataFrame)
            The DataFrame containing the residue's atoms.
        _residue_length (int)
            The number of atoms in the residue.
        _residue_name (str or None)
            The name of the residue.

    Properties
    ----------
        num (int)
            The residue number. Can be set to update the residue number in the DataFrame.

    Methods
    -------
        to_dataframe() -> DataFrame:
            Returns the DataFrame representing the residue.
        substitute(old_base: str, new_base: str, update_residue: bool = True) -> DataFrame:
            Substitutes the base in the residue with a new base. If update_residue is True, modifies the current residue in place.
            Otherwise, returns a new DataFrame with the substitution.
        __repr__() -> str:
            Returns the string representation of the residue DataFrame.
        _repr_html_() -> str | None:
            Returns the HTML representation of the residue DataFrame for display in Jupyter notebooks.
    """

    def __init__(self, df: DataFrame, residue_number: int):
        """
        Residue Constructor

        Parameters
        ----------
        df (DataFrame)
            The DataFrame containing all atom residues.
        residue_number (int)
            The residue number to extract from the DataFrame and build the Residue.
        """
        self._residue_number: int = residue_number
        self._residue: DataFrame = _get_residue_from_dataframe(df, self._residue_number)
        self._residue_length: int = self._residue.shape[0]
        self._residue_name = (
            self._residue["residue_name"].iloc[0]
            if "residue_name" in self._residue.columns
            else None
        )

    # ---------------------------------------------------------------------------- #
    #                              Getters and Setters                             #
    # ---------------------------------------------------------------------------- #

    @property
    def num(self) -> int:
        """
        Getter property for residue number.

        Returns
        -------
        int
            residue number corresponding to current Residue object.
        """
        return self._residue_number

    @num.setter
    def num(self, value: int) -> None:
        """
        Set residue number of current object.

        Parameters
        ----------
        value (int)
            New residue number for Residue object.
        """
        self._update_number(value)

    def to_dataframe(self) -> DataFrame:
        """
        Returns the DataFrame representing the residue.

        Returns
        -------
        DataFrame
            Residue DataFrame representation.
        """
        return self._residue

    # ---------------------------------------------------------------------------- #
    #                                Modify Residue                                #
    # ---------------------------------------------------------------------------- #

    def substitute(
        self,
        old_base: "str",
        new_base: "str",
        update_residue: bool = True,
    ) -> DataFrame:
        """
        Substitutes the base in the residue with a new base. If update_residue is True, modifies the current residue in place.

        Otherwise, returns a new DataFrame with the substitution.

        Parameters
        ----------
        old_base (str)
            Current base to resequence from.
        new_base (str)
            New base for current Residue object.
        update_residue (bool, optional)
            Modify the current residue object if True otherwise return new DataFrame, by default True.

        Returns
        -------
        DataFrame
            Post-substitution residue DataFrame.
        """
        if update_residue:
            # Remove the base atoms before conversion
            _delete_base_atoms(self._residue, old_base, update_residue)
            # Change remaining base atoms to new base
            _change_base_atoms(self._residue, old_base, new_base, update_residue)
            # Update the residue name
            self._update_residue_name(old_base, new_base)
            return self._residue
        else:
            # Remove the base atoms before conversion
            trimmed_residue = _delete_base_atoms(
                self._residue, old_base, update_residue
            )
            # Change remaining base atoms to new base
            new_residue = _change_base_atoms(
                trimmed_residue, old_base, new_base, update_residue
            )
            # Update the residue name
            new_residue.loc[new_residue["residue_name"] == old_base, "residue_name"] = (
                new_base
            )
            return new_residue

    # ---------------------------------------------------------------------------- #
    #                               Helper Functions                               #
    # ---------------------------------------------------------------------------- #

    def _update_number(self, value: int) -> None:
        """
        Update Residue number internally and modify the DataFrame to reflect the change.

        Parameters
        ----------
        value (int)
            New residue number for Residue object.
        """
        self._residue_number = value
        self._residue.loc[:, "residue_number"] = value

    def _update_residue_name(self, old_base: str, new_base: str) -> None:
        """
        Update Residue base name and modify the DataFrame to reflect change.

        Parameters
        ----------
        old_base (str)
            Current base value
        new_base (str)
            New base value
        """
        self._residue.loc[self._residue["residue_name"] == old_base, "residue_name"] = (
            new_base
        )
        self._residue_name = new_base

    # ---------------------------------------------------------------------------- #
    #                                Magic Functions                               #
    # ---------------------------------------------------------------------------- #

    def __repr__(self) -> str:
        """
        Returns the string representation of the residue DataFrame
        """
        return self._residue.__repr__()

    def _repr_html_(self) -> str | None:
        """
        Returns the HTML representation of the residue DataFrame for display in Jupyter notebooks.
        """
        return self._residue._repr_html_()  # type: ignore


# ---------------------------------------------------------------------------- #
#                            Global Helper Functions                           #
# ---------------------------------------------------------------------------- #


def _get_residue_from_dataframe(df: DataFrame, residue_number: int) -> DataFrame:
    """
    Extract residue sub-region from entire DataFrame

    Parameters
    ----------
    df (DataFrame)
        Full PDB atom DataFrame.
    residue_number (int)
        The residue number to extract from the DataFrame and build the Residue.

    Returns
    -------
    DataFrame
        The DataFrame containing the residue's atoms.
    """
    return df[df["residue_number"] == residue_number].copy()


def _delete_base_atoms(
    df: DataFrame, base: "NucleotideBase | str | int", update_df: bool = True
) -> DataFrame:
    """
    Removes atoms corresponding to a specified nucleotide base from a DataFrame.

    Parameters
    ----------
        df (DataFrame)
            The input DataFrame containing atom information, with an "atom_name" column.
        base (NucleotideBase | str | int)
            The nucleotide base to remove atoms for. Can be a NucleotideBase enum, string, or integer.
        update_df (bool, optional)
            If True, updates and returns the input DataFrame with the base atoms removed.

            If False, returns a new DataFrame with the base atoms removed. Default is True.

    Returns
    -------
        DataFrame
            The DataFrame with atoms of the specified nucleotide base removed.

    Raises
    ------
        ValueError
            If the provided nucleotide base is not recognized or out of range.
    """
    nucleotide_base: int = _base_to_int(base)

    if nucleotide_base not in deleted_atoms.keys():
        raise ValueError(
            f"Nucleotide base out of range; got {nucleotide_base}",
        )

    # Remove rows where "atom_name" matches one of the base atoms
    trimmed_df: DataFrame = df[~df["atom_name"].isin(deleted_atoms[nucleotide_base])]

    if update_df:
        df = trimmed_df
        return df

    return trimmed_df


def _change_base_atoms(
    df: DataFrame,
    old_base: "NucleotideBase | str | int",
    new_base: "NucleotideBase | str | int",
    update_df: bool = True,
) -> DataFrame:
    """
    Change the base atoms in a nucleotide DataFrame from one base type to another.
    This function swaps the atomic structure of nucleotides in the given DataFrame
    from an old base to a new base, handling conversions between purines and pyrimidines.
    If the conversion is not between purine and pyrimidine (or vice versa), the original
    DataFrame is returned unchanged.

    Parameters
    ----------
        df (DataFrame)
            The DataFrame containing nucleotide atom information.
        old_base (NucleotideBase | str | int)
            The original nucleotide base to be changed.
        new_base (NucleotideBase | str | int)
            The target nucleotide base to change to.
        update_df (bool, optional)
            If True, updates the original DataFrame with the new base atoms.
            Defaults to True.
    Returns
    -------
        DataFrame
            The DataFrame with updated base atoms if a conversion occurred, otherwise the original DataFrame.
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

    # Check if conversion is between different compounds
    if input_base in PURINE and output_base in PYRIMIDINE:
        new_df: DataFrame = _swap_base_atoms(df, purine_to_pyrimidine, False)
    elif input_base in PYRIMIDINE and output_base in PURINE:
        new_df: DataFrame = _swap_base_atoms(df, pyrimidine_to_purine, False)
    else:
        return df

    if update_df:
        df = new_df

    return new_df


def _swap_base_atoms(
    df: DataFrame, mapping: dict[str, str], update_df: bool = True
) -> DataFrame:
    """
    Swap atom names in a DataFrame according to a provided mapping.
    This function replaces values in the 'atom_name' column of the given DataFrame
    based on the specified mapping dictionary. The replacement can be performed
    in-place or on a copy of the DataFrame.

    Parameters
    ----------
        df (DataFrame)
            The input DataFrame containing an 'atom_name' column.
        mapping (dict[str, str])
            A dictionary mapping original atom names to new atom names.
        update_df (bool, optional)
            If True, modify the input DataFrame in-place.

            If False, perform the operation on a copy and return it. Defaults to True.

    Returns
    -------
        DataFrame
            The DataFrame with updated atom names, either modified in-place or as a new copy.
    """
    if update_df:
        df.loc[df["atom_name"].isin(mapping.keys()), "atom_name"] = df.loc[
            df["atom_name"].isin(mapping.keys()), "atom_name"
        ].map(mapping)
        return df
    else:
        new_df: DataFrame = df.copy()
        new_df.loc[new_df["atom_name"].isin(mapping.keys()), "atom_name"] = new_df.loc[
            new_df["atom_name"].isin(mapping.keys()), "atom_name"
        ].map(mapping)
        return new_df


def _base_to_int(base: "NucleotideBase | str | int") -> int:
    """
    Converts a nucleotide base representation to its corresponding integer value.

    Parameters
    ----------
        base (NucleotideBase | str | int)
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
