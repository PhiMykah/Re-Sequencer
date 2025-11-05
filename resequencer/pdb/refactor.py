from typing import TYPE_CHECKING

import pandas as pd
from biopandas.pdb import PandasPdb
from pandas import DataFrame, Series

if TYPE_CHECKING:
    from typing import Any


def refactor_column(
    df: DataFrame, column: str, start: "Any | int", length: "Any | int"
):
    if isinstance(start, int) and type(start) is type(length):
        new_range = list(range(start, start + length))

        df = df.reset_index(drop=True)
        if len(new_range) != len(df):
            raise ValueError(
                f"Length mismatch: new_range ({len(new_range)}) "
                f"!= Data Frame ({len(df)})"
            )
        df[column] = new_range
        df[column] = df[column].astype(int)

    else:
        raise ValueError(
            f"Unable to refactor start and length of types "
            f"`{type(start)}` and `{type(length)}`"
        )

    return df


def update_ter(pdb: PandasPdb, new_chain: Series, chain_type: str):
    others = pdb.df["OTHERS"]
    termini_mask = others["record_name"] == "TER"

    # choose a plausible column that holds the TER entry text
    entry_col = (
        "entry"
        if "entry" in others.columns
        else ("text" if "text" in others.columns else others.columns[-1])
    )

    # get the TER entries as strings and check for 'A' in the whitespace-separated tokens of each entry
    termini_entries = others.loc[termini_mask, entry_col].astype(str)
    chain_mask = termini_entries.apply(
        lambda s: chain_type.upper()
        in [token.strip() for token in s.replace(",", " ").split()]
    ).tolist()

    # terminus entries for matching chain (not used directly)

    # If we found a matching TER entry, update it. Otherwise append a new TER entry.
    # Pull values from new_chain Series (fall back to chain_type where appropriate)
    def _get_val(series, *keys, default=None):
        for k in keys:
            if k in series.index:
                return series.loc[k]
        return default

    atom_num = _get_val(new_chain, "atom_number", "atom_num", "serial", "serial_number")
    res_name = _get_val(new_chain, "residue_name", "res_name", "residue")
    chain_id = _get_val(new_chain, "chain_id", default=chain_type)
    res_num = _get_val(new_chain, "residue_number", "res_num", default=None)
    line_idx_val = _get_val(new_chain, "line_idx", default=None)

    # Coerce types to safe defaults
    try:
        atom_num = int(atom_num) if atom_num is not None else ""
    except Exception:
        atom_num = ""

    res_name = str(res_name).strip() if res_name is not None else ""
    chain_id = str(chain_id).strip() if chain_id is not None else chain_type
    try:
        res_num = int(res_num) if res_num is not None else ""
    except Exception:
        res_num = ""

    # Construct a standard TER line similar to PDB formatting. Keep it conservative so it
    # matches likely formats stored in the 'entry'/'text' column.
    # Example: 'TER   1234      ALA A  45'
    atom_field = f"{atom_num:>5}" if isinstance(atom_num, int) else f"{atom_num}"
    resnum_field = f"{res_num:>4}" if isinstance(res_num, int) else f"{res_num}"
    new_entry = f"{atom_field}      {res_name:>3} {chain_id}{resnum_field}"

    # locate the index in the original 'others' DataFrame corresponding to the match
    others_term_df = others.loc[termini_mask]
    matched_indices = others_term_df.index[chain_mask]

    if len(matched_indices) > 0:
        # update the first matching TER (there should be one per chain)
        idx = matched_indices[0]
        others.at[idx, entry_col] = new_entry
        if "line_idx" in others.columns and line_idx_val is not None:
            try:
                others.at[idx, "line_idx"] = int(line_idx_val) + 1
            except Exception:
                others.at[idx, "line_idx"] = line_idx_val + 1
    else:
        # append a new TER row using existing columns where possible
        # initialize new row using column dtypes from existing DataFrame so assigned
        # values use compatible types (numbers -> 0, strings -> empty string)
        new_row = {}
        for col in others.columns:
            try:
                if pd.api.types.is_numeric_dtype(others[col]):
                    new_row[col] = 0
                else:
                    new_row[col] = ""
            except Exception:
                new_row[col] = ""
        new_row["record_name"] = "TER"
        new_row[entry_col] = new_entry
        if "line_idx" in others.columns and line_idx_val is not None:
            try:
                new_row["line_idx"] = int(line_idx_val)
            except Exception:
                new_row["line_idx"] = line_idx_val
        # append preserving index
        others = pd.concat([others, DataFrame([new_row])], ignore_index=True)

    # write back updated others dataframe to pdb
    pdb.df["OTHERS"] = others

    # small debug prints
    print("Updated TER entry for chain", chain_type, "->", new_entry)
