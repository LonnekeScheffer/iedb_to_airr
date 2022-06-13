'''
This script converts IEDB receptor csv files to AIRR tsv files.
The user can either choose to use data from the 'Calculated' or 'Curated' columns of the IEDB file.

Usage:
    python3 iedb_to_airr.py <input_file_path> <output_file_path> {Curated|Calculated}

'''

import sys
import pandas as pd
import airr


def get_single_chain_column(df, chain_nr, mode):
    chain_df = df.filter(regex=f"cell_id|^Receptor ID$|^Chain {chain_nr}|^{mode} Chain {chain_nr}", axis=1)
    chain_df.columns = chain_df.columns.str.replace(f"{mode} Chain {chain_nr} ", "")
    chain_df.columns = chain_df.columns.str.replace(f"Chain {chain_nr} ", "")

    return chain_df

def check_mode(mode):
    assert mode in ["Curated", "Calculated"], f"mode must either be 'Curated' or 'Calculated', found '{mode}'"

def get_column_mapping(mode):
    column_mapping = {"Receptor ID": "iedb_receptor_id",
                      "Type": "locus",
                      "Full Sequence": "sequence_aa",
                      "Nucleotide": "sequence",
                      "V Gene": "v_call",
                      "D Gene": "d_call",
                      "J Gene": "j_call"}

    for cdr_nr in [1, 2, 3]:
        column_mapping = {**column_mapping, **{f"CDR{cdr_nr} {mode}": f"cdr{cdr_nr}_aa",
                                               f"CDR{cdr_nr} Start {mode}": f"cdr{cdr_nr}_start",
                                               f"CDR{cdr_nr} End {mode}": f"cdr{cdr_nr}_end"}}

    return column_mapping

def columns_to_airr(chain_df, mode):
    column_mapping = get_column_mapping(mode)
    chain_df.rename(columns=column_mapping, inplace=True)
    chain_df = chain_df[list(column_mapping.values()) + ["cell_id"]]

    return chain_df

def get_locus_from_v_call(v_call):
    if pd.notna(v_call):
        if v_call.startswith("IGK"):
            return "IGK"
        else:
            return "IGL"

def map_locus(row):
    mapping = {"alpha": "TRA",
               "beta": "TRB",
               "gamma": "TRG",
               "delta": "TRD",
               "heavy": "IGH"}

    if row.locus == "light":
        return get_locus_from_v_call(row.v_call)
    else:
        return mapping[row.locus]

def is_junction(cdr3):
    return cdr3[0] == "C" and cdr3[-1] in ["Y", "F"]

def get_junction_if_present(cdr3):
    if pd.notna(cdr3) and is_junction(cdr3):
        return cdr3

def trim_cdr3_if_junction(cdr3):
    if pd.notna(cdr3) and is_junction(cdr3):
        return cdr3[1: -1]
    else:
        return cdr3

def number_as_string(position):
    if position is not None:
        if float(position).is_integer():
            return str(int(position))
        else:
            return str(position)

def format_cdr_positions(chain_df):
    for cdr_nr in [1,2,3]:
        for pos in ["start", "end"]:
            chain_df[f"cdr{cdr_nr}_{pos}"] = chain_df[f"cdr{cdr_nr}_{pos}"].apply(number_as_string)

def get_used_iedb_columns(mode):
    used_cols = ["Receptor ID"]
    for chain_nr in [1,2]:
        used_cols += [f"Chain {chain_nr} Type", f"Chain {chain_nr} Nucleotide", f"Chain {chain_nr} Full Sequence"]
        used_cols += [f"{mode} Chain {chain_nr} {gene} Gene" for gene in ["V", "D", "J"]]
        used_cols += [f"Chain {chain_nr} CDR{cdr_nr} {mode}" for cdr_nr in [1, 2, 3]]
        used_cols += [f"Chain {chain_nr} CDR{cdr_nr} Start {mode}" for cdr_nr in [1, 2, 3]]
        used_cols += [f"Chain {chain_nr} CDR{cdr_nr} End {mode}" for cdr_nr in [1, 2, 3]]

    return used_cols

def iedb_to_airr(df, mode):
    df.drop_duplicates(inplace=True, subset=get_used_iedb_columns(mode))

    df["cell_id"] = list(range(1, len(df)+1))

    df = pd.concat([get_single_chain_column(df, 1, mode),
                    get_single_chain_column(df, 2, mode)],
                   axis=0)

    df.dropna(subset=["Type"], inplace=True)
    df = columns_to_airr(df, mode)

    df["locus"] = df.apply(map_locus, axis=1)
    df["junction_aa"] = df["cdr3_aa"].apply(get_junction_if_present)
    df["cdr3_aa"] = df["cdr3_aa"].apply(trim_cdr3_if_junction)
    df["productive"] = True
    df = df.astype(object).where(df.notna(), None)
    format_cdr_positions(df)

    df.sort_values(by="cell_id", inplace=True)

    return df

def main(input_file_path, output_file_path, mode):
    check_mode(mode)

    df = iedb_to_airr(pd.read_csv(input_file_path), mode)
    airr.dump_rearrangement(df, output_file_path)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3].title())
