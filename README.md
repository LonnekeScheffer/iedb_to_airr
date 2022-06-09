# IEDB to AIRR
Converter from csv files in IEDB format to AIRR format. 

Notes:
- IEDB files contain both 'Curated' and 'Calculated' columns, which may contain contradicting information 
([documentation](https://help.iedb.org/hc/en-us/articles/4565781888027-What-is-the-difference-between-the-curated-and-calculated-receptor-sequence-information-)).
In this converter, the user must specify which of these two column types to use. Thus, the information
in the AIRR file either originates from the Curated or from the Calculated columns, but the two are 
never combined. 
- IEDB files may contain multiple records sharing the same Receptor ID. These receptors always 
share the same CDR3 sequence, but other fields (e.g., CDR1, CDR2) may differ. In this converter, all
such receptors are exported to the AIRR file as individual records. 
- The custom column 'iedb_receptor_id' was added to the AIRR file which should allow users to look up
additional information about the relevant receptor in the IEDB again if needed, or choose to remove 
(near-)duplicate entries sharing the same iedb_receptor_id.
- All IEDB columns that did not have an AIRR format counterpart were omitted from the AIRR output file. 
- In AIRR format, the only column that links the two individual chains that make up a receptor 
is the 'cell_id' column. This script adds a unique cell_id to each individual record in the IEDB file, 
meaning that there may be multiple receptors sharing the same iedb_receptor_id (Receptor ID in IEDB 
format) but have a different cell_id, and for a given cell_id there will always be either 1 (single chain)
or 2 (paired chain) rearrangements present in the AIRR output file. 
- The IEDB CDR3 records sometimes include the conserved leading C and trailing Y/F positions, and sometimes
omit these positions. When converting to AIRR format, the CDR3 sequence without conserved positions is
stored in the 'cdr3_aa' field, and if conserved positions were present, this version of the CDR3 is 
additionally stored in the 'cdr3_junction' field.
- The IEDB 'light' chain type can be mapped to either the IGL and IGK locus in AIRR format. For light 
chains, the locus is determined based on the V gene if available. If the V gene is absent, the locus
field for this light chain is left empty in the AIRR output format. For other chains, the IEDB 'Chain Type'
columns are used.


