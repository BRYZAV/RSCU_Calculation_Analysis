# RSCU_Calculation_Analysis

The files located in this repository are the foundation of calculating and interpreting the RSCU values for the paper (insert here)

Below is a descripution of each file and its purpose (Each file is explained in sequential order to determine which files need to be ran first):

**codonTable_script.sh** - This shell script uses emboss to create codontables for each coding sequence from the CDS yeast files and appends it to a single file of all its codon tables.

**csv_conversion.sh** - This is another shell script that just converts the codon tables into a csv file for further analysis (RSCU_analysis.py).

**RSCU_analysis.py** - This python file is a script used to calculate the RSCU values for every coding sequence in every species. This accounts certain yeast species codon reassignments (CUG-Ser1 clade, CUG-Ser2 clade, and CUG-Ala clade). The files required for this script are the codon table csv files and the 1154yeasts_21outrgoups_info_20220408.csv profile to calculate the genome wide RSCU values.




