# python scripts for terra workflows

This repository contains the custom python scripts used within our wdl workflows on Terra for processing COVID-19 sequencing data. 

# preprocess_nanopore_illumina_python_script.py
This script is used both with our preprocess-illumina and pre-process-nanopore genome assembly workflows on terra. This script uses the assemblied genome to determine the percent of non ambigous bases (i.e. percent coverage) and saves the percent coverage as a csv file


# classifyCOVIDlineages_nextclade_json_parser.py
This script is used for our classifyCOVIDlineages workflow. It parses the clade and variant information from the json file output from nextclade. The clade output and vairant output are saved as two seperate csvs.

#classifyCOVIDlineages_python_script.py
This script is used in our classifyCOVIDlienages workflow. It produces the output summary tables saved as seperate csv files  including 1) sequence assemlby metrics (percent coverage, depth etc), 2)wgs horizon report (for our internal use to parse results into our LIMS system), and 3) sequencing results (which includes the sequence assemlby metrics (percent coverage, depth), viral lineage and clade, list of mutations occuring within the RBD and Polybasic Cleavage Site of the Spike protein as well as additional mutations of interest.

 
