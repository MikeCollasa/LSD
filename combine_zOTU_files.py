#! /usr/bin/env python3

import sys, re


print("---------- combine_zOTU_files.py v. 1.2, Piotr Łukasik, 21-March-2021 ----------\n")

if len(sys.argv) != 5:
	sys.exit("This script combines info about zOTUs from several files produced by the Symbiosis Evolution Group's amplicon analysis pipeline.\n"
	         "It is intended as a part of the workflow described at https://github.com/symPiotr/amplicon_analysis_pipeline.\n\n"
	         'Usage: ./combine_zOTU_files.py <zOTU_count_table> <zOTU_taxonomy> <zOTU_fasta> <zOTU_OTU_relationships>\n'
	         'For example: ./combine_zOTU_files.py zotu_table.txt zotus.tax zotus.fasta zotu_otu_relationships.txt\n')

Script, zOTU_counts, zOTU_tax, zOTU_fasta, zOTU_OTU_relationships = sys.argv

##### Setting names of output files
Output_table = "zotu_table_expanded.txt"

##### Setting up the key arrays --- LIST for keeping sequences in order, and DICT for managing sequence info
zOTU_list = []
zOTU_dict = {}


##### Opening zOTU table
print("Opening zOTU table.................... ", end="")
COUNTS = open(zOTU_counts, "r")

for line in COUNTS:
    if line.startswith("#"):
        COUNT_headings = line.strip().split()[1:]    ### Keeping the names of libraries
    else:
        LINE = line.strip().split()
        zOTU_list.append(LINE[0])
        zOTU_dict[LINE[0]] = [LINE[1:]]

COUNTS.close()
print("OK!")


##### Adding taxonomy info to DICT
print("Adding taxonomy info.................. ", end="")
TAX = open(zOTU_tax, "r")

for line in TAX:
    LINE = line.strip().split()
    if LINE[0] in zOTU_list:
        if len(LINE) > 1:
            zOTU_dict[LINE[0]].append(LINE[1])
        else:
            zOTU_dict[LINE[0]].append("unassigned")
    else:
        print('FATAL ERROR! Taxonomy file contains zOTUs that are not in zOTU count table! ---', LINE[0])
        sys.exit()

TAX.close()
print("OK!")

##### Adding sequences from the FASTA file to DICT
print("Adding sequences...................... ", end="")
FASTA = open(zOTU_fasta, "r")
Sequence = ''
Seq_heading = FASTA.readline().strip().strip(">")

for line in FASTA:   # Copying the sequence (potentially spread across multiple lines) to a single line
    if line.startswith('>'):    # if the line contains the heading
        if Seq_heading not in zOTU_list and Seq_heading != "":     # EXIT if the previous Seq_heading is not in a list!
            print('FATAL ERROR! Fasta file contains zOTUs that are not in zOTU count table! ---', Seq_heading)
            sys.exit()
            
        zOTU_dict[Seq_heading].append(Sequence) # save the existing Seq_heading and Sequence to a DICT
        Sequence = ''    # clear sequence
        Seq_heading = line.strip().strip(">")  # use the current line as the new heading!

    else:
        Sequence = Sequence + line.strip().upper()

zOTU_dict[Seq_heading].append(Sequence) # Saves the final sequence (Seq_heading and Sequence) to a list

FASTA.close()
print("OK!")

##### Adding zOTU - OTU relationship info to DICT
print("Adding zOTU classification info....... ", end="")
RELS = open(zOTU_OTU_relationships, "r")

for line in RELS:
    LINE = line.strip().split()
    
    zOTU = re.search("^Zotu\d+", LINE[0])[0]
    if zOTU not in zOTU_list:
        print('FATAL ERROR! Relationship file contains zOTUs that are not in zOTU count table! --- ', zOTU)
        sys.exit()
    
    if LINE[1].startswith("otu"):
        zOTU_dict[zOTU].append(LINE[1])
    
    elif  LINE[1] == "noisy_chimera" or LINE[1] == "perfect_chimera" or LINE[1] == "match_chimera" or re.search("Chimera", LINE[2]) != None:
        zOTU_dict[zOTU].append("Chimera")

    elif (LINE[1] == "match" or LINE[1] == "perfect") and re.search("OTU\d+", LINE[2]) != None:
        OTU_ID = re.search("OTU\d+", LINE[2])[0].lower()
        zOTU_dict[zOTU].append(OTU_ID)
        
    else:
        print(f"Relationship file contains a term that I have not considered:\n{line}")
        sys.exit()

RELS.close()
print("OK!")


##### Outputting the Expanded Count Table
print("Outputting data....................... ", end = "")
OUTPUT_TABLE = open(Output_table, "w")

print("OTU_ID", "OTU_assignment", "Taxonomy", "Sequence", "Total", sep = "\t", end = "\t", file = OUTPUT_TABLE)
for item in COUNT_headings:
    print(item, end = "\t", file = OUTPUT_TABLE)
print("", file = OUTPUT_TABLE)

for zOTU in zOTU_list:
    Total = 0
    for no in zOTU_dict[zOTU][0]:
        Total += int(no)
    
    # Terms in DICT: 'Zotu32': [['0', '1', '100'], 'd:Bacteria(1.00)...', 'TACGT...', 'otu8']
    # I want to export: "OTU_ID", "OTU_assignment"[3], "Taxonomy"[1], "Sequence"[2], "Total"
    print(zOTU, zOTU_dict[zOTU][3], zOTU_dict[zOTU][1], zOTU_dict[zOTU][2], Total, sep = "\t", end = "\t", file = OUTPUT_TABLE)
    
    for no in zOTU_dict[zOTU][0]:
        print(no, end = "\t", file = OUTPUT_TABLE)
    
    print("", file = OUTPUT_TABLE)

OUTPUT_TABLE.close()

print("OK!\n")
print(f"Script executed successfully. Output --- {Output_table}\nEnjoy :)")
