#!/usr/bin/env python3

import sys


Script, zOTU_counts, zOTU_fasta, output = sys.argv


zOTU_list = []
zOTU_dict = {}


##### Opening zOTU table
COUNTS = open(zOTU_counts, "r")

for line in COUNTS:
    if line.startswith("#"):
        COUNT_headings = line.strip().split()[2:]    ### Keeping the names of libraries
    else:
        LINE = line.strip().split()
        zOTU_list.append(LINE[0])
        zOTU_dict[LINE[0]] = [LINE[1:]]

COUNTS.close()


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


##### Outputting the Expanded Count Table
OUTPUT_TABLE = open(output, "w")

print("OTU_ID", "Sequence", sep = "\t", end = "\t", file = OUTPUT_TABLE)
for item in COUNT_headings:
    print(item, end = "\t", file = OUTPUT_TABLE)
print("", file = OUTPUT_TABLE)

for zOTU in zOTU_list:

    print(zOTU, zOTU_dict[zOTU][1], sep = "\t", end = "\t", file = OUTPUT_TABLE)
    
    for no in zOTU_dict[zOTU][0]:
        print(no, end = "\t", file = OUTPUT_TABLE)
    
    print("", file = OUTPUT_TABLE)

OUTPUT_TABLE.close()
