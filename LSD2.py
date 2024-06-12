#! /usr/bin/env python3

import sys
import os
import questionary
import re
import subprocess

print("This is L(ibrary)S(ingle)D(enoise) 2.0 script created in Symbiosis Evolution Research Group® and updated by me\n")

question = questionary.select(
    """\nDo you know what we are doing here, or do you need some clarification?""",
    choices=["I have it all in my pinky finger!", "<°))))><"]
).ask()

if question == "I have it all in my pinky finger!":
    print("Great! Let's proceed.")
elif question == "<°))))><":
    print("""This is the core amplicon analysis workflow.\n
    It uses some basic tools implemented in LINUX (such as rename and basename) or Python and:\n
    PEAR - Paired-End reAd mergeR (https://cme.h-its.org/exelixis/web/software/pear/)\n
    USEARCH (https://www.drive5.com/usearch/)\n
    VSEARCH (https://github.com/torognes/vsearch)\n
    !!!BEFORE RUNNING, CHECK IF YOU HAVE ALL THE DEPENDENCIES INSTALLED!!!\n
    
    First, this script joins F and R (R1 and R2) reads, passing only high-quality ones outputting joined amplicons.\n
    Next, it converts FASTQ to FASTA file, dereplicate (saving memory), and denoise (deleting erronous ones) sequences in each library separately.\n
    Joins all the libraries into one table and assigns all the sequences to taxonomy.\n
    IF you still don't understand, then WELCOME TO THE CLUB!\n
    
    Visit our GitHub page for more instructions: https://github.com/Symbiosis-JU
    """)
    sys.exit()


input_dir = questionary.path("""Please enter the path to the directory containing your FASTQ files.\n
DON'T BE LAZY!!! DON'T EVEN THINK OF PUTTING A DOT THERE! GIVE ME A FULL PATH""").ask()
input_dir = os.path.expanduser(input_dir)

if not os.path.exists(input_dir):
    print("Error: Input directory does not exist.")
    exit()

output_dir = os.path.join(input_dir, "OUTPUT")
os.makedirs(output_dir, exist_ok=True)

os.chdir(input_dir)

type_of_data = questionary.select(
    """\nSelect the type of data.\n 
This will indicate which type of database will be used for the taxonomy assignment.\n
If you choose 'Custom', you'll be asked to provide the path
to your custom database':""",
    choices=["Mitochondrial", "Bacterial", "Fungal", "Custom"]
).ask()

if type_of_data == "Mitochondrial":
    n = 400
    m = 470
elif type_of_data == "Bacterial":
    n = 250
    m = 400
elif type_of_data == "Fungal":
    n = 250
    m = 560
else:
    n = questionary.text("""\nPlease enter the minimum length of your amplicons after joining.\n
For different target genes, we recommend the following values:
16S:250
ITS:250
COI:400:""").ask()
    m = questionary.text("""\nPlease enter the maximum length of your amplicon after joining.\n
For different target genes, we recommend the following values:
16SV1-2 and V4: 400 
16SV3-V4:470 
ITS:560 
COI:470

!!!REMEMBER!!! 
Depending on your sequencing lane specification, you will obtain different read lengths. 
In our Lab, we typically sequence 2x300 bp. That means that after trimming primers (during Multisplit), 
we end up with max ~545 bp joined reads:""").ask()
    n = int(n)
    m = int(m)

q = questionary.text("""\nPlease enter the minimum phread score for your sequences. 
The higher the phread score, the more reads will be discarded, but the results will be more trustworthy. 
We recommend no less than 30:""").ask()
j = questionary.text("\nPlease enter the number of cores you want to allocate in this task:").ask()

n = int(n)
m = int(m)
j = int(j)
q = int(q)


files = os.listdir(input_dir)
r1_files = [f for f in files if '_R1_' in f and f.endswith('.fastq')]
r2_files = [f for f in files if '_R2_' in f and f.endswith('.fastq')]

for r1_file in r1_files:
    sample_name = re.sub(r'_R1_.*\.fastq$', '', r1_file)
    r2_file = f"{sample_name}_R2_{r1_file.split('_R1_')[1]}"
    
    if r2_file in r2_files:
        command = f"pear -f {os.path.join(input_dir, r1_file)} -r {os.path.join(input_dir, r2_file)} -o {os.path.join(output_dir, sample_name)} -v 15 -n {n} -m {m} -q {q} -j {j}"
        print(f"Running command: {command}")
        subprocess.run(command, shell=True)

print("OK!\n")

os.chdir(output_dir)

print("Removing unassembled and discarded sequences..................... ", end="")        
### Removing unassembled and discarded sequences along with renaming assembled ones:
os.system("rm *unassembled* *discarded*")
os.system("rename -f 's/.assembled//' *fastq")
print("OK!")

print("Fastq to Fasta formatting..................... ", end="")
### Fastq to Fasta formatting:
os.system("""for file in *.fastq; do
    SampleName=`basename $file .fastq`
    vsearch -fastq_filter $SampleName.fastq -fastaout $SampleName.fasta -relabel "$SampleName"._ -fasta_width 0
done""")

os.system("mkdir fastq && mv *.fastq fastq/")
print("OK!")   
  
os.system("rename -f 's/.fasta/_raw.fasta/' *fasta ")

print("Dereplicating data..................... ", end="")
os.system("""for file in *.fasta; do
    SampleName=`basename $file _raw.fasta`
    vsearch -derep_fulllength "$SampleName"_raw.fasta -output "$SampleName".derep.fasta -sizeout -uc "$SampleName".derep_info.txt
done""")

os.system("mkdir raw_fasta && mv *_raw.fasta raw_fasta/")

###Dereplicating data - picking representative sequences:
os.system("""for file in *derep.fasta; do
    SampleName=`basename $file .derep.fasta`
    vsearch -sortbysize $SampleName.derep.fasta --output $SampleName.sorted.fasta -minsize 2
done""")

os.system("mkdir derep && mv *derep* derep/")
print("OK!")

print("Denoising..................... ", end="")

###Denoising:

os.system("""for file in *sorted.fasta; do
    SampleName=`basename $file .sorted.fasta`
    usearch -unoise3 $SampleName.sorted.fasta -zotus $SampleName.zotus.fasta -tabbedout $SampleName.denoising.summary.txt -minsize 1
done""")

os.system("mkdir sorted && mv *sorted.fasta sorted/")
os.system("mkdir denoising_summary && mv *denoising.summary.txt denoising_summary/")

os.system(f"""for file in *.fasta; do
    SampleName=$(basename $file .zotus.fasta)
    usearch -otutab ./raw_fasta/"$SampleName"_raw.fasta -zotus $SampleName.zotus.fasta -otutabout "$SampleName"_zotu_table.txt -threads {j}
done""")

print("OK!") 

print("Joinign data from all the libraries into one tabel..................... ", end="")

###Adding sequence to zOTU_table with add_seq_to_zotu.py:
os.system("""for file in *.fasta; do
    SampleName=`basename $file .zotus.fasta`
    /mnt/qnap/users/symbio/software/Informative_indexes_script/add_seq_to_zotu.py "$SampleName"_zotu_table.txt "$SampleName".zotus.fasta "$SampleName"_zotu_table_with_seq.txt
done""")

os.system("mkdir raw_zotu && mv *_zotu_table.txt raw_zotu && mv *zotus.fasta raw_zotu")
os.system("mkdir zotu_tables_with_sequences && mv *zotu_table_with_seq.txt zotu_tables_with_sequences")



###Creating a dictionary with sequences as keys and nested dictionary as values. Nested dictionary uses names of libraries as keys and number of reads as values:

seq_dict = {}

def get_headings(file):
   global headings
   headings = file.readline().strip('\t\n').split('\t')
   
def get_lib_name(heading):
    global lib_name
    for i in range(0, len(headings)):
       lib_name = headings[2]

def create_library(Lib_info):
   global lib_dict
   for line in Lib_info:
        LINE = line.strip('\t\n').split('\t')
        key = LINE[1]
        counts = LINE[-1]
        if not key in seq_dict.keys():
        #Creating nested dictionary with library names as keys and counts as values
            seq_dict[key] = {lib_name : counts}
        else:
            seq_dict[key].update({lib_name : counts})

#Changing working directory to one with our zotu_tables_with_sequences
os.chdir(os.path.join(input_dir, "OUTPUT", "zotu_tables_with_sequences"))
cwd = os.getcwd()  # Get the current working directory (cwd)
files = os.listdir(cwd) 

for filename in files:
   with open(filename, "r") as Lib_info:
       get_headings(Lib_info)
       get_lib_name(headings)
       create_library(Lib_info)

os.chdir(output_dir)

# Creating a list with the names of the libraries
libs = []
for k1 in seq_dict.keys():
    libs += [lib for lib in seq_dict[k1] if lib not in libs] #List comprehension
    
### Creating an empty list which will be our final table
data = []
index = 0 #index of a table in data, basically 0 will be raw with the first key, 1 - with the second ect
for seq in seq_dict.keys():
    data.append(["",seq, 0]) ### append  empty zotu_ID, sequence, total (which is zero at the beginning)
    for lib in libs: #For every library in our list
        data[index].append(0 if lib not in seq_dict[seq].keys() else seq_dict[seq][lib]) #append 0 if library is not in keys of nested dictionary for given sequence
        if lib in seq_dict[seq].keys(): ### summing Total for every library in the keys of sequence
            data[index][2] += int(seq_dict[seq][lib])
    index += 1 #increase index of one and go to the next sequence
 

def byTotal(Total): ###  Function which will indicate what to sort by
    return Total[2]
 
headers = ["#OTU_ID"] + libs ### Creating a list with headers
with open("all_libraries_zotu_table.txt", 'w') as bigFile:
    element = ""
    for header in headers:
        element += header + "\t" 
    bigFile.write(element[:-1] + '\n')
 
    data.sort(key=byTotal, reverse = True) ### sorting our table by Total with decreasing number of reads for each zOTU
    ###appending new zotu name for each sequence starting with the most abundant
    counter = 1
    for zotu in data:
        zotu[0] = "Zotu" + str(counter) 
        counter += 1

    for zotu in data:
        line = "" ###creating an empty string
        newData = [zotu[0]] + zotu[3:] ### adding ZOTU_ID and all the libraries to the list
        for element in newData:
            line += str(element) + "\t"
        bigFile.write(line[:-1] + '\n')
###Creating fasta file with sorted sequences of zOTUs:   
with open ("zotus.fasta", 'w') as fasta:
    for zotu in data:
        seq = ">" + zotu[0] + ";size=" + str(zotu[2]) + '\n' + zotu[1] + '\n'
        fasta.write(seq)
print("OK!") 

print("OTU picking and chimeras assignment..................... ", end="")
###OTU picking and chimeras removal using ASV as an input:
command = f"usearch -cluster_otus zotus.fasta -otus otus.fasta -relabel OTU -uparseout zotu_otu_relationships.txt"
os.system(command)
print("OK!") 

### Creating a new fasta file of zOTUs without information about the size:
os.system("sed -E 's/;size=[0-9].{0,}//g' zotus.fasta > new_zotus.fasta")

print("Assigning taxonomy..................... ", end="")
###Assigning taxonomy:
#Please pay attention to what cutoff value you want to use!

if type_of_data == "Mitochondrial":
    db_path = "/mnt/qnap/users/symbio/software/databases/MIDORI_with_tax_spikeins_endo_RDP.fasta"
elif type_of_data == "Bacterial":
    db_path = "/mnt/qnap/users/symbio/software/databases/SILVA_endo_spikeins_RDP.fasta"
elif type_of_data == "Fungal":
    db_path = "/mnt/qnap/users/symbio/software/databases/utax_reference_dataset_all_10.05.2021.fasta"
elif type_of_data == "Custom":
    db_path = questionary.text("Please enter the path to the custom database (-db):").ask()

os.system(f"""vsearch --sintax new_zotus.fasta -db {db_path} -tabbedout zotus.tax -strand both -sintax_cutoff 0.8 -threads {j}
vsearch --sintax otus.fasta -db {db_path} -tabbedout otus.tax -strand both -sintax_cutoff 0.8 -threads {j}""")

###Removing redundant info from out taxonomy files:
os.system("""sed -i 's/[dpcofgs]\://g' zotus.tax
sed -i 's/[dpcofgs]\://g' otus.tax""")
print("OK!") 

print("Outputting OTU and zOTU Tables..................... ", end="")

##### Setting names of output files
Output_table = "zotu_table_expanded.txt"

##### Setting up the key arrays --- LIST for keeping sequences in order, and DICT for managing sequence info
zOTU_list = []
zOTU_dict = {}


##### Opening zOTU table

COUNTS = open("all_libraries_zotu_table.txt", "r")

for line in COUNTS:
    if line.startswith("#"):
        COUNT_headings = line.strip().split()[1:]    ### Keeping the names of libraries
    else:
        LINE = line.strip().split()
        zOTU_list.append(LINE[0])
        zOTU_dict[LINE[0]] = [LINE[1:]]

COUNTS.close()


##### Adding taxonomy info to DICT

TAX = open("zotus.tax", "r")

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

##### Adding sequences from the FASTA file to DICT
FASTA = open("new_zotus.fasta", "r")
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

##### Adding zOTU - OTU relationship info to DICT

RELS = open("zotu_otu_relationships.txt", "r")

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
        print("Relationship file contains a term that I have not considered")
        sys.exit()

RELS.close()


##### Outputting the Expanded Count Table
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

print("zOTU_Table_expanded has been created!")


### Creating OTU_Table:

OTU = open("zotu_table_expanded.txt", "r")
OTU_TABLE = []
for line in OTU:
    LINE = line.strip().split()
    if line.startswith("OTU_ID"):
        COUNT_headings = line.strip().split()[4:]    ### Keeping the names of libraries
    else:
        OTU_TABLE.append(LINE)   
OTU.close()

otu_dict = {}
for row_no in range(0, len(OTU_TABLE)):
    otu_key = OTU_TABLE[row_no][1]
    if not otu_key in otu_dict.keys():
         otu_dict[otu_key] = OTU_TABLE[row_no][4:]
    else:
        otu_dict[otu_key] = [sum(map(int, i)) for i in list(zip(otu_dict[otu_key], OTU_TABLE[row_no][4:]))]

##### Adding taxonomy info to DICT
TAX = open("otus.tax", "r")
OTU_TAX = []
for line in TAX:
    LINE = line.strip().split()
    OTU_TAX.append(LINE)

### Lowering the #OTU in Taxonomy file:
for list in OTU_TAX:
    list[0] = list[0].lower()       
        

for row_no in range(0, len(OTU_TAX)):   
    if OTU_TAX[row_no][0] in otu_dict.keys():
        if len(OTU_TAX[row_no]) > 1:
            otu_dict[OTU_TAX[row_no][0]].append(OTU_TAX[row_no][1])
        else:
            otu_dict[OTU_TAX[row_no][0]].append("unassigned")
TAX.close()


                
###We are adding 1 to the end of our dictionary to 
for row_no in range(0, len(OTU_TABLE)):
    if otu_dict[OTU_TABLE[row_no][1]][-1] != 1 and OTU_TABLE[row_no][1] in otu_dict.keys():
        otu_dict[OTU_TABLE[row_no][1]].append(OTU_TABLE[row_no][3])
        otu_dict[OTU_TABLE[row_no][1]].append(1)      

                
COUNT_headings.insert(0,"#OTU")
COUNT_headings.insert(1,"Taxonomy")
COUNT_headings.insert(2,"Sequence")
data = []
data.append(COUNT_headings)      
for otu in otu_dict.keys():
    data.append([otu] + [otu_dict[otu][-3]] + [otu_dict[otu][-2]] + otu_dict[otu][:-3])
 


with open("OTU_Table.txt", "w") as bigFile:
    for LINE in data:
        for item in LINE[:-1]:
            print(item, end="\t", file = bigFile)
        for item in LINE [-1:]:
            print(item, end='\n',file = bigFile)
bigFile.close()
print("OTU_Table has been created!")


print("Symbio® Na zdrowie! Salud! Gānbēi (干杯)! Skål!")
