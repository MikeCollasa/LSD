#! /usr/bin/env python3

import sys, os
if len(sys.argv) != 4:
    sys.exit("""ERROR! CHECK YOUR INPUT PARAMETERS!
Please provide:
1) sample list with information about your libraries created in following manner:
Sample_name Sample_name_R1.fastq Sample_name_R2.fastq
2) path to the directory with R1 and R2 fiels for all the amplicon libraries that you want to analyse e.g.:
/home/Data/For/Nature/Publication/)
3) type of data, please indicate if you are going to analyse 16SV4 (Bacterial) or COI.""")
Script, sample_list, path_to_your_raw_data, type_of_data = sys.argv

SAMPLE_LIST = open(sample_list, "r")
input = os.listdir(path_to_your_raw_data)


for line in SAMPLE_LIST:
    LINE = line.strip().split()
    if type_of_data == "COI":
        os.system("pear -f %s -r %s -o %s -v 15 -n 400 -m 470 -q 30 -j 20" % (LINE[1], LINE[2], LINE[0]))
    elif type_of_data == "16SV4":
        os.system("pear -f %s -r %s -o %s -v 15 -n 250 -m 400 -q 30 -j 20" % (LINE[1], LINE[2], LINE[0]))
### Removing unassebles and discarded sequences along with renaming assembled ones:
os.system("rm *unassembled* *discarded*")
os.system("rename -f 's/.assembled//' *fastq")
os.system ("mkdir reads && mv *_R?.fastq reads/")

### Fastq to Fasta formatting:
os.system("""for file in *.fastq; do
    SampleName=`basename $file .fastq`
    vsearch -fastq_filter $SampleName.fastq -fastaout $SampleName.fasta -relabel "$SampleName"._ -fasta_width 0
done""")

os.system("mkdir fastq && mv *.fastq fastq/")

### Trimming primers and size filtering:
if type_of_data == "COI":
    os.system("""for file in *.fasta; do
        SampleName=`basename $file .fasta`
        egrep -B 1 "CC[ACT]GA[TC]AT[AG]GC[ACT]TT[TC]CC[ACT]CG.{400,430}TG[AG]TT[TC]TT[TC]GG[ACGT]CA[TC]CC[ACT]G" $SampleName.fasta | egrep -v "^\-\-$" |
        sed -E "s/.*CC[ACT]GA[TC]AT[AG]GC[ACT]TT[TC]CC[ACT]CG//; s/TG[AG]TT[TC]TT[TC]GG[ACGT]CA[TC]CC[ACT]G.*//" > "$SampleName".trimmed.fasta
    done""")
elif type_of_data == "16SV4":
   os.system("""for file in *.fasta; do
        SampleName=`basename $file .fasta`
        egrep -B 1 "GTG[TC]CAGC[CA]GCCGCGGTAA.{250,260}ATTAGA[AT]ACCC[CGT][ACGT]GTAGTCC" all_samples.fasta | egrep -v "^\-\-$" |
        sed -E 's/.*GTG[TC]CAGC[CA]GCCGCGGTAA//; s/ATTAGA[AT]ACCC[ACGT][ACGT]GTAGTCC.*//' > all_samples_trimmed.fasta
    done""")
   
###Deleting untrimmed sequences:
os.system("find . -maxdepth 1 -not -name '*trimmed.fasta' -name '*.fasta' -delete")

os.system("""for file in *.fasta; do
    SampleName=`basename $file .trimmed.fasta`
    vsearch -derep_fulllength $SampleName.trimmed.fasta -output "$SampleName".derep.fasta -sizeout -uc "$SampleName".derep_info.txt
done""")

os.system("mkdir trimmed && mv *trimmed.fasta trimmed/")

###Dereplicating data - picking representative sequences:
os.system("""for file in *derep.fasta; do
    SampleName=`basename $file .derep.fasta`
    vsearch -sortbysize $SampleName.derep.fasta --output $SampleName.sorted.fasta -minsize 2
done""")

os.system("mkdir derep && mv *derep* derep/")

###Denoising:
os.system("""for file in *sorted.fasta; do
    SampleName=`basename $file .sorted.fasta`
    usearch -unoise3 $SampleName.sorted.fasta -zotus $SampleName.zotus.fasta -tabbedout $SampleName.denoising.summary.txt -minsize 1
done""")

os.system("mkdir sorted && mv *sorted.fasta sorted/")
os.system("mkdir denoising_summary && mv *denoising.summary.txt denoising_summary/")

os.system("""for file in *.fasta; do
    SampleName=`basename $file .zotus.fasta`
    usearch -otutab ./trimmed/"$SampleName".trimmed.fasta -zotus $SampleName.zotus.fasta -otutabout "$SampleName"_zotu_table.txt
done""")

###Adding sequence to zOTU_table with add_seq_to_zotu.py:
os.system("""for file in *.fasta; do
    SampleName=`basename $file .zotus.fasta`
    ./add_seq_to_zotu.py "$SampleName"_zotu_table.txt "$SampleName".zotus.fasta "$SampleName"_zotu_table_with_seq.txt
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
os. chdir(path_to_your_raw_data + "/zotu_tables_with_sequences")
cwd = os.getcwd()  # Get the current working directory (cwd)
files = os.listdir(cwd) 

for filename in files:
   with open(filename, "r") as Lib_info:
       get_headings(Lib_info)
       get_lib_name(headings)
       create_library(Lib_info)

os. chdir(path_to_your_raw_data)

# Creating a list with the names of the libraries
libs = []
for k1 in seq_dict.keys():
    libs += [lib for lib in seq_dict[k1] if lib not in libs] #List comprehension
    
### Creating an empty list which will be our final table
data = []
index = 0 #index of a table in data, basically 0 will raw with the first key, 1 - with the second ect
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
        zotu[0] = "Zotu_" + str(counter) 
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
os.system ("mv zotu.fasta ../ && mv all_libraries_zotu_table.txt ../ && cd ..")
###Creating fasta file with all trimmed sequences from all the libraries:
os.system("cd trimmed && cat *.fasta > all_samples_trimmed.fasta && mv all_samples_trimmed.fasta ../ && cd ..")
###OTU picking and chimeras removal using ASV as an input:
os.system("usearch -cluster_otus zotus.fasta -otus otus.fasta -relabel OTU -uparseout zotu_otu_relationships.txt")
os.system("usearch -usearch_global all_samples_trimmed.fasta -db otus.fasta -strand plus -id 0.97 -otutabout otu_table.txt")
### Creating a new fasta file of zOTUs without information about size:
os.system("sed -E 's/;size=[0-9].{0,}//g' zotus.fasta > new_zotus.fasta")
###Assigning taxonomy:
if type_of_data == "COI":
    os.system("""vsearch --sintax new_zotus.fasta -db /mnt/matrix/symbio/db/MIDORI/MIDORI_with_tax_spikeins_endo_RDP.fasta -tabbedout zotus.tax -strand both -sintax_cutoff 0.5
    vsearch --sintax otus.fasta -db /mnt/matrix/symbio/db/MIDORI/MIDORI_with_tax_spikeins_endo_RDP.fasta -tabbedout otus.tax -strand both -sintax_cutoff 0.5""")
elif type_of_data == "16SV4":
    os.system("""search --sintax new_zotus.fasta -db /mnt/matrix/symbio/db/SILVA_138/SILVA_endo_spikeins_RDP.fasta -tabbedout zotus.tax -strand both -sintax_cutoff 0.5
vsearch --sintax otus.fasta -db /mnt/matrix/symbio/db/SILVA_138/SILVA_endo_spikeins_RDP.fasta -tabbedout otus.tax -strand both -sintax_cutoff 0.5""")
###Removing redundant info from out taxonomy files:
os.system("""sed -i 's/[dpcofgs]\://g' zotus.tax
sed -i 's/[dpcofgs]\://g' otus.tax""")
###Putting it all together with Piotr's scripts:
os.system("""./combine_zOTU_files.py all_libraries_zotu_table.txt zotus.tax new_zotus.fasta zotu_otu_relationships.txt
combine_OTU_files.py otu_table.txt otus.tax otus.fasta""")

print("Symbio® Na zdrowie! Salud! Gānbēi (干杯)! Skål!")
