#! /usr/bin/env python3

import sys, os
if len(sys.argv) != 2:
    sys.exit("""ERROR! CHECK YOUR INPUT PARAMETERS!
Please provide:
1) sample list with information about your libraries created in following manner:
Sample_name Sample_name_R1.fastq Sample_name_R2.fastq
2) path to the directory with R1 and R2 fiels for all the amplicon libraries that you want to analyse e.g.:
/home/Data/For/Nature/Publication/""")
Script, sample_list, path_to_your_raw_data = sys.argv

SAMPLE_LIST = open(sample_list, "r")
input = os.listdir(path_to_your_raw_data)


for line in SAMPLE_LIST:
    LINE = line.strip().split()
    os.system("pear -f %s -r %s -o %s -v 15 -n 400 -m 470 -q 30 -j 20" % (LINE[1], LINE[2], LINE[0]))

os.system("rm *unassembled* *discarded*")
os.system("rename -f 's/.assembled//' *fastq")
os.system ("mkdir reads && mv *_R?.fastq reads/")

os.system("""for file in *.fastq; do
    SampleName=`basename $file .fastq`
    vsearch -fastq_filter $SampleName.fastq -fastaout $SampleName.fasta -relabel "$SampleName"._ -fasta_width 0
done""")

os.system("mkdir fastq && mv *.fastq fastq/")


os.system("""for file in *.fasta; do
    SampleName=`basename $file .fasta`
    egrep -B 1 "CC[ACT]GA[TC]AT[AG]GC[ACT]TT[TC]CC[ACT]CG.{400,430}TG[AG]TT[TC]TT[TC]GG[ACGT]CA[TC]CC[ACT]G" $SampleName.fasta | egrep -v "^\-\-$" |
    sed -E "s/.*CC[ACT]GA[TC]AT[AG]GC[ACT]TT[TC]CC[ACT]CG//; s/TG[AG]TT[TC]TT[TC]GG[ACGT]CA[TC]CC[ACT]G.*//" > "$SampleName".trimmed.fasta
done""")

os.system("find . -maxdepth 1 -not -name '*trimmed.fasta' -name '*.fasta' -delete")

os.system("""for file in *.fasta; do
    SampleName=`basename $file .trimmed.fasta`
    vsearch -derep_fulllength $SampleName.trimmed.fasta -output "$SampleName".derep.fasta -sizeout -uc "$SampleName".derep_info.txt
done""")

os.system("mkdir trimmed && mv *trimmed.fasta trimmed/")

os.system("""for file in *derep.fasta; do
    SampleName=`basename $file .derep.fasta`
    vsearch -sortbysize $SampleName.derep.fasta --output $SampleName.sorted.fasta -minsize 2
done""")

os.system("mkdir derep && mv *derep.* derep/")


os.system("""for file in *sorted.fasta; do
    SampleName=`basename $file .sorted.fasta`
    usearch -unoise3 $SampleName.sorted.fasta -zotus $SampleName.zotus.fasta -tabbedout $SampleName.denoising.summary.txt -minsize 1
done""")

os.system("mkdir sorted && mv *sorted.fasta sorted/")
os.system("mkdir denoising_summary && mv *denoising.summary.txt denoising_summary/")
os.system("""for file in *.fasta; do
    SampleName=`basename $file .zotus.fasta`
    usearch -otutab /home/michal.kolasa/macrosteles_laevis_COI/trimmed/"$SampleName".trimmed.fasta -zotus $SampleName.zotus.fasta -otutabout "$SampleName"_zotu_table.txt
done""")
