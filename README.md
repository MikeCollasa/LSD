# LSD
LSD: Library Single Denoising


Here we provide the details of the updated amplicon data analysis workflow developed in February 2021 by [Michał Kolasa](https://symbio.eko.uj.edu.pl/staff/michal-kolasa), [Sandra Ahlén Mulio](https://symbio.eko.uj.edu.pl/phd-students/sandra-ahlen-mulio) and [Catalina Valdivia](https://symbio.eko.uj.edu.pl/phd-students/catalina-valdivia). Our pipeline is heavely based on [one](https://github.com/symPiotr/amplicon_analysis_pipeline/blob/main/COI_and_16S_rRNA_amplicon_bioinformatic_pipeline.md) developed in the [Symbiosis Evolution Group](https://symbio.eko.uj.edu.pl/en_GB/) in December 2020. 

# Differences
The main difference between previous approach and new one is that instead of analysing **all** sequences from all the libraries at once our script analases **each library separately** and then joines them in one table.

You can run it with one click and **sip your favourite drink waiting for the results**

## Note
We are aware that our script is not perfect, it should be simplified and more user friendly. **We will get there soon!**

The most important thing for now is that:



![alt text](https://media.makeameme.org/created/its-working-oyy433.jpg)


# Usage

To run the script use:
LSD.py <sample_list.txt> <path_to_your_raw_data> <type_of_the_data>

**LSD.py** - our shiny, new script,
**sample_list** - list with libraries that you want to analyse created in following manner:
Sample_name Sample_name_R1.fastq Sample_name_R2.fastq
**path_to_your_raw_data** - path to the directory with R1 and R2 fiels for all the amplicon libraries from sample_list.txt e.g.:
/home/Data/For/Nature/Publication/)
**type_of_date** - please indicate if you are going to analyse 16SV4 or COI (only optiopns for now).
