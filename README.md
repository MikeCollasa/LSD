# LSD
LSD: Library Single Denoising


Here we provide the details of the updated amplicon data analysis workflow developed in February 2021 by Michał Kolasa, Catalina Valdivia and Sandra Ahlén Mulio. Our pipeline is [strongly based on one developed](https://github.com/symPiotr/amplicon_analysis_pipeline/blob/main/COI_and_16S_rRNA_amplicon_bioinformatic_pipeline.md) in the [Symbiosis Evolution Group](https://symbio.eko.uj.edu.pl/en_GB/) in December 2020. 

# Differences
The main difference between previous approach and new one is that instead of analysing **all** sequences from all the libraries at once our script analases **each library separately** and then joines them in one table.
Also our apprach is Python-based, so you can just customize our script with info about your directories, run it and **sip your favourite drink waiting for the results**

## Note
We are aware that our script is not perfect, it should be simplified and more user friendly. **We will get there soon**
