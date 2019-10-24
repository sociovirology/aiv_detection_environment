# Analysis of MinION sequence data from water and sediment samples for avian influenza detection in California wetlands

# Code for "Linking remote sensing for targeted surveillance of Avian Influenza virus via tangential flow ultra-filtration and whole segment amplification in California wetlands"

This github repository includes code (and links to data) from the manuscript: 
Code for "Linking remote sensing for targeted surveillance of Avian Influenza virus via tangential flow ultra-filtration and whole segment amplification in California wetlands"
Madeline M. McCuen | Maurice E. Pitesky | Ana Paula da Silva | Rodrigo A. Gallardo | Jeff J. Buler | Sarai Acosta | Alex Wilcox | Ronald F. Bond | Samuel L. Díaz-Muñoz 

If you are reading or using this, let us know how these data were useful for you. Always open to collaborate! Please contact us!

If you use these data and code, please cite the repository or the paper.

##Quick Start
1. Make sure packages are installed (see #2 below)
2. git clone https://github.com/sociovirology/aiv_detection_environment.git
3. chmod +x minion_demultiplexing_flu_assignment.sh
4. ./minion_demultiplexing_flu_assignment.sh
5. Rscript aiv_detection_environment_analysis.R (or load interactively in R)
6. chmod +x env_samples_qPCR_primer_analysis.sh ; ./env_samples_qPCR_primer_analysis.sh (depends on 4)
7. chmod +x AIV_qPCR_primer_analysis.sh ; ./AIV_qPCR_primer_analysis.sh (can be run independently of above)

### CONTENTS
### 1. Project Description
### 2. Packages and software used to test code
### 3. Data
### 4. Code

### 1. Project Description
Birds are the "natural" host for influenza A viruses (especially duck, geese, and other waterfowl). These avian influenza viruses can spread to commercial poultry operations and have been involved in creating global pandemics after mixing with human influenza strains. Influenza is a digestive infection in birds, which leads to sheding considerable amounts of influenza virus into bodies of water. Therefore, having a good idea of this environmental reservior of viruses is important for fod safety and human health. We used water filtration and single molecule long-read sequencing (Oxford Nanopore Technologies MinION sequencer) to improve detection of avian influenza viruses over current, standard quantitative PCR method on unfiltered water. This repository includes all data (or links to) and analyses in the paper.

Abstract:
Migratory waterfowl, including geese and ducks, are indicated as the primary reservoir of avian influenza viruses (AIv) which can be subsequently spread to commercial poultry. The US Department of Agriculture’s (USDA) surveillance efforts of waterfowl for AIv have been largely discontinued in the contiguous United States. Consequently, the use of technologies to identify areas of high waterfowl density and detect the presence of AIv in habitat such as wetlands has become imperative. Here we identified two high waterfowl density areas in California using processed NEXt generation RADar (NEXRAD) and collected water samples to test the efficacy of two tangential flow ultrafiltration methods and two nucleic acid based AIv detection assays. Whole-segment M-RTPCR yielded more positive samples than standard M-segment qPCR methods (50.0% vs 2.6%, p<0.0001). We determined that this difference in positivity was due to mismatches in published primers to our samples and that these mismatches would result in failing to detect in the vast majority of currently sequenced AIv genomes in public databases. The whole segment sequences were subsequently used to provide subtype and potential host information of the AIv environmental reservoir. There was no statistically significant difference in sequencing reads recovered from the RexeedTM filtration compared to the unfiltered surface water. This overall approach combining remote sensing, filtration, and sequencing provides a novel and potentially more effective, surveillance approach for AIv.

### 2. Data
Data consists of sequencing output from the Oxford Nanopore Technologies MinION sequencer platform (FASTQ files), sample information incl ph/temperature/salinity measurements, gels, and Influenza Research Database FASTA files, and sample barcodes

1) Sequencing file is available in the Sequencing Read Archive

2) Sample information is in data/AIV_filtration_sample_key_9_26_18.csv 

3) Data file with temperature, pH, and salinity data is in data/Temp_pH_Salinitt_for_B&C.csv

4) Primer and barcode information is in data/illumina_fwd.fasta and data/illumina_rev.fasta

5) Database file of all avian influenza sequences in FluDB in data/all_avian_flu.fasta.gz (used in minion_demultiplexing_flu_assignment.sh)

6) Database file of all avian influenza sequences in FluDB updated on August 12, 2019 data/avian_complete_genomes_aug12_2019.fasta.gz (used in AIV_qPCR_primer_analysis.sh)

### 3. Code
Below are descriptions of the code files used to generate the tables, figures, and statistics in the paper.

1) minion_demultiplexing_flu_assignment.sh: This file is shell script that downloads and processes raw sequencing reads 

2) aiv_detection_environment_analysis.R: This file is an R script that couples sequencing information with sample information to arive at the main conclusions in the manuscript

3) env_samples_qPCR_primer_analysis.sh: This file is a shell script that evaluates whether the sequences derived from samples will yield positives using the Spackman et al 2003 M-segment qPCR primers 

4) AIV_qPCR_primer_analysis.sh: This shell script uses published sequences in the Influenza Research Database (http://fludb.org) to verify if Spackman et al 2003 primers (for M-segment qPCR) are expected to work on M segments of avian influenza virus