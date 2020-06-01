# Analysis of MinION sequence data from water and sediment samples for avian influenza detection in California wetlands

This github repository includes code (and links to data) from the manuscript:  
"A comparison of amplification methods to detect Avian Influenza viruses in California wetlands targeted via remote sensing of waterfowl"  
Madeline M. McCuen | Maurice E. Pitesky | Jeff J. Buler | Sarai Acosta | Alexander Wilcox | Ronald F. Bond | Samuel L. Díaz-Muñoz

If you are reading or using this, let us know how these data were useful for you. If you use these data and code, please cite the repository or the paper. Always open to collaborate! Please contact us!

### Quick Start
1. Make sure packages are installed (see #2 below)
2. git clone https://github.com/sociovirology/aiv_detection_environment.git
3. chmod +x minion_demultiplexing_flu_assignment.sh
4. ./minion_demultiplexing_flu_assignment.sh
5. Rscript aiv_detection_environment_analysis.R (or load interactively in R)
6. chmod +x env_samples_qPCR_primer_analysis.sh ; ./env_samples_qPCR_primer_analysis.sh (depends on 4)
7. chmod +x AIV_qPCR_primer_analysis.sh ; ./AIV_qPCR_primer_analysis.sh (can be run independently of above)
8. chmod +x positive_control_analysis.sh; ./positive_control_analysis.sh (depends on 4)

### CONTENTS
1. Project Description
2. Packages and software used to test code
3. Data
4. Code

### 1. Project Description
Birds are the "natural" host for influenza A viruses (especially duck, geese, and other waterfowl). These avian influenza viruses can spread to commercial poultry operations and have been involved in creating global pandemics after mixing with human influenza strains. Influenza is a digestive infection in birds, which leads to shedding considerable amounts of influenza virus into bodies of water. Therefore, having a good idea of this environmental reservior of viruses is important for food safety and human health. We used water filtration and single molecule long-read sequencing (Oxford Nanopore Technologies MinION sequencer) to improve detection of avian influenza viruses over the current, standard quantitative PCR method on unfiltered water. This repository includes all data (or links to) and analyses in the paper.

Abstract:
Migratory waterfowl, including geese and ducks, are indicated as the primary res- ervoir of avian influenza viruses (AIv) which can be subsequently spread to com- mercial poultry. The US Department of Agriculture's (USDA) surveillance efforts of waterfowl for AIv have been largely discontinued in the contiguous United States. Consequently, the use of technologies to identify areas of high waterfowl density and detect the presence of AIv in habitat such as wetlands has become imperative. Here we identified two high waterfowl density areas in California using processed NEXt generation RADar (NEXRAD) and collected water samples to test the efficacy of two tangential flow ultrafiltration methods and two nucleic acid based AIv detection as- says. Whole-segment amplification and long-read sequencing yielded more positive samples than standard M-segment qPCR methods (57.6% versus 3.0%, p < .0001). We determined that this difference in positivity was due to mismatches in published primers to our samples and that these mismatches would result in failing to detect in the vast majority of currently sequenced AIv genomes in public databases. The whole segment sequences were subsequently used to provide subtype and potential host information of the AIv environmental reservoir. There was no statistically significant difference in sequencing reads recovered from the RexeedTM filtration compared to the unfiltered surface water. This overall approach combining remote sensing, filtra- tion and sequencing provides a novel and potentially more effective, surveillance approach for AIv.

### 2. Packages and software used to test code
This code was tested using the following software packages:

* NanoFilt (2.3.0)
* cutadapt (2.4)
* BLAST (2.2.29)
* R (3.5.0 (2018-04-23) -- "Joy in Playing") with packages:
    + readr, tidyr, ggplot2, dplyr, ggthemes, gridExtra, reshape2
* tntblast (2.04)
* seqtk (1.3)
* minimap2 (2.17)
* GNU Parallel 2018 (needs to be installed on OS X: brew install parallel)
* usearch (8.1.1861, only for checking_similar_sequences.sh)

### 3. Data
Data consists of sequencing output from the Oxford Nanopore Technologies MinION sequencer platform (FASTQ files), sample information incl ph/temperature/salinity measurements, gels, and Influenza Research Database FASTA files, and sample barcodes

1) Sequencing file is available in the Sequence Read Archive (Accession SRX7014890)

2) Sample information is in data/AIV_filtration_sample_key_9_26_18.csv 

3) Data file with temperature, pH, and salinity data is in data/Temp_pH_Salinitt_for_B&C.csv

4) Primer and barcode information is in data/illumina_fwd.fasta and data/illumina_rev.fasta

5) Database file of all avian influenza sequences in FluDB in data/all_avian_flu.fasta.gz (used in minion_demultiplexing_flu_assignment.sh)

6) Database file of all avian influenza sequences in FluDB updated on August 12, 2019 data/avian_complete_genomes_aug12_2019.fasta.gz (used in AIV_qPCR_primer_analysis.sh)

7) Reference sequence for the positive control sample (A/Puerto Rico/8/1934 H1N1) is in data/PR8_Mt_Sinai_NYC.fasta

### 4. Code
Below are descriptions of the code files used to generate the tables, figures, and statistics in the paper.

1) minion_demultiplexing_flu_assignment.sh: This file is shell script that downloads and processes raw sequencing reads 

2) aiv_detection_environment_analysis.R: This file is an R script that couples sequencing information with sample information to arive at the main conclusions in the manuscript

3) env_samples_qPCR_primer_analysis.sh: This file is a shell script that evaluates whether the sequences derived from samples will yield positives using the Spackman et al 2003 M-segment qPCR primers 

4) AIV_qPCR_primer_analysis.sh: This shell script uses published sequences in the Influenza Research Database (http://fludb.org) to verify if Spackman et al 2003 primers (for M-segment qPCR) are expected to work on M segments of avian influenza virus

5) positive_control_analysis.sh: This shell script analyzes positive control sample against the reference sequence. Note that this file depends on minion_demultiplexing_flu_assignment.sh

6) checking_similar_sequences.sh: This shell script checks if there are identical/very similar reads in the sequencing dataset