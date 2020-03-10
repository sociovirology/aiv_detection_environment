#!/bin/bash

# Script: AIV_qPCR_primer_analysis.sh
# This shell script uses published sequences in the Influenza Research Database (http://fludb.org) to verify if Spackman et al 2003 primers (for M-segment qPCR) are expected to work on M segments of avian influenza virus

#This script is part of the following manuscript:
#Analysis of MinION sequence data from water and sediment samples for avian influenza detection in California wetlands
#This github repository includes code (and links to data) from the manuscript:
#"Linking remote sensing for targeted surveillance of Avian Influenza virus via tangential flow ultra-filtration and whole segment amplification in California wetlands"
#Madeline M. McCuen | Maurice E. Pitesky | Ana Paula da Silva | Rodrigo A. Gallardo | Jeff J. Buler | Sarai Acosta | Alexander Wilcox | Ronald F. Bond | Samuel L. Díaz-Muñoz

################ Avian Influenza Virus M-segment qPCR Evaluation ################
#### 1. Script Description, Data Sources, and Preparing Data
#### 2. Method A: Exact sequence match analysis using GREP
#### 3. Method B: Biochemical simulation of qPCR using primers against sequence database using ThermonucleotideBLAST 
#### 4. Conclusions 

#### 1. Script Description, Data Sources, and Preparing Data ####

# This shell script uses published sequences in the Influenza Research Database (http://fludb.org) to verify if Spackman et al 2002 (see below) primers are expected to work on M segments of avian influenza virus 
#Spackman, Erica, et al. "Development of a real-time reverse transcriptase PCR assay for type A influenza virus and the avian H5 and H7 hemagglutinin subtypes." Journal of clinical microbiology 40.9 (2002): 3256-3260.
# The M primers are (from Table 1) :
# Influenza A virus	
#  M + 25	      AGA TGA GTC TTC TAA CCG AGG TCG
#	 M − 124	    TGC AAA AAC ATC TTC AAG TCT CTG
#	 M + 64	  FAM-TCA GGC CCC CTC AAA GCC GA-TAMRA 

# Without spaces: 
#   M+25:   AGATGAGTCTTCTAACCGAGGTCG
#   M−124:	TGCAAAAACATCTTCAAGTCTCTG REVCOMPL> CAGAGACTTGAAGATGTTTTTGCA
#   M+64:	  FAM-TCAGGCCCCCTCAAAGCCGA-TAMRA

#The primary data is a FASTA file of all avian influenza virus complete genomes available in the Influenza Research Database (http://fludb.org) - downloaded Aug 12, 2019
#Fasta available in data subfolder. Let's unzip that file
gunzip data/avian_complete_genomes_aug12_2019.fasta.gz

#FASTA file is multi-line: fix and export to working folder

wc -l data/avian_complete_genomes_aug12_2019.fasta
#2520795 data/avian_complete_genomes_aug12_2019.fasta

grep ">" -c data/avian_complete_genomes_aug12_2019.fasta
#93766
#Ok something clearly wrong, way more lines than read entries according to > character for FASTA

#Remove newlines script from https://stackoverflow.com/questions/15857088/remove-line-breaks-in-a-fasta-file
awk '!/^>/ { printf "%s", $0; n = "\n" } 
/^>/ { print n $0; n = "" }
END { printf "%s", n }
' data/avian_complete_genomes_aug12_2019.fasta > fixed_all_avian_flu.fasta

#Okay let's try again!
wc -l fixed_all_avian_flu.fasta
#187532 fixed_all_avian_flu.fasta

wc -l fixed_all_avian_flu.fasta | awk '{print $1/2}'
#93766

grep ">" -c fixed_all_avian_flu.fasta
#93766
#Okay, now we're cooking


#### 2. Method A: Exact sequence match analysis using GREP ####
#First we'll check for exact matches to the Spackman primers. Recall we have 93766 segment sequences total
#Out of those, we have the following amount of M segment complete sequences:
grep "Segment:7" -c fixed_all_avian_flu.fasta
#9855
#This will be the relevant denominator for the analyses below 

#Forward Primer Exact Matches
grep "AGATGAGTCTTCTAACCGAGGTCG" -c fixed_all_avian_flu.fasta
#6348

#Reverse Primer Exact Matches
grep "CAGAGACTTGAAGATGTTTTTGCA" -c fixed_all_avian_flu.fasta
#538

#Probe Primer Exact Matches
grep "TCAGGCCCCCTCAAAGCCGA" -c fixed_all_avian_flu.fasta
#9630

#Obviously we need these primers/probe to be on the same molecule, so adjusting the grep search
grep AGATGAGTCTTCTAACCGAGGTCG.*TCAGGCCCCCTCAAAGCCGA.*CAGAGACTTGAAGATGTTTTTGCA -c fixed_all_avian_flu.fasta
#232

#Check if these are the expected size amplicon size of ~100bp
grep AGATGAGTCTTCTAACCGAGGTCG.*TCAGGCCCCCTCAAAGCCGA.*CAGAGACTTGAAGATGTTTTTGCA -o fixed_all_avian_flu.fasta | awk '{ print length }'
#101 (for every line), so as expected

#This analysis would suggest a overall success rate of 232/9855 = 0.02354135 ~ 2.35%, based on exact primer matches  

#We might want to make the search more flexible by excluding the 5' ends of the primers that are less important for annealing 

#Forward Primer Only first 16nt from 3'
grep "AGATGAGTCTTCTAAC" -c fixed_all_avian_flu.fasta
#6572

#Reverse Primer Only first 16nt from 3'
grep "CAGAGACTTGAAGATG" -c fixed_all_avian_flu.fasta
#6876

#Probe Only first 16nt from 3'
grep "TCAGGCCCCCTCAAAG" -c fixed_all_avian_flu.fasta
#9634

#Now all together, adjusting the grep search for primers/probe to be on the same molecule, so 
grep AGATGAGTCTTCTAAC.*TCAGGCCCCCTCAAAG.*CAGAGACTTGAAGATG -c fixed_all_avian_flu.fasta
#5067


#Quick check if we're getting expected amplicon sizes
grep AGATGAGTCTTCTAAC.*TCAGGCCCCCTCAAAG.*CAGAGACTTGAAGATG -o fixed_all_avian_flu.fasta | awk '{ print length }'
#93 for every line, so as expected

#This analysis would suggest a overall success rate of 5067/9855 = 0.5141553 ~ 51.16%, based on exact matches to the first 16nt of the primers/probe

#So based on this analysis at best half of known avian influenza viruses would be detectable 
#However, sucess rate could be higher given that PCR is able to tolerate mismatches. So examining a second method 


#### 3. Method B: Biochemical simulation of qPCR using primers against sequence database using ThermonucleotideBLAST ####
#This method uses ThermonucleotideBLAST (https://public.lanl.gov/jgans/tntblast/tntblast_doc.html) 
#Citation for ThermonucleotideBLAST: J. D. Gans and M. Wolinsky "Improved assay-dependent searching of nucleic acid sequence databases" Nucleic Acid Res. 2008 Jul;36(12):e74
#ThermonucleotideBLAST is a software program for searching a target database of nucleic acid sequences using an assay specific query. ThermonucleotideBLAST queries are based on biochemical assays (i.e. a pair of oligonucleotide sequences representing PCR primers or Padlock probes, a triplet of oligos representing PCR primers and a TaqManTM probe or a single oligo representing a hybridization probe). Unlike existing programs (e.g. BLAST) which use heuristic measures of sequence similarity for identifying matches between a query and target sequence, ThermonucleotideBLAST uses physically relevant measures of sequence similarity -- free energy and melting temperature. For example, given a pair of PCR primers, a database of DNA targets and an annealing temperature, ThermonucleotideBLAST will return a list of predicted amplicons that will (ideally) match experimental PCR results.

#Code below assumes you have ThermonucleotideBLAST installed as a module. 
module load tntblast
#Module tntblast/2.04 loaded

#First step is to create files with the primer sequences and probes
echo "TaqmanqPCR-msegment    AGATGAGTCTTCTAACCGAGGTCG    TGCAAAAACATCTTCAAGTCTCTG    TCAGGCCCCCTCAAAGCCGA" > taqman_qPCR_query.txt

#Create second file with just the amplicon information (regular qPCR or PCR)
echo "qPCR-msegment AGATGAGTCTTCTAACCGAGGTCG    TGCAAAAACATCTTCAAGTCTCTG" > qPCR_query.txt

#First evaluate the Taqman assay, using Spackman et al 2002 conditions of 60 degrees C annealing temp; 10 pmol of each primer and 0.3 μM probe in a 20uL reaction 
tntblast -i taqman_qPCR_query.txt -d fixed_all_avian_flu.fasta -e 60 -E 60 -m 1 -t 0.0000003 -T 0.0000005 -o taqman_qPCR_query_60.fasta
#No matches

#Check with more permissive temperature 50 degrees C
tntblast -i taqman_qPCR_query.txt -d fixed_all_avian_flu.fasta -e 50 -E 50 -m 1 -t 0.0000003 -T 0.0000005 -o taqman_qPCR_query_50.fasta
#Found 5144 (total) matches
#Amplicon:
#	101 <= Amplicon length <= 101

#Even more permissive conditions, 40 degrees C
tntblast -i taqman_qPCR_query.txt -d fixed_all_avian_flu.fasta -e 40 -E 40 -m 1 -t 0.0000003 -T 0.0000005 -o taqman_qPCR_query_40.fasta
#Found 9050 (total) matches
#Amplicon:
#	101 <= Amplicon length <= 101

#This analysis would suggest no amplification based on published PCR conditions.
#Under more permissive conditions (50 degree C annealing), the analysis suggests an overall success rate of 5144/9855 = 0.5219685 ~ 52.20%
#At the an annealing temperature to 40 degrees C, the analysis suggests an overall success rate of 9050/9855 = 0.9183156 ~ 91.83%. Despite low annealing temperature predicted amplicons were still in the expected size

#So based on this analysis, M-Segment Taqman qPCR (as published by Spackman et al 2002) will detect very few AIV M segments. Reducing annealing temperatures could increase detection success substantially


#Second let's evaluate the amplification of the amplicon (no probe), but otherwise using Spackman et al 2002 conditions of 60 degrees C annealing temp; 10 pmol of each primer and 0.3 μM probe in a 20uL reaction 
tntblast -i qPCR_query.txt -d fixed_all_avian_flu.fasta -e 60 -m 1 -t 0.0000003 -o qPCR_query_60.fasta
#Found 0 (total) matches

tntblast -i qPCR_query.txt -d fixed_all_avian_flu.fasta -e 50 -m 1 -t 0.0000003 -o qPCR_query_50.fasta
#Found 5145 (total) matches
#Amplicon:
#	101 <= Amplicon length <= 101

tntblast -i qPCR_query.txt -d fixed_all_avian_flu.fasta -e 40 -m 1 -t 0.0000003 -o qPCR_query_40.fasta
#Found 9051 (total) matches
#Amplicon:
#	101 <= Amplicon length <= 101

#So based on this analysis, the qPCR without the probe, would have identical success to the Taqman assay

#### 4. Conclusions ####
#Regarding the Spackman et al. 2003 M-segment Taqman assay:
#About 2% of sequenced AIV strains have exact matches to the primers and probe
#About 50% of sequenced AIV strains have exact matches to the first 16nt of the primers and probe
#Very few if any published sequenced strains are expected to yield positives with published conditions (60ºC annealing temp) based on assay simulation 
#Under less stringent annealing conditions, the assay simulation is expected to generate more positives (~50% for 50C and ~90% for 40C)
