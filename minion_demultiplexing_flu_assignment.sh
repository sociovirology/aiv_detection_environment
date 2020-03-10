#!/bin/bash

# Script: minion_demultiplexing_flu_assignment.sh
# This file is shell script that downloads and processes raw sequencing reads

#This script is part of the following manuscript:
#Analysis of MinION sequence data from water and sediment samples for avian influenza detection in California wetlands
#This github repository includes code (and links to data) from the manuscript:
#"Linking remote sensing for targeted surveillance of Avian Influenza virus via tangential flow ultra-filtration and whole segment amplification in California wetlands"
#Madeline M. McCuen | Maurice E. Pitesky | Ana Paula da Silva | Rodrigo A. Gallardo | Jeff J. Buler | Sarai Acosta | Alexander Wilcox | Ronald F. Bond | Samuel L. Díaz-Muñoz

#Download FASTQ sequencing file from SRA via EBI (Guppy HAC basecalling)
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/013/SRR10302213/SRR10302213_1.fastq.gz -O data/aiv_env_all.fastq.gz

#Unzip 
gunzip data/aiv_env_all.fastq.gz

#Also unzip the FluDB FASTA which we will use to create the BLAST database downstream
gunzip data/all_avian_flu.fasta.gz

#First let's filter low quality reads
cat data/aiv_env_all.fastq | NanoFilt -q 9 > aiv_env_all_filtered.fastq

#Then convert to FASTA
#next line doesn't work on OS X
#sed -n '1~4s/^@/>/p;2~4p' aiv_env_all_filtered.fastq  > aiv_env_all_filtered.fasta
seqtk seq -A aiv_env_all_filtered.fastq > aiv_env_all_filtered.fasta

cutadapt -g ATCAGTTATGGTGGTCCTACCAGCAAAAGCAGG --untrimmed-output untrimmed_aiv_env_pass_forward.fasta -e 0.3 -o aiv_env_pass_forward.fasta aiv_env_all_filtered.fasta > cutadapt_ONT_forward_report.txt
#Total reads processed:                 143,716
#Reads with adapters:                    89,571 (62.3%)

#Save barcode demultiplexing info into report
head -n 11 cutadapt_ONT_forward_report.txt >> demultiplexing_assignment_report.txt

wc -l aiv_env_pass_forward.fasta
#179142 aiv_env_pass_forward.fasta

wc -l untrimmed_aiv_env_pass_forward.fasta
#179142untrimmed_aiv_env_pass_forward.fasta

#179142+108290 = 287432
#Matches total!
wc -l aiv_env_all_filtered.fasta 
#287432 aiv_env_all_filtered.fasta

cutadapt -a CCTGCTTTTGCTGGTAGGACCACCATAACTGAT --untrimmed-output untrimmed_reverse_aiv_env_pass.fasta -e 0.3 -o aiv_env_pass_reverse.fasta untrimmed_aiv_env_pass_forward.fasta > cutadapt_ONT_reverse_report.txt
#Total reads processed:                  54,145
#Reads with adapters:                    39,737 (73.4%)

#Save barcode demultiplexing info into report
head -n 11 cutadapt_ONT_reverse_report.txt >> demultiplexing_assignment_report.txt

wc -l aiv_env_pass_reverse.fasta
#79474 aiv_env_pass_reverse.fasta

#Trim the Illumina Barcode from Forward reads
cutadapt -a file:data/illumina_rev.fasta -e 0.05 -O 21 --untrimmed-output untrimmed_barcodes_forward_aiv_env_pass.fasta -o forward_{name}.fasta aiv_env_pass_forward.fasta > cutadapt_forward_barcode_report.txt
#Total reads processed:                  89,571
#Reads with adapters:                    38,331 (42.8%)

head -n 11 cutadapt_forward_barcode_report.txt >> demultiplexing_assignment_report.txt

wc -l untrimmed_barcodes_forward_aiv_env_pass.fasta 
#102480 untrimmed_barcodes_forward_aiv_env_pass.fasta

#Trim the Illumina Barcode from Reverse reads
cutadapt -g file:data/illumina_fwd.fasta -e 0.0 -O 8 --untrimmed-output untrimmed_barcodes_reverse_aiv_env_pass.fasta -o reverse_{name}.fasta aiv_env_pass_reverse.fasta > cutadapt_reverse_barcode_report.txt
#Total reads processed:                  39,737
#Reads with adapters:                     9,631 (24.2%)

head -n 11 cutadapt_reverse_barcode_report.txt >> demultiplexing_assignment_report.txt

wc -l untrimmed_barcodes_reverse_aiv_env_pass.fasta
#60212 untrimmed_barcodes_reverse_aiv_env_pass.fasta

#Demultiplexing complete. Now export text files to serve as input for downstream scripts
#Export the reads assigned to each sample, remembering to divide by 2 because it's FASTA
wc -l reverse_*.fasta | awk '{print $1/2" " $2}' | sed '$d'  > demultiplexing_by_sample.txt
wc -l forward_*.fasta | awk '{print $1/2" " $2}' | sed '$d' >> demultiplexing_by_sample.txt

#Now looking up sequences against all avian influenza virus whole genome sequences from FluDB

#Need to make a BLAST database using the FluDB
#File is at data/all_avian_flu.fasta
makeblastdb -in data/all_avian_flu.fasta -dbtype nucl -out all_avian_flu

#This command uses GNU parallel to conduct BLAST searches in parallel for each demultiplexed file
#First looking at forward read files
ls forward_*.fasta | awk -F'[.]' '{print $1}'| parallel -j+0 --eta 'blastn -db all_avian_flu -query {.}.fasta -out {.}.out  -max_target_seqs 1 -outfmt 6'

#Export info to our running report
wc -l forward_*.out >> demultiplexing_assignment_report.txt

##Second, let's do the same for reverse read files
ls reverse_*.fasta | awk -F'[.]' '{print $1}'| parallel -j+0 --eta 'blastn -db all_avian_flu -query {.}.fasta -out {.}.out  -max_target_seqs 1 -outfmt 6'

#Export info to our running report
wc -l reverse_*.out >> demultiplexing_assignment_report.txt

#Now export text files to serve as input for downstream scripts

#Now add the number of reads matching the avian flu database to the file avian_blast_matches.txt
wc -l reverse_*.out | awk '{print $1" " $2}' | sed '$d' > avian_blast_matches.txt
wc -l forward_*.out | awk '{print $1" " $2}' | sed '$d' >> avian_blast_matches.txt

#Finally look up metadata for the closest blast matches for each read
#Now use the .out files to generate a list of GI's to lookup against database
for file in *.out; 
do
  awk '{print $2}' $file | awk -F ":" '{print $2}' | awk -F "|" '{print $1}' > ${file%.out}.tmp
done

#Use that file to get matches from the database
for file in *.tmp; 
do
  grep -f $file data/all_avian_flu.fasta > "${file%.out}.tab"
done
