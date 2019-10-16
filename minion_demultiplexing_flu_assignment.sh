#!/bin/bash

# New Demultiplexing with cutadapt

#First let's filter low quality reads
cat aiv_env_all.fastq | NanoFilt -q 9 > aiv_env_all_filtered.fastq

#Then convert to FASTA
sed -n '1~4s/^@/>/p;2~4p' aiv_env_all_filtered.fastq  > aiv_env_all_filtered.fasta

cutadapt -g ATCAGTTATGGTGGTCCTACCAGCAAAAGCAGG --untrimmed-output untrimmed_aiv_env_pass_forward.fasta -e 0.3 -o aiv_env_pass_forward.fasta aiv_env_all.fasta > cutadapt_ONT_forward_report.txt
#Reads with adapters: 83,654 (48.8%)

head -n 11 cutadapt_ONT_forward_report.txt >> demultiplexing_assignment_report.txt

wc -l aiv_env_pass_forward.fasta
#167308 forward_aiv_env_pass.fasta

wc -l untrimmed_aiv_env_pass_forward.fasta
#175822 untrimmed_forward_aiv_env_pass.fasta

#167308+175822 = 343130

#Matches total!
wc -l aiv_env_pass.fasta 
#343130 aiv_env_pass.fasta

cutadapt -a CCTGCTTTTGCTGGTAGGACCACCATAACTGAT --untrimmed-output untrimmed_reverse_aiv_env_pass.fasta -e 0.3 -o aiv_env_pass_reverse.fasta untrimmed_aiv_env_pass_forward.fasta > cutadapt_ONT_reverse_report.txt
#This should be -a flag (YES! AUG 22)
#Reads with adapters: 53,081 (60.4%)

head -n 11 cutadapt_ONT_reverse_report.txt >> demultiplexing_assignment_report.txt

wc -l aiv_env_pass_reverse.fasta
#106162 aiv_env_pass_reverse.fasta
#Perhaps a little fewer than expected

#Trim the Illumina Barcode from Forward reads
cutadapt -a file:illumina_rev.fasta -e 0.0 -O 21 --untrimmed-output untrimmed_barcodes_forward_aiv_env_pass.fasta -o forward_{name}.fasta aiv_env_pass_forward.fasta > cutadapt_forward_barcode_report.txt
#Reads with adapters: 11,689 (14.0%)

head -n 11 cutadapt_forward_barcode_report.txt >> demultiplexing_assignment_report.txt

wc -l untrimmed_barcodes_forward_aiv_env_pass.fasta 
#143930 untrimmed_barcodes_forward_aiv_env_pass.fasta
#Looks like a lot of reads missed!!

#Trim the Illumina Barcode from Reverse reads
cutadapt -g file:illumina_fwd.fasta -e 0.0 -O 8 --untrimmed-output untrimmed_barcodes_reverse_aiv_env_pass.fasta -o reverse_{name}.fasta aiv_env_pass_reverse.fasta > cutadapt_reverse_barcode_report.txt
#Reads with adapters: 5,280 (9.9%)

head -n 11 cutadapt_reverse_barcode_report.txt >> demultiplexing_assignment_report.txt

wc -l untrimmed_barcodes_reverse_aiv_env_pass.fasta
#95602 untrimmed_barcodes_reverse_aiv_env_pass.fasta
#Not too many...

#This was attempt to join forward and reverse reads into one file. REVISIT!! AUG 21
#for file in forward_*.fasta; 
#do
#  cat $file > reverse_${file%.fasta}.fasta
#done

#for file in forward_*.fasta; 
#do
#  echo reverse_${file%.fasta}.fasta
#done

ls forward_*.fasta | awk -F'[.]' '{print $1}'| parallel -j+0 --eta 'blastn -db ../all_avian_flu -query {.}.fasta -out {.}.out  -max_target_seqs 1 -outfmt 6'

wc -l forward_*.out >> demultiplexing_assignment_report.txt

ls reverse_*.fasta | awk -F'[.]' '{print $1}'| parallel -j+0 --eta 'blastn -db ../all_avian_flu -query {.}.fasta -out {.}.out  -max_target_seqs 1 -outfmt 6'

wc -l reverse_*.out >> demultiplexing_assignment_report.txt

#Export as text and build a script!
#Now export text files to serve as input for scripts

#Export the reads assigned to each sample, remembering to divide by 4 because it's FASTQ not FASTA
wc -l reverse_*.fasta | awk '{print $1/2" " $2}' | sed '$d'  > demultiplexing_by_sample.txt
wc -l forward_*.fasta | awk '{print $1/2" " $2}' | sed '$d' >> demultiplexing_by_sample.txt

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
  grep -f $file ~/aiv_environment/aiv_detection_environment/all_avian_flu.fasta > ${file%.out}.tab
done
