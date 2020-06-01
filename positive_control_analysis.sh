#!/bin/bash

# Script: positive_control_analysis.sh
# This file is shell script that analyzes positive control sample against the reference sequence
# NOTE: This file depends on minion_demultiplexing_flu_assignment.sh. Run that script first!

#This script is part of the following manuscript:
#Analysis of MinION sequence data from water and sediment samples for avian influenza detection in California wetlands
#This github repository includes code (and links to data) from the manuscript:
#"Linking remote sensing for targeted surveillance of Avian Influenza virus via tangential flow ultra-filtration and whole segment amplification in California wetlands"
#Madeline M. McCuen | Maurice E. Pitesky | Ana Paula da Silva | Rodrigo A. Gallardo | Jeff J. Buler | Sarai Acosta | Alexander Wilcox | Ronald F. Bond | Samuel L. Díaz-Muñoz

#Mt Siani PR8 (A/Puerto Rico/8/1934 H1N1) as a FASTA sequence to serve as a reference:
#less data/PR8_Mt_Sinai_NYC.fasta

#Let's install minimap2, which is a long-read aligner
#conda install minimap2
#minimap2-2.17

#Now let's run a command to align Oxford Nanopore reads to a reference sequence
minimap2 -ax map-ont data/PR8_Mt_Sinai_NYC.fasta forward_PR8_RNA_B.fasta > forward_PR8_RNA_B.sam
#[M::mm_idx_gen::0.004*1.59] collected minimizers
#[M::mm_idx_gen::0.006*1.97] sorted minimizers
#[M::main::0.006*1.96] loaded/built the index for 8 target sequence(s)
#[M::mm_mapopt_update::0.007*1.91] mid_occ = 2
#[M::mm_idx_stat] kmer size: 15; skip: 10; is_hpc: 0; #seq: 8
#[M::mm_idx_stat::0.007*1.87] distinct minimizers: 2496 (100.00% are singletons); average occurrences: 1.000; average spacing: 5.443
#[M::worker_pipeline::0.352*2.82] mapped 2915 sequences
#[M::main] Version: 2.17-r941
#[M::main] CMD: minimap2 -ax map-ont data/PR8_Mt_Sinai_NYC.fasta forward_PR8_RNA_B.fasta
#[M::main] Real time: 0.354 sec; CPU: 0.995 sec; Peak RSS: 0.010 GB

#Let's examine the relevant parts of the SAM alignment file
samtools stats forward_PR8_RNA_B.sam | grep ^SN | cut -f 2-
#raw total sequences:	2915
#filtered sequences:	0
#sequences:	2915
#is sorted:	0
#1st fragments:	2915
#last fragments:	0
#reads mapped:	2671
#reads mapped and paired:	0	# paired-end technology bit set + both mates mapped
#reads unmapped:	244
#reads properly paired:	0	# proper-pair bit set
#reads paired:	0	# paired-end technology bit set
#reads duplicated:	0	# PCR or optical duplicate bit set
#reads MQ0:	0	# mapped and MQ=0
#reads QC failed:	0
#non-primary alignments:	0
#total length:	2420518	# ignores clipping
#bases mapped:	2350428	# ignores clipping
#bases mapped (cigar):	2292408	# more accurate
#bases trimmed:	0
#bases duplicated:	0
#mismatches:	362012	# from NM fields
#error rate:	1.579178e-01	# mismatches / bases mapped (cigar)
#average length:	830
#maximum length:	2335
#average quality:	255.0
#insert size average:	0.0
#insert size standard deviation:	0.0
#inward oriented pairs:	0
#outward oriented pairs:	0
#pairs with other orientation:	0
#pairs on different chromosomes:	0

#So 2671/2915 = 0.9162950 mapping with 0.1579178, per base pair error rate
#NOTE: This is just from forward reads

#Now we will do the same for the reverse reads
minimap2 -ax map-ont data/PR8_Mt_Sinai_NYC.fasta reverse_PR8_RNA_B.fasta > reverse_PR8_RNA_B.sam
#[M::mm_idx_gen::0.004*1.65] collected minimizers
#[M::mm_idx_gen::0.007*2.03] sorted minimizers
#[M::main::0.007*2.02] loaded/built the index for 8 target sequence(s)
#[M::mm_mapopt_update::0.007*1.97] mid_occ = 2
#[M::mm_idx_stat] kmer size: 15; skip: 10; is_hpc: 0; #seq: 8
#[M::mm_idx_stat::0.007*1.93] distinct minimizers: 2496 (100.00% are singletons); average occurrences: 1.000; average spacing: 5.443
#[M::worker_pipeline::0.125*2.65] mapped 786 sequences
#[M::main] Version: 2.17-r941
#[M::main] CMD: minimap2 -ax map-ont data/PR8_Mt_Sinai_NYC.fasta reverse_PR8_RNA_B.fasta
#[M::main] Real time: 0.127 sec; CPU: 0.332 sec; Peak RSS: 0.006 GB

#Let's examine the relevant parts of the SAM alignment file
samtools stats reverse_PR8_RNA_B.sam | grep ^SN | cut -f 2-
#raw total sequences:	786
#filtered sequences:	0
#sequences:	786
#is sorted:	0
#1st fragments:	786
#last fragments:	0
#reads mapped:	739
#reads mapped and paired:	0	# paired-end technology bit set + both mates mapped
#reads unmapped:	47
#reads properly paired:	0	# proper-pair bit set
#reads paired:	0	# paired-end technology bit set
#reads duplicated:	0	# PCR or optical duplicate bit set
#reads MQ0:	0	# mapped and MQ=0
#reads QC failed:	0
#non-primary alignments:	0
#total length:	665575	# ignores clipping
#bases mapped:	651427	# ignores clipping
#bases mapped (cigar):	635008	# more accurate
#bases trimmed:	0
#bases duplicated:	0
#mismatches:	99861	# from NM fields
#error rate:	1.572594e-01	# mismatches / bases mapped (cigar)
#average length:	846
#maximum length:	2040
#average quality:	255.0
#insert size average:	0.0
#insert size standard deviation:	0.0
#inward oriented pairs:	0
#outward oriented pairs:	0
#pairs with other orientation:	0
#pairs on different chromosomes:	0

#So 739/786 = 0.9402036 mapping with 0.1572594, per base pair error rate
#NOTE: This is just from reverse reads

#Now all PR8_RNA_B reads together
cat forward_PR8_RNA_B.fasta reverse_PR8_RNA_B.fasta > PR8_RNA_B.fasta
minimap2 -ax map-ont data/PR8_Mt_Sinai_NYC.fasta PR8_RNA_B.fasta > PR8_RNA_B.sam
samtools stats PR8_RNA_B.sam | grep ^SN | cut -f 2-
#In total 92.137% mapping with 0.158 per base error rate

#Now we will do a quick sanity check with an actual sample, which we don't expect to map very well
cat forward_B1_SW_Unfltrd.fasta reverse_B1_SW_Unfltrd.fasta > B1_SW_Unfltrd.fasta
minimap2 -ax map-ont data/PR8_Mt_Sinai_NYC.fasta forward_B1_SW_Unfltrd.fasta > B1_SW_Unfltrd.sam
samtools stats B1_SW_Unfltrd.sam | grep ^SN | cut -f 2-
#So 17.452% mapping with 0.1405933 error rate (which is irrelevant because this is not correct reference)
#This lower mapping percentage reflects the fact we are mapping environmental influenza samples (that match the avian influenza virus genome database) against a human influenza virus genome

#In summary, we have very high mapping percentage (>92%) of our control sample to the reference genome from which they were derived
#Mapping an experimental sample to the control reference genome predictably leads to a drop in the percent reads aligning (i.e. we don't expect human influenza virus necessarily),
  #but we do get some alignment as we do have influenza virus in our sample 
#This analysis provides another line of evidence (in addition to negative controls) to boost confidence that demultiplexing was conducted properly.