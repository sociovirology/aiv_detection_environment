#!/bin/bash

# Script: checking_similar_sequences.sh
# This file checks if there are identical/very similar reads in the dataset
# NOTE: This file depends on minion_demultiplexing_flu_assignment.sh. Run that script first!

#This script is part of the following manuscript:
#Analysis of MinION sequence data from water and sediment samples for avian influenza detection in California wetlands
#This github repository includes code (and links to data) from the manuscript:
#"Linking remote sensing for targeted surveillance of Avian Influenza virus via tangential flow ultra-filtration and whole segment amplification in California wetlands"
#Madeline M. McCuen | Maurice E. Pitesky | Ana Paula da Silva | Rodrigo A. Gallardo | Jeff J. Buler | Sarai Acosta | Alexander Wilcox | Ronald F. Bond | Samuel L. Díaz-Muñoz

#Start with quality filtered, demultiplexed, and trimmed reads
cat forward_*.fasta > forward_pooled_seqs.fasta
cat reverse_*.fasta > reverse_pooled_seqs.fasta
cat forward_pooled_seqs.fasta reverse_pooled_seqs.fasta > pooled_seqs.fasta

#Quick check that it worked
wc -l *pooled_seqs.fasta
#   76662 forward_pooled_seqs.fasta
#   95924 pooled_seqs.fasta
#   19262 reverse_pooled_seqs.fasta
#  191848 total
#It did!!!

#Load usearch module
module load usearch/8.1.1861 

#Sort 
usearch -sortbylength pooled_seqs.fasta -fastaout pooled_seqs_sorted.fasta -minseqlength 400
#00:01  72Mb  100.0% Reading pooled_seqs.fasta
#00:01  39Mb Sorting by length                
#00:01  39Mb Length min 1, median 484, max 9027
#00:01  39Mb  100.0% Writing output 12375 short, 0 long
#00:01  39Mb Write done, closing file and exiting

usearch -cluster_fast pooled_seqs_sorted.fasta -id 0.98 -centroids nr.fasta -uc clusters.uc
#35578 seqs, 35578 uniques, 35578 singletons (100.0%)
#00:01  78Mb Min size 1, median 1, max 1, avg 1.00
#00:01  74Mb done.
#00:36 258Mb  100.0% 35578 clusters, max size 1, avg 1.0
#00:36 259Mb  100.0% Writing centroids to nr.fasta      
                                                 
#      Seqs  35578 (35578 ()
#  Clusters  35578 (35578 (35578 (3)
#  Max size  1
#  Avg size  1.0
#  Min size  1
#Singletons  35578 (35578 (35578 (35578 (35), 100.0% of seqs, 100.0% of clusters
#   Max mem  275Mb
#      Time  36.0s
#Throughput  988.3 seqs/sec.

#First usearch tells us there are no identical sequences (35578 seqs, 35578 uniques, 35578 singletons)
#Furthermore after clustering at 98% identity all clusters are singletons, so no sequences are 98% similar

usearch -cluster_fast pooled_seqs_sorted.fasta -id 0.97 -centroids nr.fasta -uc clusters.uc
#35578 seqs, 35578 uniques, 35578 singletons (100.0%)
#00:00  78Mb Min size 1, median 1, max 1, avg 1.00
#00:00  74Mb done.
#00:41 253Mb  100.0% 35575 clusters, max size 3, avg 1.0
#00:41 254Mb  100.0% Writing centroids to nr.fasta      
                                                 
#      Seqs  35578 (35578 ()
#  Clusters  35575 (35575 (35575 (3)
#  Max size  3
#  Avg size  1.0
#  Min size  1
#Singletons  35573 (35573 (35573 (35573 (35), 100.0% of seqs, 100.0% of clusters
#   Max mem  262Mb
#      Time  41.0s
#Throughput  867.8 seqs/sec

#Clustering at 97% identity nearly all clusters are singletons, the max has 3 sequences

usearch -cluster_fast pooled_seqs_sorted.fasta -id 0.96 -centroids nr.fasta -uc clusters.uc
#35578 seqs, 35578 uniques, 35578 singletons (100.0%)
#00:01  78Mb Min size 1, median 1, max 1, avg 1.00
#00:01  74Mb done.
#00:47 253Mb  100.0% 35528 clusters, max size 7, avg 1.0
#00:48 253Mb  100.0% Writing centroids to nr.fasta      
                                                 
#      Seqs  35578 (35578 ()
#  Clusters  35528 (35528 (35528 (3)
#  Max size  7
#  Avg size  1.0
#  Min size  1
#Singletons  35497 (35497 (35497 (35497 (35), 99.8% of seqs, 99.9% of clusters
#   Max mem  274Mb
#      Time  47.0s
#Throughput  757.0 seqs/sec.

#Clustering at 96% identity nearly all clusters are singletons, the max has 7 sequences

#Okay let's cluster at 95% ID 
usearch -cluster_fast pooled_seqs_sorted.fasta -id 0.95 -centroids nr.fasta -uc clusters.uc
#35578 seqs, 35578 uniques, 35578 singletons (100.0%)
#00:01  78Mb Min size 1, median 1, max 1, avg 1.00
#00:01  74Mb done.
#00:52 259Mb  100.0% 35351 clusters, max size 20, avg 1.0
#00:53 259Mb  100.0% Writing centroids to nr.fasta
                                          
#      Seqs  35578 (35578 ()
#  Clusters  35351 (35351 (35351 (3)
#  Max size  20
#  Avg size  1.0
#  Min size  1
#Singletons  35227 (35227 (35227 (35227 (35), 99.0% of seqs, 99.6% of clusters
#   Max mem  259Mb
#      Time  52.0s
#Throughput  684.2 seqs/sec.

#Clustering at 95% identity nearly all clusters are singletons, the max has 20 sequences

usearch -cluster_fast pooled_seqs_sorted.fasta -id 0.90 -centroids nr.fasta -uc clusters.uc
#35578 seqs, 35578 uniques, 35578 singletons (100.0%)
#00:00  78Mb Min size 1, median 1, max 1, avg 1.00
#00:00  74Mb done.
#00:59 219Mb  100.0% 25106 clusters, max size 927, avg 1.4
#00:59 219Mb  100.0% Writing centroids to nr.fasta        
                                                 
#      Seqs  35578 (35578 ()
#  Clusters  25106 (25106 (25106 (2)
#  Max size  927
#  Avg size  1.4
#  Min size  1
#Singletons  22961 (22961 (22961 (22961 (22), 64.5% of seqs, 91.5% of clusters
#   Max mem  227Mb
#      Time  59.0s
#Throughput  603.0 seqs/sec.

#Even at 90% identity (so not even near identical) 91.5% of clusters are singletons, the max has 927 sequences

#Above wasn't considering coverage. uclust is semi global alignment so doesn't even account for coverage
#If we consider coverage, and assign a conservative 80% query coverage, 91.8% of clusters are singletons, max 872 
usearch -cluster_fast pooled_seqs_sorted.fasta -id 0.90 -query_cov 0.80 -centroids nr.fasta -uc clusters.uc
#35578 seqs, 35578 uniques, 35578 singletons (100.0%)
#00:01  78Mb Min size 1, median 1, max 1, avg 1.00
#00:01  74Mb done.
#01:00 216Mb  100.0% 24471 clusters, max size 872, avg 1.5
#01:01 216Mb  100.0% Writing centroids to nr.fasta        
                                                 
#      Seqs  35578 (35578 ()
#  Clusters  24471 (24471 (24471 (2)
#  Max size  872
#  Avg size  1.5
#  Min size  1
#Singletons  22473 (22473 (22473 (22473 (22), 63.2% of seqs, 91.8% of clusters
#   Max mem  222Mb
#      Time  01:00
#Throughput  593.0 seqs/sec.

#Conclusion: There are no identical or near identical sequences within or across samples