#!/bin/bash
#### ---------------------------------------------------------------------------------------------//
#### This is a helper script to "quantify_bedpe_HiC.R".
#### Overview: This script takes a single .bedpe input file of accessible HiC interactions and creates two datatables of interaction anchors
# and annotates each of these anchors with unique ATAC-seq identifiers in two separate files (i.e. 'p' and 'q' anchors)
# 
#### Args (passed from command line call): 1) .bedpe file of all (Rec-1 Ib-sensitive) HiC reads including natural duplicates,
# 2) .bed file corresponding to (REC-1) size-adjusted ATAC-seq peaks created as intermediate output by "1_intersect_anchor_peak.sh"
#
#### Dependencies: bedtools suite
#
#### Final Outputs: 2 BED-like files containing all observed, accessible enhancer/promoter anchors annotated with unique id 
# consisting of size-adjusted ATAC peak interval with 
# columns "chrom[HiC]" "start[HiC]" "end[HiC]" "uniq_read_id" "chrom[ATAC]" "start[ATAC]" "end[ATAC]" "ATAC_peak_id"
#### ---------------------------------------------------------------------------------------------//

#### Load modules/dependencies
module load bedtools-2.29.1

####  Create two .bed4 files from .bedpe of (Rec-1 Ib-sensitive) HiC interactions for 'p' and 'q' anchors (respectively)
awk ' {print $1"\t"$2"\t"$3"\t"$7}' $1 > HiC_p.bed 
awk ' {print $4"\t"$5"\t"$6"\t"$7}' $1 > HiC_q.bed 

# Intersect each of the two anchor files separately with the .bed file corresponding to 5kb size-adjusted ATAC-seq
# peaks created as intermediate output by "intersect_anchor_peak.sh"
bedtools intersect -a $2 -b HiC_p.bed -wa -wb > HiC_ATAC_2500_p.bed  
bedtools intersect -a $2 -b HiC_q.bed -wa -wb > HiC_ATAC_2500_q.bed 

