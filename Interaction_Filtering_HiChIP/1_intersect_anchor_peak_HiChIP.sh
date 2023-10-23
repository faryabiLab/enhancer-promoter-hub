#!/bin/bash  
#### --------------------------------------------------------------------------------------------------------------------------//
#### This is the first of 4 scripts used to create a list of valid HiChIP EE/PP/EP interactions as input for hub calling. 
#### Overview: This script ultimately creates two .bed files containing the coordinates for all enhancer and promoter loop anchors 
# from a HiChIP experiment (e.g. Rec-1), where enhancers are determined from H3K27ac ChIP-seq peaks
# and promoters are determined from known hg19 transcription start sites (TSS)                                                                                                                                   
#
#### Inputs: hg19 assembly information, .bed4 file of MACS-filtered H3K27ac peak summits (i.e. center of peak), .bedpe file of Mango pipeline
# output for significant SMC1 HiChIP loops (PETS >= 4), .bed4 file of all hg19 TSS, integers representing the length of desired size adjustment for all TSS, HiChIP, and H3K27ac peaks
#
#### Dependency: bedtools suite
#
#### Intermediate Outputs: H3K27ac summit .bed files containing centers of H3K27ac peaks (respectively), 
# H3K27ac/TSS ext .bed file containing size-adjusted, centered H3K27ac/TSS peaks (10kb & 5kb in size, respectively), .bed file containing 
# union of all H3K27ac and TSS size-adjusted, centered peaks, and .bed file containing the intersection of all H3K27ac and TSS size-adjusted, centered peaks
#
#### Final Output:.bed file containing intervals of all accessible, enhancer and promoter loop anchors (from Rec-1 SMC1 HiChIP)
#### -------------------------------------------------------------------------------------------------------------------------//

#### Load server modules
module() { eval `/usr/bin/modulecmd bash $*`; }
module use /mnt/data1/modulefiles
module load bedtools-2.25.0

#### Set directories
WORKDIR="/mnt/nas3/users/brent/analysis/Finalized_Hub_Scripts/Interaction_Filtering_HiChIP/example_input_data"
cd $WORKDIR

#### Inputs: Load files containing hg19 chromosome intervals, H3K27ac ChIP-seq summits, SMC1 HiChIP significant loops, and hg19 TSS
GENOME=GRCh37.75_ChromInfo.txt # format: "chrom" "length"
CHIPEAK=150615_Rec1_H3K27ac_4hWO_merged_sorted-MACS2-FDR_1e-5-SHIFT_85_summits.bed # center of H3K27ac peaks, format: "chrom" "start" "end" "peak_id"
SIGLOOP=s99_180201_Rec1_DMSO_Smc1_FDR_0.05_PETS_4.mango # Mango output of significant SMC1 HiChIP loops (>4 PETs), format: "chrom1" "start1" "end1"  "chrom2" "start2" "end2" "pair_id" "length" "peak1_id" "peak2_id" "PETs"... 
TSS=170709_hg19_ensembl_ensg_longest_TSS.txt # format: "chrom" "start" "end" "ENSG"

#### Inputs: Initialize variable containing # base pairs to adjust/extend all peaks
EXT=5000
TSSLOP=2500

#### Inputs: base names of output files
CHIP_BASE=REC1_WO_H3K27ac
TSS_BASE=hg19_ensembl_ensg_TSS
OUT_BASE=180408_REC1

#### [Optional] Sort chip-seq peaks and TSS
#bedSort ${CHIPEAK} ${CHIPEAK} 
#bedSort ${TSS} ${TSS}

#### Adjust size of H3K27ac ChIP-seq peaks
# Amount of size adjustment (i.e. 5kb vs. 10kb can be chosen based on data resolution/quality)
CHIPEXT=${CHIP_BASE}_10k.bed
cut -f 1-4 ${CHIPEAK} | bedtools slop -i stdin -b ${EXT} -g ${GENOME} > ${CHIPEXT}

#### Adjust size of TSS regions 
TSSEXT=${TSS_BASE}_5k.bed
bedtools slop -i ${TSS} -b ${TSSLOP} -g ${GENOME} > ${TSSEXT}

#### Intermediate output: create union file of all size-adjusted H3K27ac ChIP-seq peaks and TSS intervals
cat ${TSSEXT} ${CHIPEXT} | sort -k1,1 -k2,2n > ${CHIP_BASE}_TSS_all.bed

#### Intermediate output: create file containing the intersection of all size-adjusted H3K27ac ChIP-seq peaks and TSS intervals
# This file will be used in the subsequent script to determine if an anchor coinciding with a H3K27ac peak AND TSS should be 
# annotated as an enhancer or promoter. 
bedtools intersect -a ${CHIPEXT} -b ${TSSEXT} -wa -wb > ${OUT_BASE}_peak_tss_intersect.bed

#### Divide SMC1 HiChIP loops (i.e. anchor pairs) into 2 separate files of p and q anchors, 
# find the center/summit of the SMC1 HiChIP anchor peaks, and adjust their size
cut -f 1,2,3,7,10 ${SIGLOOP} | awk  -F "\t" 'BEGIN{OFS=FS}{print $0,int(($3-$2)/2)}' | awk -F "\t" 'BEGIN{OFS=FS}{print $1,$2+$6,$2+$6+1,$4,$5}' | bedtools sort -i stdin | bedtools slop -i stdin -b ${EXT} -g ${GENOME} | sort -k 4 > p.temp.bed
cut -f 4,5,6,7,10 ${SIGLOOP} | awk  -F "\t" 'BEGIN{OFS=FS}{print $0,int(($3-$2)/2)}' | awk -F "\t" 'BEGIN{OFS=FS}{print $1,$2+$6,$2+$6+1,$4,$5}' | bedtools sort -i stdin | bedtools slop -i stdin -b ${EXT} -g ${GENOME} | sort -k 4 > q.temp.bed

#### Intersect size-adjusted SMC1 HiChIP loop anchors with union file of all size-adjusted H3K27ac and TSS peaks
# to annotate anchors for enhancers and promoters (respectively)
# header of output file: "anchor_chrom" "anchor_start" "anchor_end" "pair_id" "PETs" "peak_chrom" "peak_start" "peak_end" "peak_id (MACS* or ENSG*)"
bedtools intersect -a p.temp.bed -b ${CHIP_BASE}_TSS_all.bed -wa -wb > ${OUT_BASE}_anchor_intersect_p.bed
bedtools intersect -a q.temp.bed -b ${CHIP_BASE}_TSS_all.bed -wa -wb > ${OUT_BASE}_anchor_intersect_q.bed
