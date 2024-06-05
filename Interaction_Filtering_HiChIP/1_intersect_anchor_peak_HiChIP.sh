#!/bin/bash  
#### --------------------------------------------------------------------------------------------------------------------------//
#### This is the first of 4 scripts used to create a list of valid HiChIP EE/PP/EP interactions as input for hub calling. 
#### Overview: This script ultimately creates two .bed files containing the coordinates for all enhancer and promoter interaction anchors 
# from a HiChIP experiment (e.g. MB157), where enhancers are determined from H3K27ac ChIP-seq peaks
# and promoters are determined from known hg19 transcription start sites (TSS)                                                                                                                                   
#
#### Inputs: hg19 assembly information, .bed4 file of MACS-filtered H3K27ac peaks/summits, .bedpe file of FitHiChIP pipeline
# output for significant SMC1 HiChIP interactions (q < 0.05), .bed4 file of all hg19 TSS, integers representing the length of 
# desired size adjustment for all TSS, HiChIP, and H3K27ac peaks
#
#### Dependency: bedtools suite
#
#### Intermediate Outputs: H3K27ac/TSS ext .bed file containing size-adjusted, centered H3K27ac/TSS peaks (10kb & 5kb in size, respectively), .bed file containing 
# union of all H3K27ac and TSS size-adjusted, centered peaks, and .bed file containing the intersection of all H3K27ac and TSS size-adjusted, centered peaks
#
#### Final Output:.bed file containing intervals of all accessible, enhancer and promoter interaction anchors (from MB157 SMC1 HiChIP)
#### -------------------------------------------------------------------------------------------------------------------------//

# Load modules 
module load bedtools-2.29.1
module load ucsc369

#### Set directory
WKDIR=test_HiChIP
DATADIR=Interaction_Filtering_HiChIP/example_input_data
mkdir ${WKDIR}
cd ${WKDIR}

#### Inputs: Load file path containing hg19 chromosome intervals
GENOME=${DATADIR}/GRCh37.75_ChromInfo.txt

#### Inputs: Load path to file containing all hg19 Transcription Start Sites (TSS) and name output file prefix
TSS=170709_hg19_ensembl_ensg_longest_TSS.txt
TSS_BASE=hg19_ensembl_TSS

#### Inputs: Load path to .bed file containing union of H3K27ac-seq peaks (from MACS peak calling)
#ENH=MB157_WO_H3K27ac_merged_peaks.bed
ENH_SUMMIT=${DATADIR}/MB157_WO_H3K27ac_merged_summit.bed # supply here if already created by previous pipeline 
ENH_BASE=MB157_H3K27ac

#### Inputs: Initialize variables containing # base pairs to adjust/extend peaks
EXT=5000
TSSLOP=2500

#### Inputs: Load path to file containing all significant FitHiChIP interactions (output from FitHiChIP)
SIGLOOP=MB157_DMSO_SMC1-HiChIP.interactions_FitHiC_0.05Q.bed
SIGLOOP2=MB157_DMSO_SMC1-HiChIP.interactions_FitHiC_0.05Q_id.bed

#### Create unique identifiers for each FitHiChIP interaction
awk  -F "\t" 'BEGIN{OFS=FS}{print $0,$1"_"$2"_"$3"_"$4"_"$5"_"$6}' ${DATADIR}/${SIGLOOP} > ${SIGLOOP2}

#### Extend the TSS peak summits by $TSSLOP both upstream and downstream
bedSort ${DATADIR}/${TSS} ${TSS}
TSSEXT=${TSS_BASE}_${TSSLOP}.bed
bedtools slop -i ${TSS} -b ${TSSLOP} -g ${GENOME} > ${TSSEXT}

#### If necessary, generate .bed file of the center/summit of each peak interval
# bedSort ${DATADIR}/${ENH} ${ENH}
# ENH_SUMMIT=${ENH_BASE}_summit.txt
# cut -f 1-3 ${ENH} | awk  -F "\t" 'BEGIN{OFS=FS}{print $0,int(($3-$2)/2)}' | awk -F "\t" 'BEGIN{OFS=FS}{print $1,$2+$4,$2+$4+1,$1"_"$2"_"$3}' > ${ENH_SUMMIT}

#### Extend the H3K27ac peak summits by $EXT both upstream and downstream
ENHEXT=${ENH_BASE}_summit_ext_${EXT}.bed
bedtools slop -i ${ENH_SUMMIT} -b ${EXT} -g ${GENOME} > ${ENHEXT}

#### Intermediate output: create union file of all size-adjusted H3K27ac ChIP-seq peaks and TSS intervals
cat ${TSSEXT} ${ENHEXT} | sort -k1,1 -k2,2n > MB157_peak_tss_all.bed

#### Intermediate output: create file containing the intersection of all size-adjusted H3K27ac ChIP-seq peaks and TSS intervals
# This file will be used in the subsequent script to determine if an anchor coinciding with a H3K27ac peak AND TSS should be 
# annotated as an enhancer or promoter. 
bedtools intersect -a ${ENHEXT} -b ${TSSEXT} -wa -wb > MB157_peak_tss_intersect.bed

#### Divide SMC1 HiChIP interactions into 2 separate files of p and q anchors, 
# find the center/summit of the SMC1 HiChIP anchor peaks, and adjust their size
cut -f 1,2,3,25,26,27 ${SIGLOOP2} | awk  -F "\t" 'BEGIN{OFS=FS}{print $0,int(($3-$2)/2)}' | awk -F "\t" 'BEGIN{OFS=FS}{print $1,$2+$7,$2+$7+1,$6,$4,$5}' | bedtools sort -i stdin | bedtools slop -i stdin -b ${EXT} -g ${GENOME} | sort -k 4 > p.temp.bed
cut -f 4,5,6,25,26,27 ${SIGLOOP2} | awk  -F "\t" 'BEGIN{OFS=FS}{print $0,int(($3-$2)/2)}' | awk -F "\t" 'BEGIN{OFS=FS}{print $1,$2+$7,$2+$7+1,$6,$4,$5}' | bedtools sort -i stdin | bedtools slop -i stdin -b ${EXT} -g ${GENOME} | sort -k 4 > q.temp.bed

#### Intersect size-adjusted SMC1 HiChIP anchors with union file of all size-adjusted H3K27ac and TSS peaks
# to annotate anchors for enhancers and promoters (respectively)
# header of output file: "anchor_chrom" "anchor_start" "anchor_end" "pair_id" "p_value" "q_value" "peak_chrom" "peak_start" "peak_end" "peak/TSS_id (MACS* or ENSG*)"
bedtools intersect -a p.temp.bed -b MB157_peak_tss_all.bed -wa -wb > MB157_anchor_intersect_p_0.05Q.bed
bedtools intersect -a q.temp.bed -b MB157_peak_tss_all.bed -wa -wb > MB157_anchor_intersect_q_0.05Q.bed

