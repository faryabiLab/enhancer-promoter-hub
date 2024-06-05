#!/bin/bash  
#### --------------------------------------------------------------------------------------------------------------------------//
#### This is the first of 4 scripts used to create a list of valid HiC EE/PP/EP interactions as input for hub calling. 
#### Overview: This script ultimately creates a .bed file containing the union of all accessible enhancer and promoter
# anchors across 2 conditions (i.e. Ib-sensitive and Ib-resistant Rec-1), where accessible 
# enhancers are determined from the intersection of H3K27ac ChIP-seq peaks and ATAC-seq peaks and accessible 
# promoters are determined from the intersection of known hg19 transcription start sites (TSS) and ATAC-seq peaks.                                                                                                                                  
#
#### Inputs: hg19 assembly information, .bed file of significant MACS-filtered H3K27ac peaks, .bed file of significant 
# MACS-filtered ATAC-seq peaks, BED4 file of all hg19 TSS, integer representing 1/2 length of desired size (in bp) for all 
# ATAC, TSS, and H3K27ac peaks
#
#### Dependency: bedtools suite
#
#### Intermediate Outputs: ATAC/H3K27ac summit .bed4 files containing centers of ATAC/H3K27ac peaks (respectively), 
# ATAC/H3K27ac/TSS ext .bed4 file containing size-adjusted, centered ATAC/H3K27ac/TSS peaks (all 5kbp in size, respectively)
#
#### Final Output: .bed (8) file containing intervals of all size-adjusted, accessible, enhancer/promoter anchors
#### -------------------------------------------------------------------------------------------------------------------------//

# Load modules 
module load bedtools-2.29.1
module load ucsc369

#### Set directory
WKDIR=test_HiC
DATADIR=Interaction_Filtering_HiC/example_input_data
cd ${WKDIR}

#### Inputs: Load file path containing hg19 chromosome intervals
GENOME=${DATADIR}/GRCh37.75_ChromInfo.txt

#### Inputs: Load path to .bed file containing union of ATAC-seq peaks (from MACS peak calling) from Ib-sens, Ib-res, & Ib-sens treated cells and name output file prefix
ATAC=REC1_DMSO_IBR_RES_ATAC_merged_id.bed
ATAC_BASE=REC1_ATAC

#### Inputs: Load path to file containing all hg19 Transcription Start Sites (TSS) and name output file prefix
TSS=170709_hg19_ensembl_ensg_longest_TSS.txt
TSS_BASE=hg19_ensembl_TSS

#### Inputs: Load path to .bed file containing union of H3K27ac-seq peaks (from MACS peak calling) from Ib-sens, Ib-res, & Ib-sens treated cells and name output file prefix
ENH=REC1_DMSO_IBR_RES_H3K27ac_merged_id.bed
ENH_BASE=REC1_H3K27ac

#### Inputs: Name of final output file with intersected, size-adjusted peaks
OUT=230306_REC1_TSS_ENH_ATAC_intersect.bed

#### Inputs: Initialize variable containing # base pairs to adjust/extend all peaks
EXT=2500

#### Sort ATAC peaks and generate .bed of just the center/summit of each peak interval 
ATAC_SUMMIT=${ATAC_BASE}_summit.txt
bedSort ${DATADIR}/${ATAC} ${ATAC}
awk -F "\t" 'BEGIN{OFS=FS}{print $1,$2,$3,int(($3-$2)/2),$4}' ${ATAC} | awk -F "\t" 'BEGIN{OFS=FS}{print $1,$2+$4,$2+$4+1,$5}' > ${ATAC_SUMMIT}

#### Extend the ATAC peak summits by $EXT both upstream and downstream
ATACEXT=${ATAC_BASE}_summit_ext_${EXT}.txt
bedtools slop -i ${ATAC_SUMMIT} -b ${EXT} -g ${GENOME} > ${ATACEXT}

#### Extend the TSS peak summits by $EXT both upstream and downstream
TSSEXT=${TSS_BASE}_${EXT}.bed
bedtools slop -i ${DATADIR}/${TSS} -b ${EXT} -g ${GENOME} > ${TSSEXT}

#### Sort H3K27ac peaks and generate .bed of just the center/summit of each peak interval 
ENH_SUMMIT=${ENH_BASE}_summit.txt
bedSort ${DATADIR}/${ENH} ${ENH}
cut -f 1-3 ${ENH} | awk  -F "\t" 'BEGIN{OFS=FS}{print $0,int(($3-$2)/2)}' | awk -F "\t" 'BEGIN{OFS=FS}{print $1,$2+$4,$2+$4+1,$1"_"$2"_"$3}' > ${ENH_SUMMIT}

#### Extend the H3K27ac peak summits by $EXT both upstream and downstream
ENHEXT=${ENH_BASE}_summit_ext_${EXT}.bed
bedtools slop -i ${ENH_SUMMIT} -b ${EXT} -g ${GENOME} > ${ENHEXT}

#### Output: Intersect intermediate files to create final output .bed file containing all accessible, enhancer and promoter anchors where each
# anchor is annotated with a unique identifier corresponding to its size-adjusted ATAC peak
cat ${TSSEXT} ${ENHEXT} | bedtools sort -i stdin | bedtools intersect -a stdin -b ${ATACEXT} -wa -wb > ${OUT}
