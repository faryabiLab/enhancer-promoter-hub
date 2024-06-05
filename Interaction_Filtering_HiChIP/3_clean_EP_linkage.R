#### --------------------------------------------------------------------------------------------------------------------------//
#### This is the third of 4 scripts used to create a list of valid HiChIP EE/PP/EP interactions as input for hub calling. 
#### Overview: This script is used to perform further filtering (if desired) and clean up the dataframe of spatial interactions
# between regulatory elements to facilitate the final step of creating hub calling input.                                                                                                                        
#
#### Inputs: Datatable of SMC1 HiChIP PP/EP/EE interactions for the given condition (output from "2_make_filtered_EP_linkage.R"),
# (optional) .bed file of H3K27ac peaks, (optional) .bed file of H3K27ac summits, (optional) .bed file of expressed genes with TSS
#
#### Dependencies: See R libraries below
#
#### Final Output: .rda file containing a dataframe of all valid SMC1 HiChIP PP/EE/EP interactions in the given condition.
#### -------------------------------------------------------------------------------------------------------------------------//

#### Load dependencies
rm(list = ls())
library("plyr")
library("dplyr")
library("cluster")
library("data.table")                                                                                                                   
library("stringr")

#### ----------------------------------------------------------------------------------------------------------//
#### Inputs  --------------------------------------------------------------------------------------------------//
#### ----------------------------------------------------------------------------------------------------------//

#### Specify desired working directory
datadir <- "Interaction_Filtering_HiChIP/example_input_data"
wkdir <- "test_HiChIP"
setwd(wkdir)

#### Load .rda containing SMC1 HiChIP PP/EP/EE interactions (output from "2_make_filtered_EP_linkage.R")
load("MB157_SMC1_DMSO_linkage.rda")

#### Define output file name 
out_name <- "MB157_SMC1_regulatory_int.rda"

#### [Optional] Read in enhancer positions for additional enhancer filtering by size (MACS-called, significant H3K27ac peaks, NOT uniformly size-adjusted)
enh <- fread(paste(datadir, "MB157_WO_H3K27ac_merged_peaks.bed", sep = "/"), header = F, stringsAsFactors = F)

#### [Optional] Load expressed gene TSS to correct for few overlapping TSSes
load(paste(datadir, "MB157_expressed_TSS.rda", sep = "/"))

#### [Optional] Load H3K27ac summits to correct for few overlapping TSSes
enh_summit <- fread(paste(datadir, "MB157_WO_H3K27ac_merged_summit.bed", sep = "/"), header = F, stringsAsFactors = F)


#### ----------------------------------------------------------------------------------------------------------//
#### Main Script -------------------------------------------------------------------------------------//
#### ----------------------------------------------------------------------------------------------------------//

#### Create comprehensive, annotated table of all PP/EE/EP interactions
# Eliminate duplicate entries, including those from anchors overlapping with expressed promoters with H3K27ac peaks
DMSO.keep <- data.table(subset(DMSO.anchor.keep, (grepl("ENSG", p.int.id) & grepl("ENSG", q.int.id) & linkage.type == "PP") | 
                                 (grepl("MACS", p.int.id) & grepl("MACS", q.int.id) & linkage.type == "EE") |  
                                 ((grepl("ENSG", p.int.id) | grepl("ENSG", q.int.id)) & linkage.type == "PE"))) 

#### Create unique ids for each enhancer and promoter in format of "chrX_0000000000_000000000"
DMSO.keep$p_id <- paste(DMSO.keep$p_chr, DMSO.keep$p.int.start, DMSO.keep$p.int.end, sep="_")
DMSO.keep$q_id <- paste(DMSO.keep$q_chr, DMSO.keep$q.int.start, DMSO.keep$q.int.end, sep="_")

#### [Optional] Characterize each interaction by its smallest FitHiChIP p/q values (facilitates further filtering, if desired)
DMSO.score <- DMSO.keep[,.(min(pval), min(qval)), by = .(p.int.id, p_id, q.int.id, q_id, linkage.type)] #original for hub input
colnames(DMSO.score)[c(6,7)] <- c("MB157_loop_dmso_count", "MB157_loop_dmso_norm")
DMSO.score <- data.frame(DMSO.score)

#### [Extra filtering] Calculate enhancer size (from before size adjustment) to ensure that only valid enhancers 
# that are >500bp in length are kept, V1 = chrom, V2 = start, V3 = end, V4 = MACS_peak_id
enh[, `:=`(size = V3 - V2, summit = gsub("peak", "summit", V4))]
enh.keep <- enh[size > 500,]

#### [Extra filtering] Filter interactions with final list of valid enhancers to ensure that only the relevant 
# interactions are carried forward and make final dataframe
DMSO.01 <- DMSO.score[(DMSO.score$linkage.type == "PP"), ]
DMSO.02 <- DMSO.score[(DMSO.score$linkage.type == "PE" & (DMSO.score$p.int.id%in%enh.keep$summit | DMSO.score$q.int.id%in%enh.keep$summit)), ]
DMSO.03 <- DMSO.score[(DMSO.score$linkage.type == "EE" & (DMSO.score$p.int.id%in%enh.keep$summit & DMSO.score$q.int.id%in%enh.keep$summit)), ]

DMSO.filt <- rbind(DMSO.01, DMSO.02, DMSO.03)
DMSO.filt <- data.frame(DMSO.filt)

#### ----------------------------------------------------------------------------------------------------------//
#### keep proper EE/EP/PP -------------------------------------------------------------------------------------//
#### ----------------------------------------------------------------------------------------------------------//

#### [Optional filtering] Correct for overlapping TSSes, not strictly necessary for downstream analysis 
ss <- MB157_expressed_TSS
s1 <- ss[duplicated(ss$tss), ]
s1$tss <- s1$tss + 1
pos1 <- rbind(ss[!duplicated(ss$tss), ], s1)
pos2 <- unique(enh_summit[,c("V4", "V2", "V1")])
colnames(pos2) <- colnames(pos1)
ep.pos <- rbind(pos1, pos2)

#### [Optional filtering] Update interaction positions based on overlapping TSS correction
DMSO.filt <- merge(DMSO.filt, ep.pos, by.x = "p.int.id", by.y = "id", all.x = T)
colnames(DMSO.filt)[c(ncol(DMSO.filt) - 1, ncol(DMSO.filt))] <- c("p.start", "p.chr")
DMSO.filt <- merge(DMSO.filt, ep.pos, by.x = "q.int.id", by.y = "id", all.x = T)
colnames(DMSO.filt)[c(ncol(DMSO.filt) - 1, ncol(DMSO.filt))] <- c("q.start", "q.chr")
DMSO.filt <- data.frame(DMSO.filt)

#### Mutate dataframe to get final table for input into interaction filtration script
#### annotate each p and q element
DMSO.filt2 <- mutate(DMSO.filt, p_annot = case_when(
    grepl("ENSG", p.int.id) == TRUE ~ "promoter", 
    grepl("MACS", p.int.id) == TRUE ~ "enhancer",
    TRUE ~ 'NA'))

DMSO.filt3 <- mutate(DMSO.filt2, q_annot = case_when(
  grepl("ENSG", q.int.id) == TRUE ~ "promoter", 
  grepl("MACS", q.int.id) == TRUE ~ "enhancer",
  TRUE ~ 'NA'))

#### Remove unnecessary columns from dataframe
DMSO.filt.clique = subset(DMSO.filt3, select=c(p_id, q_id, linkage.type, MB157_loop_dmso_norm, p_annot, q_annot))

#### Assign NA to 0 value promoter/enhancer_id cols in order to remove them
DMSO.filt.clique[DMSO.filt.clique == 0] <- NA

#### [Optional, not illustrated] Generate unique score for each interaction (e.g., based on read count or statistical significance)
# if desired, to further filter dataframe before hub-calling and/or to influence interaction clustering

#### Save dataframe of significant regulatory spatial interactions for final step in subsequent script
save(DMSO.filt.clique, file = out_name)

