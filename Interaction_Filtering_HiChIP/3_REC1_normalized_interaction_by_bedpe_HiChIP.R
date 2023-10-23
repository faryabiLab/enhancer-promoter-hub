#### --------------------------------------------------------------------------------------------------------------------------//
#### This is the third of 4 scripts used to create a list of valid HiChIP EE/PP/EP interactions as input for hub calling. 
#### Overview: This script calculates normalized contact frequency score for each loop in the input dataframe(s), which contains all 
# valid EE/PP/EP SMC1 HiChIP interactions for given condition(s). This script is currently written to evaluate PP/EE/EP hubs from 
# 2 conditions to demonstrate how it can be used to generate differential hub input, if desired. However, it can be easily modified
# to only generate output for 1 condition (see "optional" comments throughout).                                                                                                                              
#
#### Inputs: Datatable of SMC1 HiChIP PP/EP/EE interactions from each condition (output files from "2_make_filtered_EP_linkage.R"),
# .bed file of MACS-called, significant H3K27ac peaks that are NOT yet uniformly size-adjusted, and the number of rows in the 
# .bedpe Mango file (PETS >= 4) of all SMC1 HiChIP loops for normalization purposes. 
#
#### Dependencies: See R libraries below
#
#### Final Output: .rda file containing a dataframe of all valid SMC1 HiChIP PP/EE/EP loops with normalized contact frequency 
# (across 2 conditions) such that the loops can be filtered and then serve as input into the hub calling pipeline.
#### -------------------------------------------------------------------------------------------------------------------------//

#### Load dependencies
library("plyr")
library("dplyr")
library("cluster")
library("data.table")                                                                                                                   
library("stringr")

#### ----------------------------------------------------------------------------------------------------------//
#### Inputs  --------------------------------------------------------------------------------------------------//
#### ----------------------------------------------------------------------------------------------------------//

#### Specify desired working directory
mydir <- "/mnt/data1/yeqiao/analysis/Breast_cancer_notch/idea39_DMSO_GSI_linkage_REC1_500bp"
wd <- file.path(mydir, "/analysis_Fig3_EP")
setwd(wd)

mydir <- "/mnt/nas3/users/brent/analysis/Finalized_Hub_Scripts/Interaction_Filtering_HiChIP/example_input_data"
setwd(mydir)


#### Load .rda containing SMC1 HiChIP PP/EP/EE loops for condition 1 (output from "2_make_filtered_EP_linkage.R")
load("180408_REC1_SMC1_DMSO_linkage.rda")

#### (Optional) Load .rda containing SMC1 HiChIP PP/EP/EE loops for condition 2 (output from opt. second run of "2_make_filtered_EP_linkage.R")
load("180408_REC1_SMC1_GSI_linkage.rda")

#### [Optional] Read in enhancer positions for additional enhancer filtering by size (MACS-called, significant H3K27ac peaks, NOT uniformly size-adjusted)
# This file is in the precursor to the "*summits.bed" file and is in the format: "chr" "start" "end" "MACS_peak_id"
enh <- fread("150615_Rec1_H3K27ac_4hWO_merged_sorted-MACS2-FDR_1e-5-SHIFT_85_peaks.bed", header = F, stringsAsFactors = F)

#### Use number of rows in original (pre-filtered) SMC1 HiChIP loop .bedpe file for normalization
dmso.bd <- 318292691 #rows in .bedpe for condition 1
gsi.bd <- 356725566 # Optional: rows in .bedpe for condition 2
scale <- 1E8 

# Define output file name 
out_name <- "180612_REC1_smcLoop_DMSO_GSI_221108_BP_clique_input.rda"

#### ----------------------------------------------------------------------------------------------------------//
#### Main Script -------------------------------------------------------------------------------------//
#### ----------------------------------------------------------------------------------------------------------//

#### Create comprehensive, annotated table of all PP/EE/EP interactions
DMSO.keep <- data.table(subset(DMSO.anchor.keep, (grepl("ENSG", p.int.id) & grepl("ENSG", q.int.id) & linkage.type == "PP") | 
                      (grepl("MACS", p.int.id) & grepl("MACS", q.int.id) & linkage.type == "EE") |  
                      ((grepl("ENSG", p.int.id) | grepl("ENSG", q.int.id)) & linkage.type == "PE"))) 

#### Create unique ids for each enhancer and promoter in format of "chrX_0000000000_000000000"
DMSO.keep$p_id <- paste(DMSO.keep$p_chr, DMSO.keep$p.int.start, DMSO.keep$p.int.end, sep="_")
DMSO.keep$q_id <- paste(DMSO.keep$q_chr, DMSO.keep$q.int.start, DMSO.keep$q.int.end, sep="_")

#### Calculate normalized interaction score / contact frequency for loops in condition 1 
DMSO.score <- DMSO.keep[,.(sum(score), sum(score)*1E8/dmso.bd), by = .(p.int.id, p_id, q.int.id, q_id, linkage.type)]
colnames(DMSO.score)[c(6,7)] <- c("REC1_loop_dmso_count", "REC1_loop_dmso_norm")
DMSO.score <- data.frame(DMSO.score)

#### [Optional-- Cond2] Process condition 2 (i.e. GSI) using same code as above
GSI.keep <- data.table(subset(GSI.anchor.keep, (grepl("ENSG", p.int.id) & grepl("ENSG", q.int.id) & linkage.type == "PP") | 
                      (grepl("MACS", p.int.id) & grepl("MACS", q.int.id) & linkage.type == "EE") |  
                      ((grepl("ENSG", p.int.id) | grepl("ENSG", q.int.id)) & linkage.type == "PE")))

#### [Optional-- Cond2] Create unique ids for each enhancer and promoter in format of "chrX_0000000000_000000000"
GSI.keep$p_id <- paste(GSI.keep$p_chr, GSI.keep$p.int.start, GSI.keep$p.int.end, sep="_")
GSI.keep$q_id <- paste(GSI.keep$q_chr, GSI.keep$q.int.start, GSI.keep$q.int.end, sep="_")

#### [Optional-- Cond2] Calculate normalized interaction score / contact frequency for loops in condition 2 
GSI.score <- GSI.keep[,.(sum(score), sum(score)*1E8/gsi.bd), by = .(p.int.id, p_id, q.int.id, q_id, linkage.type)]
colnames(GSI.score)[c(6,7)] <- c("REC1_loop_gsi_count", "REC1_loop_gsi_norm")
GSI.score <- data.frame(GSI.score)

#### [Optional-- Cond2] Combine condition 1 and 2 dataframes for easier analysis
DMSO.GSI.raw <- data.table(merge(DMSO.score, GSI.score, by = c("p.int.id", "q.int.id", "linkage.type"), all.x = T, all.y = T))
DMSO.GSI.raw[is.na(DMSO.GSI.raw)] <- 0

#### [Extra filtering] Calculate enhancer size (from before size adjustment) to ensure that only valid enhancers 
# that are >500bp in length are kept, V1 = chrom, V2 = start, V3 = end, V4 = MACS_peak_id
enh[, `:=`(size = V3 - V2, summit = gsub("peak", "summit", V4))]
enh.keep <- enh[size > 500,]

#### [Extra filtering] Filter interactions with final list of valid enhancers to ensure that only the relevant 
# interactions are carried forward and make final dataframe
DMSO.GSI.01 <- DMSO.GSI.raw[(linkage.type == "PP") ,]
DMSO.GSI.02 <- DMSO.GSI.raw[(linkage.type == "PE" & (p.int.id%in%enh.keep$summit | q.int.id%in%enh.keep$summit)), ]
DMSO.GSI.03 <- DMSO.GSI.raw[(linkage.type == "EE" & p.int.id%in%enh.keep$summit & q.int.id%in%enh.keep$summit), ]

DMSO.GSI.filt <- rbind(DMSO.GSI.01, DMSO.GSI.02, DMSO.GSI.03)
DMSO.GSI.filt <- data.frame(DMSO.GSI.filt)

#### Mutate dataframe to get final table for input into interaction filtration script
#### annotate each p and q element
DMSO.GSI.filt2 <- mutate(DMSO.GSI.filt, p_annot = case_when(
  grepl("ENSG", p.int.id) == TRUE ~ "promoter", 
  grepl("MACS", p.int.id) == TRUE ~ "enhancer",
  TRUE ~ 'NA'))

DMSO.GSI.filt3 <- mutate(DMSO.GSI.filt2, q_annot = case_when(
  grepl("ENSG", q.int.id) == TRUE ~ "promoter", 
  grepl("MACS", q.int.id) == TRUE ~ "enhancer",
  TRUE ~ 'NA'))

#### Remove unnecessary columns from dataframe
DMSO.GSI.filt.clique = subset(DMSO.GSI.filt3, select=c(p_id.x, q_id.x, p_id.y, q_id.y, linkage.type, REC1_loop_dmso_norm, REC1_loop_gsi_norm, p_annot, q_annot))

#### Assign NA to 0 value promoter/enhancer_id cols in order to remove them
DMSO.GSI.filt.clique[DMSO.GSI.filt.clique == 0] <- NA

#### Coalesce columns 
DMSO.GSI.filt.clique$p_id <- coalesce(DMSO.GSI.filt.clique$p_id.x,DMSO.GSI.filt.clique$p_id.y)
DMSO.GSI.filt.clique$q_id <- coalesce(DMSO.GSI.filt.clique$q_id.x,DMSO.GSI.filt.clique$q_id.y)

#### Remove unnecessary columns again and make NA --> 0
REC_DMSO.GSI.filt.clique = subset(DMSO.GSI.filt.clique, select=c(p_id, q_id, linkage.type, REC1_loop_dmso_norm, REC1_loop_gsi_norm, p_annot, q_annot))
REC_DMSO.GSI.filt.clique[is.na(REC_DMSO.GSI.filt.clique)] <- 0

# Save final dataframe for final interaction filtration step in subsequent script
save(REC_DMSO.GSI.filt.clique, file = out_name)

