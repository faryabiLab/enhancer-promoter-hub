#### ---------------------------------------------------------------------------------------------//
#### This is the fourth of 4 scripts used to create a list of valid HiChIP EE/PP/EP interactions as input for hub calling. 
#### Overview: This script filters the EE/PP/EP interactions for 2 given conditions (i.e. Rec-1 sensitive/resistant)
# using a given normalized interaction score cutoff to create a final list of valid contacts between regulatory elements for hub calling.
# Note that this script can be easily modified to work with only 1 input condition, if desired. 
#
#### Dependencies: see library below
#
#### Inputs: 1) Dataframe containing EE/PP/EP interactions for condition 1 (DMSO) with normalized interaction scores, 
# 2) [optional] Dataframe containing EE/PP/EP interactions for condition 2 (GSI) with normalized interaction scores, and 
# 3) Double representing the normalized interaction score to consider an EE/EP/PP pair as a valid, true contact
#
#### Outputs: One .csv file of valid, filtered EE/PP/EP physical contacts for each condition, which 
# can be directly used for subsequent command-line calling of hubs via hierarchical spectral clustering
# with format "enh/pro1", "enh/pro2", "normalized_loop_score_condX", "annot1", "annot2"
#### ---------------------------------------------------------------------------------------------//

#### Load library and numeric preferences
library("dplyr")
options(scipen=999) # Convert sci notation to decimal if necessary

#### Clear local environment and set workdir
rm(list = ls())
setwd("/mnt/data0/brent/3D_Cliques/221109_REC1_GSI_cf2")

setwd("/mnt/nas3/users/brent/analysis/Finalized_Hub_Scripts/Interaction_Filtering_HiChIP/example_input_data")

#### Inputs: Load file containing EE/PP/EP contacts with normalized interaction scores for each condition (output of script 2)
load("180612_REC1_smcLoop_DMSO_GSI_221108_BP_clique_input.rda")

#### Input: Define contact frequency cutoff for interaction filtering (2 for SMC1 HiChIP contacts)
contact_cutoff = as.double('2')

#### Input: Define intermediate output file names/paths, out1 = control and (optional) out2 = case/experimental
out1 <- "REC1_DMSO_PPEEEP_cf2.csv"
out2 <- "REC1_GSI_PPEEEP_cf2.csv" #optional

# Filter condition 1 (i.e. DMSO) PPEEEP interactions and clean up dataframe
REC1_DMSO_Clique_PPEEEP <- REC_DMSO.GSI.filt.clique[REC_DMSO.GSI.filt.clique$REC1_loop_dmso_norm > contact_cutoff,]
REC1_DMSO_Clique_PPEEEP = subset(REC1_DMSO_Clique_PPEEEP, select=c(q_id, p_id, REC1_loop_dmso_norm, q_annot, p_annot))
colnames(REC1_DMSO_Clique_PPEEEP) <- c("item_id1", "item_id2", "loop_DMSO_norm", "item_annot1", "item_annot2")
REC1_DMSO_Clique_PPEEEP <- unique(REC1_DMSO_Clique_PPEEEP)

# [Optional] Filter condition 2 (i.e. GSI) PPEEEP interactions and clean up dataframe
REC1_GSI_Clique_PPEEEP <- REC_DMSO.GSI.filt.clique[REC_DMSO.GSI.filt.clique$REC1_loop_gsi_norm > contact_cutoff,]
REC1_GSI_Clique_PPEEEP = subset(REC1_GSI_Clique_PPEEEP, select=c(q_id, p_id, REC1_loop_gsi_norm, q_annot, p_annot))
colnames(REC1_GSI_Clique_PPEEEP) <- c("item_id1", "item_id2", "loop_GSI_norm", "item_annot1", "item_annot2")
REC1_GSI_Clique_PPEEEP <- unique(REC1_GSI_Clique_PPEEEP)

# Save datatables to separate files as inputs to final helper bash shell script 
write.table(REC1_DMSO_Clique_PPEEEP, file = out1, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(REC1_GSI_Clique_PPEEEP, file = out2, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

#### Call helper bash shell script to reorganize .csv files to create final dataframes of EE/PP/EP interactions for input into hub-calling pipeline 
# in format (chr_A_start1_stop1,chrA_start2_stop2,contactFreq) as well as enhancer/promoter annotation files for anchors 
cmd <- paste("../generate_ct_input.sh", out1, out2)
system(cmd)

