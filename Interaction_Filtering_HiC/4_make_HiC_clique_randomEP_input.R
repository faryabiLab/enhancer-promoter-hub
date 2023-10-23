#### ---------------------------------------------------------------------------------------------//
#### This is the fourth of 4 scripts used to create a list of valid HiC EE/PP/EP interactions as input for hub calling. 
#### Overview: This script filters the EE/PP/EP interactions for 2 given conditions (i.e. Rec-1 Ib-sensitive/Ib-resistant)
# using a given normalized interaction score cutoff to create a final list of valid regulatory element contacts for hub calling.
#
#### Inputs: 1) Dataframe containing EE/PP/EP interactions for condition 1 (Ib-sens) with normalized interaction scores (script 3 output), 
# 2) Dataframe containing EE/PP/EP interactions for condition 2 (Ib-res) with normalized interaction scores (script 3 output), and 
# 3) Double representing the normalized interaction score to consider an EE/EP/PP pair as a valid, true contact
#
#### Outputs: One .csv files of valid, filtered EE/PP/EP physical contacts for each condition (2 files total), which 
# can be directly used for subsequent command-line calling of hubs via hierarchical spectral clustering
# with format "enh/pro1", "enh/pro2", "normalized_loop_score_condX", "annot1", "annot2"
#### ---------------------------------------------------------------------------------------------//

#### Load library and numeric preferences
library("dplyr")
options(scipen=999) # Convert sci notation to decimal if necessary

#### Clear local environment and set workdir
rm(list = ls())
setwd("/mnt/data0/brent/3D_Cliques/230306_REC1_Arima_IBR_Random_EP/cf3")

#### Inputs: Load file containing EE/PP/EP contacts with normalized interaction scores for each condition
load("/mnt/data0/brent/3D_Cliques/230306_REC1_Arima_IBR_Random_EP/clique_input/230306_REC1_DMSO_Arima_clique_input_random_EP.rda")
load("/mnt/data0/brent/3D_Cliques/230306_REC1_Arima_IBR_Random_EP/clique_input/230306_REC1_RES_Arima_clique_input_random_EP.rda")

#### Input: Define contact frequency cutoff for interaction filtering (3 for HiC contacts and 2 for SMC1 HiChIP contacts)
contact_cutoff = as.double('3')

#### Input: Define intermediate output file names/paths, out1 = control and out2 = case/experimental
out1 <- "REC1_DMSO_PPEEEP_cf3.csv"
out2 <- "REC1_RES_PPEEEP_cf3.csv"

#### Filter condition 1 (Rec-1 Ib-sens) PPEEEP interactions and clean up dataframe
Control_Clique_PPEEEP <- Control_HiC_clique_input[Control_HiC_clique_input$HiC_over_random_EP_norm_2 > contact_cutoff,]
Control_Clique_PPEEEP = subset(Control_Clique_PPEEEP, select=c("p_int_id", "q_int_id", "HiC_over_random_EP_norm_2", "p_annotation", "q_annotation"))
colnames(Control_Clique_PPEEEP) <- c("item_id1", "item_id2", "loop_DMSO_norm", "item_annot1", "item_annot2")
Control_Clique_PPEEEP <- unique(Control_Clique_PPEEEP)

#### Filter condition 2 (Rec-1 Ib-res) PPEEEP interactions and clean up dataframe
Case_Clique_PPEEEP <- Case_HiC_clique_input[Case_HiC_clique_input$HiC_over_random_EP_norm_2 > contact_cutoff,]
Case_Clique_PPEEEP = subset(Case_Clique_PPEEEP, select=c("p_int_id", "q_int_id", "HiC_over_random_EP_norm_2", "p_annotation", "q_annotation"))
colnames(Case_Clique_PPEEEP) <- c("item_id1", "item_id2", "loop_RES_norm", "item_annot1", "item_annot2")
Case_Clique_PPEEEP <- unique(Case_Clique_PPEEEP)

#### Save both dataframes to separate files as inputs to final helper bash shell script 
write.table(Control_Clique_PPEEEP, file = out1, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(Case_Clique_PPEEEP, file = out2, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

#### Call helper bash shell script to reorganize .csv files to create final dataframes of EE/PP/EP interactions for input into hub-calling pipeline 
# in format (chr_A_start1_stop1,chrA_start2_stop2,contactFreq) as well as enhancer/promoter annotation files for anchors 
cmd <- paste("./generate_ct_input.sh", out1, out2)
system(cmd)

