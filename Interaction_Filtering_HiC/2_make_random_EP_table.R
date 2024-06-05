#### ---------------------------------------------------------------------------------------------//
#### This is the second of 4 scripts used to create a list of valid HiC EE/PP/EP interactions as input for hub calling. 
#### Overview: This script creates a datatable of all possible combinations of accessible enhancer-
# enhancer (EE), enhancer-(active) promoter (EP), and (active) promoter-promoter (PP) interactions across the given
# conditions (i.e. Rec-1 Ib-sens/Ib-res cells) that will allow identification and normalization of 
# EE/EP/PP spatial interactions in the subsequent script to create hub calling input. 
# 
#### Inputs: 1) output .bed8 file of "1_intersect_anchor_peak.sh" containing intervals of all accessible, 
# enhancers and promoters from Ib-sens/Ib-res Rec-1, 2) .bed4 file containing the union of expressed 
# hg19 genes within the given conditions, 3) integer representing the maximum possible distance 
# (in base pairs) between any 2 connected anchors
#
#### Dependencies: see R libraries below
#
#### Outputs: 1) datatable of all possible combinations of (Rec-1) cis-EE/EP/PP interactions with columns 
# "chr_start_stop(pro/enh1)" "chr_start_stop(pro/enh2)" "start(pro/enh1)" "start(pro/enh2)" "length" "loop_id"
# 2) datatable containing all accessible anchors with annotations as enhancer or promoter with columns
# "chrom_enh/pro" "start_enh/pro" "end_enh/pro" "chrom_atac" "start_atac" "end_atac" "all_ann[ignore]" "annotation"
#### ---------------------------------------------------------------------------------------------//

#### Load dependencies
library("dplyr")
library("cluster")
library('ggplot2')
library('data.table')
# install.packages("stringr")
library('stringr')

#### Set working directory
mydir <- "test_HiC"
setwd(mydir)

#### Input 1: Table containing all hg19 genes annotated with expression status in relevant conditions (i.e. Ib-sens/Ib-treated/ Ib-res) from RNA-seq
# with columns: Geneid	symbol	DMSO_express(boolean)	IBR_express(boolean)	RES_express(boolean)
REC1_RNA <- read.table("Interaction_Filtering_HiC/example_input_data/REC1_ENSG_all_exp.txt", header=TRUE, sep = "\t")

#### Input 2: List of enhancer & promoter anchors intersected with ATAC peaks from "*intersect_anchor_peak.sh" script
enh_tss_all <- fread("230306_REC1_TSS_ENH_ATAC_intersect.bed")

#### Input 3: Integer representing the longest possible cis interaction distance, likely not longer than 2MB based on the longest observed
# input interactions
max_anchor_distance <- 2000000

#### Create union list of genes that are expressed in any of the relevant conditions of the cell model (e.g. Ib-sens / Ib-treated sens / res)
expressed <- REC1_RNA[REC1_RNA$DMSO_express == TRUE | REC1_RNA$IBR_express == TRUE | REC1_RNA$RES_express == TRUE, ]

#### Annotate all accessible enhancer and promoter anchors and filter promoters by expression status
# Annotate all expressed anchors as promoters and non-expressed anchors as enhancers
enh_tss_all[,`:=`(all_ann = ifelse(V4%in%expressed$Geneid, "promoter", "enhancer")),]
# Eliminate all non-expressed anchors that are not enhancers (i.e. eliminate non-expressed promoter anchors)
enh_tss <- enh_tss_all[!(grepl("ENSG", V4) & all_ann == "enhancer"),]
enh_tss[, `:=`(annotation = ifelse("promoter"%in%all_ann, "promoter", "enhancer")), by = .(V8)]

#### Make dataframe containing all enhancer & promoter anchors annotated by atac_id
all_EP_by_atac_df <- unique(subset(enh_tss, select=c("V5", "V8")))

#### Iterate through anchors by chromosome & find all (cis) EE/EP/PP combinations within max_anchor_distance
for (chrom in unique(all_EP_by_atac_df$V5)){

  # Create dataframe holding all combinations of anchors (excl. trans interactions)
  all_EP_by_atac_id <- unique(all_EP_by_atac_df[V5 == chrom, ])
  all_EP_by_atac_id1 <- all_EP_by_atac_id$V8
  all_EP_by_atac_id2 <- all_EP_by_atac_id1
  temp1 <- expand.grid(all_EP_by_atac_id1, all_EP_by_atac_id2)

  # Filter out self-self interactions, duplicates in reverse orientation (ex: p,q <--> q,p), 
  # and interactions between regulatory elements that are > max_anchor_distance apart
  temp2 <- temp1 %>%
    mutate(strtoi(str_split_i(Var1, "_", 2)),
           strtoi(str_split_i(Var2, "_", 2)))
  colnames(temp2) <- c("reg_atac_id1", "reg_atac_id2", "start1", "start2")
  
  temp2_filt <- temp2[(temp2$start1 < temp2$start2), ]
  
  temp2_filt <- temp2_filt %>%
    mutate(length = abs(start1 - start2))
  temp2_filt <- temp2_filt[temp2_filt$length < max_anchor_distance, ]
  
  # If first iteration of loop, create final output dataframe
  if (chrom == "chr1"){
    all_reg_int_cis_only <- temp2_filt
  }
  else{
    # If subsequent iteration, append dataframe to output dataframe:
    all_reg_int_cis_only <- rbind(all_reg_int_cis_only, temp2_filt)
  }
}

#### Annotate output dataframe with unique interaction_ids based on atac_ids of the 2 connecting anchors
all_reg_int_cis_only$loop_id <- paste(str_split_i(all_reg_int_cis_only$reg_atac_id1, "_", 1), str_split_i(all_reg_int_cis_only$reg_atac_id1, "_", 2), str_split_i(all_reg_int_cis_only$reg_atac_id2, "_", 3), sep = "_")

#### Outputs: 1) Dataframe of all possible accessible EE/EP/PP interactions and 2) dataframe of enhancer/promoter anchor annotations
# both for intersection with accessible HiC interactions .bedpe file in the subsequent "*quantify_bedpe_HiC.R" script
save(all_reg_int_cis_only, enh_tss, file = "230306_random_REC1_EP_table.rda")

#### [Optional] Annotate output df
# all_reg_int_expand <- merge(all_reg_int_cis_only, enh_tss, by.x = "reg_atac_id1", by.y = "V8", all.x = T)
