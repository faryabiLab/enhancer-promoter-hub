#### ---------------------------------------------------------------------------------------------//
#### This is the third of 4 scripts used to create a list of valid HiC EE/PP/EP interactions as input for hub calling. 
#### Overview: This script counts the number of observed (Rec-1 Ib-sensitive) HiC interactions connecting accessible enhancer
#  and (active) promoter anchors and normalizes raw counts by the total number of observed input HiC interations 
# to create hub calling input. This script should be run twice, once with condition 1 (i.e. drug-sens) Hi-C data as input and 
# again with condition 2 (i.e. drug-res) Hi-C data as input. 
# 
#### Inputs: 1) .bedpe file of all (Rec-1 Ib-sensitive) reads, 2) .bed file corresponding to 
# (REC-1) size-adjusted ATAC-seq peaks created as intermediate output by "1_intersect_anchor_peak.sh", 3) .rda file containing
# a) datatable with all possible (Rec-1) cis EE/EP/PP linkages (output from "2_make_random_EP_table.R") and b) datatable with 
# enhancer/promoter annotations for all anchors (see script 2 of 4 for more details on these 2 input files)
# 
#### Dependencies: see R libraries below
#
#### Output: Dataframe containing all accessible EE/EP/PP interactions in the given condition (i.e. Rec-1 Ib-sensitive) 
# as well as normalized scores for each interaction for subsequent filtering and hub-calling
#### ---------------------------------------------------------------------------------------------//

#### Load dependencies 
library('data.table')
library('parallel')
library('plyr')
library('dplyr')
library('stringr')

#### Set working directory
mainDir <- "test_HiC"
setwd(mainDir)

#### Input 1: Path to .bedpe file containing all (Ib-sensitive Rec-1) paired-end reads after HiC Pro alignment & basic filtering
loop <- "example_input_data/S01_220105_REC1_parental_DMSO_ArimaHiC_mango_ABRIDGED.rmdup.bedpe"
#loop <- "S03_220105_REC1_IBR_resistant_ArimaHiC_mango.rmdup.bedpe" #for ibrutinib-resistant loops

#### Input 2: Path to .bed file corresponding to (REC-1) size-adjusted ATAC-seq peaks created as intermediate output by "intersect_anchor_peak.sh"
atacPeak <- "test_HiC/REC1_ATAC_summit_ext_2500.txt"

#### Input 3: Load datatable of all possible cis linkages between accessible enhancers and promoters in Ib-sensitive Rec-1 and 
# datatable with enhancer/promoter annotations for all anchors from .rda output of "2_make_random_EP_table.R"
load("230306_random_REC1_EP_table.rda")

#### Execute helper bash shell script to create list of HiC interactions (from .bedpe input file) with two accessible anchors
# and annotate each of these anchors with unique ATAC-seq identifiers in two separate files (i.e. 'p' and 'q' anchors)
# can take ~30 mins depending on processing power
cmd <- paste("Interaction_Filtering_HiC/make_accessible_anchor_lists.sh", loop, atacPeak)
system(cmd)

#### Load helper script output files of anchors annotated by ATAC peaks, can take a few minutes [do not open files-- R will crash]
p <- fread("HiC_ATAC_2500_p.bed")
q <- fread("HiC_ATAC_2500_q.bed")

#### Create slimmed version of random EP table to expedite later quantification of observed HiC reads
# Extract atac id associated with each p and q anchor to create one dataframe of annotated EE/EP/PP interactions
p_slim <- subset(p, select=c("V4", "V8"))
colnames(p_slim) = c("p_atac_id", "linkage_id")
q_slim <- subset(q, select=c("V4", "V8"))
colnames(q_slim) = c("q_atac_id", "linkage_id")
pq_linked <- merge(p_slim, q_slim, by="linkage_id")

#### Re-order linkage such that p always < q to facilitate merge in next step
pq_linked_adjust <- pq_linked %>%
  mutate(start1 = strtoi(str_split_i(p_atac_id, "_", 2)),
         start2 = strtoi(str_split_i(q_atac_id, "_", 2)),
         switch = start1 > start2,
         p_atac_id2 = ifelse(switch, q_atac_id, p_atac_id), 
         q_atac_id2 = ifelse(switch, p_atac_id, q_atac_id))

pq_linked_adjust = subset(pq_linked_adjust, select=c("linkage_id", "p_atac_id2", "q_atac_id2"))

#### Prepare HiC pq linkage and random EP table for merge by creating 'pq_id' common identifier
pq_linked_adjust$pq_id = paste(pq_linked_adjust$p_atac_id2, pq_linked_adjust$q_atac_id2, sep = "_")
all_reg_int_cis_only$pq_id <- paste(all_reg_int_cis_only$reg_atac_id1, all_reg_int_cis_only$reg_atac_id2, sep = "_")

#### Merge HiC pq linkage and random EP table by common "pq_id"
# If 1 HiC linkage has anchor that maps onto multiple atac_ids, program automatically counts each overlapping linkage as separate read** 
all_reg_int_cis_only_filt <- merge(pq_linked_adjust, all_reg_int_cis_only, by="pq_id")
# Filter EP table to keep only interactions from HiC
all_reg_int_cis_only_filt <- unique(all_reg_int_cis_only_filt)
all_reg_int_cis_only_filt2 <- subset(all_reg_int_cis_only_filt, select=c("reg_atac_id1", "reg_atac_id2"))

# Count number of duplicated rows in dataframe (i.e. number of HiC reads occuring over a given random EP linkage) as raw interaction score
HiC_over_random_EP <- plyr::ddply(all_reg_int_cis_only_filt2, .(reg_atac_id1, reg_atac_id2), nrow)
colnames(HiC_over_random_EP) <- c("p_int_id", "q_int_id", "HiC_count")

#### Calculate normalized interaction score for each EPPPEE interaction by normalizing by number of detected interactions (i.e. rows) from original HiC .bedpe file 
# First, find number of loops detected from HiC
cmd1 <- paste("wc -l", loop)
loop_out <- system(cmd1, intern = TRUE)
total_loop <- strtoi(str_split_i(trimws(loop_out), " ", 1)) + 0 #total_loop <- 455394016 
# Second, normalize raw interaction counts
HiC_over_random_EP_2 <- HiC_over_random_EP %>% 
  mutate(HiC_over_random_EP_norm_2 = HiC_count/total_loop * 1E8)
#save(REC1_DMSO_Arima_HiC_over_random_EP, file="230306_REC1_DMSO_Arima_over_random_EP.rda") #optional 

#### Annotate scored interaction anchors by enhancers and promoters and create dataframe for clique-calling input
enh_tss_slim <- subset(enh_tss, select=c("V8", "annotation"))
HiC_over_random_EP_annot <- unique(merge(HiC_over_random_EP_2, enh_tss_slim, by.x="p_int_id", by.y="V8"))
colnames(HiC_over_random_EP_annot) = c("p_int_id", "q_int_id", "HiC_count", "HiC_over_random_EP_norm_2", "p_annotation")
HiC_over_random_EP_annot <- unique(merge(HiC_over_random_EP_annot, enh_tss_slim, by.x="q_int_id", by.y="V8"))
colnames(HiC_over_random_EP_annot) = c("q_int_id", "p_int_id", "HiC_count", "HiC_over_random_EP_norm_2", "p_annotation", "q_annotation")

#### Final output: dataframe containing all valid EE/PP/EP interactions, their anchor annotations, and their interaction scores
Control_HiC_clique_input <- HiC_over_random_EP_annot[,c(2,1,3,4,5,6)]
save(Control_HiC_clique_input, file="230306_REC1_DMSO_Arima_clique_input_random_EP.rda")

# Use this code instead when running this script using case/experimental data
# Case_HiC_clique_input <- HiC_over_random_EP_annot[,c(2,1,3,4,5,6)]
# save(Case_HiC_clique_input, file="230306_REC1_RES_Arima_clique_input_random_EP.rda")
