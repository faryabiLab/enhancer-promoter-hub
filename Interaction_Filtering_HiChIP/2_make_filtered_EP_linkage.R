#### --------------------------------------------------------------------------------------------------------------------------//
#### This script is the second of four scripts used to create a list of valid HiChIP EE/PP/EP interactions as input for hub calling. 
#### Overview: This script annotates and filters a list of SMC1 HiChIP loops to keep only interactions that are between valid,
# expressed promoter(s) and/or enhancer(s) for a given condition. The interactions in this dataframe will be assigned normalized
# contact frequency scores in the subsequent script to faciltate interaction filtering.                                                                                                                                
#
#### Inputs: 2 R datatables of SMC1 HiChIP anchors intersected with size-adjusted, centered H3K27ac peaks and size-adjusted TSS
# for a given condition (output files from "intersect_anchor_peak_HiChIP.sh"),.bedpe file of Mango pipeline output for 
# significant SMC1 HiChIP loops (PETS >= 4), .bed file containing the intersection of all size-adjusted H3K27ac ChIP-seq peaks 
# and TSS intervals, .rda file containing a datatable with all hg19 TSS and their observed, normalized RNA-seq signal in the 
# condition(s) of interest (including technical replicates for RNA-seq) 
#
#### Dependencies: See R libraries below.
#
#### Final Output:.rda file containing a dataframe of all valid SMC1 HiChIP PP/EE/EP loops from a single condition that will 
# be 'scored' in the subsequent script in the pipeline. 
#### -------------------------------------------------------------------------------------------------------------------------//

#### Clear R studio working environment and load in libraries/dependencies (many are optional)
rm(list = ls())
library("plyr")
library("dplyr")

#### Inputs: set relevant directories for reading in and working with files
mydir <- "/mnt/nas3/users/brent/analysis/Finalized_Hub_Scripts/Interaction_Filtering_HiChIP/example_input_data"
setwd(mydir)

#### Input: Create name of the final output .rda file from this script
out_rda_all <- "180408_REC1_SMC1_DMSO_linkage.rda"

#### Inputs: Load files containing intersection of anchors, peaks and tss; 
# these are the final output files from the "1_intersect_anchor_peak_HiChIP.sh" script
mytab1 <- file.path(mydir, "180408_REC1_anchor_intersect_p.bed")
mytab2 <- file.path(mydir, "180408_REC1_anchor_intersect_q.bed")

#### Input: Load file containing intersection of H3K27ac ChIP-seq size-adjusted peaks and TSS to be used to determine
# annotation of anchors overlapping H3K27ac peaks and TSSes, this file is intermediate output from the "1_intersect_anchor_peak_HiChIP.sh" script
mytab3 <- file.path(mydir, "180408_REC1_peak_tss_intersect.bed")

#### Input: File containing list of all SMC1 HiChIP reads occuring between validated, significant anchors (PETs >= 4),
# this is output from Mango HiChIP-processsing pipeline
mytab4 <- paste(mydir, "s99_180201_Rec1_DMSO_Smc1_FDR_0.05_PETS_4.mango", sep = "/")

#### Input: .rda file containing 'data.all.sorted' datatable in BED4 format containing all hg19 TSS and their 
# observed, normalized RNA-seq signal in the given condition(s) to be analyzed
mytab5 <- "REC1_WO_GSI_RNAseq_RPKM.rda"

#### Create list of expressed genes in (either) control (OR case condition) from RNA-seq data -----------------------------------//
# Read in R dataframe containing normalized (RPKM), processed RNA-seq data for condition(s) or interest
load(mytab5)
rownames(REC1.data.all.sorted) = REC1.data.all.sorted$id

#### Manually define "expressed" genes and create list of expressed genes
# In this example, an 'expressed' gene must have > 1 RPKM across at least 4 out of 6 RNA-seq replicates 
rpkm.col <- rowSums(REC1.data.all.sorted[, grep("rpkm", colnames(REC1.data.all.sorted))] > 1)
expressed <- REC1.data.all.sorted[rpkm.col > 3, ]

#### Read in intersection of valid Mango-filtered SMC1 HiChIP anchors and H3K27ac / TSS peaks
anchors.p <- read.table(mytab1, header = F, stringsAsFactors = F)
anchors.q <- read.table(mytab2, header = F, stringsAsFactors = F)

#### Filter out anchors that are connected to non-promoter/enhancer anchors and rename columns
# Column V4 = unique SMC1 HiChIP pairID 
anchors.p.keep <- anchors.p[anchors.p$V4%in%anchors.q$V4,]
anchors.q.keep <- anchors.q[anchors.q$V4%in%anchors.p$V4,]
colnames(anchors.p.keep) <- c("p_chr", "p_start", "p_end", "pair.id", "score", "p.int.chr", "p.int.start", "p.int.end", "p.int.id")
colnames(anchors.q.keep) <- c("q_chr", "q_start", "q_end", "pair.id", "score", "q.int.chr", "q.int.start", "q.int.end", "q.int.id")

#### Re-construct linkages of SMC1 HiChIP loops using the filtered anchors 
anchor <- merge(anchors.p.keep, anchors.q.keep, by = c("pair.id", "score"))

#### Label anchors that are promoters which are NOT expressed (or not labeled with proper enhancer/promoter IDs) with "NE"
anchor.k1 <- mutate(anchor, p.int.id = ifelse( (grepl("ENSG", p.int.id) & p.int.id%in%expressed$id) | grepl("MACS", p.int.id) , p.int.id, "NE"),
                      q.int.id = ifelse( (grepl("ENSG", q.int.id) & q.int.id%in%expressed$id) | grepl("MACS", q.int.id) , q.int.id, "NE"))
# Remove all loops between nonexpressed promoter or invalid anchor(s)
anchor.k2 <- subset(anchor.k1, p.int.id!= "NE" & q.int.id != "NE")

#### For anchors that overlap with both H3K27ac peaks and TSS, label them as promoters if they are expressed, else label them as enhancers
# Read in table containing intersection of H3K27ac size-adjusted, centered peaks and size-adjusted tss; only keep expressed anchors
intersect <- read.table(mytab3, header = F, stringsAsFactors = F)
intersect.keep <- intersect[intersect$V8%in%expressed$id, ]
# If anchor is contained within list of expressed anchors overlapping H3K27ac peaks, label it as a promoter, else label it an enhancer
anchor.k3 <- mutate(anchor.k2, p.int.ann = ifelse( p.int.id%in%intersect.keep$V4, "promoter", "enhancer"),
                      q.int.ann = ifelse( q.int.id%in%intersect.keep$V4, "promoter", "enhancer"))

#### Annotate remaining anchors as 'promoter' if it overlaps with expressed TSS or 'enhancer' if it overlaps with H3K27ac peak,
# as determined above
anchor.keep <- mutate(anchor.k3, p.int.ann = ifelse( grepl("ENSG", p.int.id), "promoter", p.int.ann),
                      q.int.ann = ifelse( grepl("ENSG", q.int.id), "promoter", q.int.ann))

#### Annotate linkages as between 2 promoters (PP), 2 enhancers (EE), or a promoter and an enhancer (PE)
anchor.keep <- mutate(anchor.keep, linkage.type = ifelse(p.int.ann == "promoter" & q.int.ann == "promoter", "PP", "PE" ))
anchor.keep <- mutate(anchor.keep, linkage.type = ifelse(p.int.ann == "enhancer" & q.int.ann == "enhancer", "EE", linkage.type ))

#### Filter out any linkages containing two of the same anchors or two anchors which overlap with one another
anchor.keep2 <- subset(anchor.keep, p.int.id != q.int.id & 
                        p.int.start < q.int.start & p.int.end < q.int.start & p_end < q_start)

#### Final output: save final list of SMC1 HiChIP loops between valid, annotated PP/EE/EP anchors as R dataframe in .rda file
# IF CONDITION 1, uncomment and run this code: 
DMSO.anchor.keep <- anchor.keep2 
save(list=c("DMSO.anchor.keep"), file = out_rda_all)

# (optional) IF CONDITION 2, uncomment and run this code: 
#GSI.anchor.keep <- anchor.keep2 
#save(list=c("GSI.anchor.keep"), file = out_rda_all)
