#### --------------------------------------------------------------------------------------------------------------------------//
#### This script is the second of four scripts used to create a list of valid HiChIP EE/PP/EP interactions as input for hub calling. 
#### Overview: This script annotates and filters a list of SMC1 HiChIP interactions to keep only interactions that are between valid,
# expressed promoter(s) and/or enhancer(s) for a given condition.                                                                                                                              
#
#### Inputs: R datatable of SMC1 HiChIP anchors intersected with size-adjusted, centered H3K27ac peaks and size-adjusted TSS
# for a given condition (output file from "1_intersect_anchor_peak_HiChIP.sh"), .bed file containing the intersection of all size-adjusted 
# H3K27ac ChIP-seq peaks and TSS intervals (output file from "1_intersect_anchor_peak_HiChIP.sh"), .rda file containing a 
# datatable with all expressed hg19 TSS for the given condition
#
#### Dependencies: See R libraries below.
#
#### Final Output:.rda file containing a dataframe of preliminary HiChIP PP/EE/EP interactions from a single condition 
#### -------------------------------------------------------------------------------------------------------------------------//

#### Clear R studio working environment and load in libraries/dependencies 
rm(list = ls())
library("plyr")
library("dplyr")

#### Inputs: set relevant directories for reading in and working with files
datadir <- "Interaction_Filtering_HiChIP/example_input_data"
wkdir <- "test_HiChIP"
setwd(wkdir)

#### Input: Create name of the final output .rda file from this script
out_rda_all <- "MB157_SMC1_DMSO_linkage.rda"

#### Inputs: Load files containing intersection of anchors, peaks and tss; 
# these are the final output files from the "1_intersect_anchor_peak_HiChIP.sh" script
mytab1 <- file.path(wkdir, "MB157_anchor_intersect_p_0.05Q.bed")
mytab2 <- file.path(wkdir, "MB157_anchor_intersect_q_0.05Q.bed")

#### Input: Load file containing intersection of H3K27ac ChIP-seq size-adjusted peaks and TSS to be used to determine
# annotation of anchors overlapping H3K27ac peaks and TSSes, this file is intermediate output from the "1_intersect_anchor_peak_HiChIP.sh" script
mytab3 <- file.path(wkdir, "MB157_peak_tss_intersect.bed")

#### Input: .rda file containing hg19 TSS and their expression status in the condition of interest
mytab4 <- paste(datadir, "MB157_ENSG_expressed.rda", sep = "/")

#### Input: read in list of expressed genes in (either) control (OR case condition) from RNA-seq data -----------------------------------//
load(mytab4)
expressed <- MB157_expressed

#### Read in intersection of significant FitHiChIP SMC1 HiChIP anchors and H3K27ac / TSS peaks
anchors.p <- read.table(mytab1, header = F, stringsAsFactors = F)
anchors.q <- read.table(mytab2, header = F, stringsAsFactors = F)

#### Filter out anchors that are connected to non-promoter/enhancer anchors and rename columns
# Column V4 = unique SMC1 HiChIP pairID 
anchors.p.keep <- anchors.p[anchors.p$V4%in%anchors.q$V4,]
anchors.q.keep <- anchors.q[anchors.q$V4%in%anchors.p$V4,]
colnames(anchors.p.keep) <- c("p_chr", "p_start", "p_end", "pair.id", "pval", "qval", "p.int.chr", "p.int.start", "p.int.end", "p.int.id")
colnames(anchors.q.keep) <- c("q_chr", "q_start", "q_end", "pair.id", "pval", "qval", "q.int.chr", "q.int.start", "q.int.end", "q.int.id")

#### Re-construct linkages of SMC1 HiChIP interactions using the filtered anchors 
anchor <- merge(anchors.p.keep, anchors.q.keep, by = c("pair.id", "pval", "qval"))

#### Label anchors that are promoters which are NOT expressed (or not labeled with proper enhancer/promoter IDs) with "NE"
anchor.k1 <- mutate(anchor, p.int.id = ifelse( (grepl("ENSG", p.int.id) & p.int.id%in%expressed$id) | grepl("MACS", p.int.id) , p.int.id, "NE"),
                    q.int.id = ifelse( (grepl("ENSG", q.int.id) & q.int.id%in%expressed$id) | grepl("MACS", q.int.id) , q.int.id, "NE"))
# Remove all interactions between nonexpressed promoter or invalid anchor(s)
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
DMSO.anchor.keep <- subset(anchor.keep, p.int.id != q.int.id & 
                         p.int.start < q.int.start & p.int.end < q.int.start & p_end < q_start)

save(list=c("DMSO.anchor.keep"), file = out_rda_all)
