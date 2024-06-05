#### --------------------------------------------------------------------------------------------------------------------------//
#### Overview: This script can be used to annotate an input list of hubs from a single condition by categorizing individual hubs 
# within the list as hyperconnected (or regularly connected) as well as a graph of cliques by interaction count (optional).                                                                     
#
#### Inputs: Single list of hubs (ex: control_clique_by_interval.txt or case_clique_by_interval.txt or could be the *unionlist*.txt file but avoid if possible b/c less accurate),
# working directory, output file names
#
#### Dependencies: see libraries below. 
#
#### Final Outputs: 1) A .tsv file of hubs that is annotated for hyperconnected cliques in format: "hub_id /t hub_interaction_count /t	i /t hubClass", 
# 2) a .txt file of only hyperconnected hubs in the format "chrom /t start /t end /t hub_id /t interaction_count /n"
# 3) OPTIONAL-- graph of hub rank vs. hub interaction count with line denoting hyperconnected hub cutoff
#### -------------------------------------------------------------------------------------------------------------------------//

#### Load dependencies 
library("plyr")
library("dplyr")
library("cluster")
library('ggplot2')
library("reshape")
library("sets")
library("amap")
library("data.table")                                                                                                                   
library("stringr")
library("parallel")
library("tidyr")
library("extrafont")
# plotting
library('Sushi')

# graph
library("igraph")

#### Inputs
# Set working directory
setwd("/mnt/data0/brent/analysis/idea32_230824_final_hub_scripts_submission/230825_Finalized_Hub_Scripts/Hub_Pipeline/example_input_data")

# Input file to this script should be control_clique_by_interval.txt OR case_clique_by_interval.txt
# or can be the hub unionlist output of the differential hub program (AVOID IF POSSIBLE-- less accurate)
input_df <- read.table("Case_hub_by_interval_sort.txt", header=FALSE, sep="\t")

# Define output file names 
output_filename = "REC1_Arima_RES"
output2 = "REC1_Arima_RES_hypercliques.txt"
graph_out <- "REC1_Arima_RES_hub_edge_rank_new.pdf" #OPTIONAL

#******************************************************************************************************************
#this is an accessory function, that determines the number of points below a diagnoal passing through [x,yPt] [DON'T CHANGE]
#******************************************************************************************************************
numPts_below_line <- function(myVector,slope,x){
  yPt <- myVector[x]
  b <- yPt-(slope*x)
  xPts <- 1:length(myVector)
  return(sum(myVector<=(xPts*slope+b)))
}

#******************************************************************************************************************
#This function calculates the cutoff by sliding a diagonal line and finding where it is tangential (or as close as possible) [DON'T CHANGE]
#******************************************************************************************************************
calculate_cutoff <- function(inputVector,...){
  inputVector <- sort(inputVector)
  inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
  slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
  xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum)
  #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
  y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.

  return(list(absolute=y_cutoff,overMedian=y_cutoff/median(inputVector),overMean=y_cutoff/mean(inputVector)))
}

#******************************************************************************************************************
# Main function for super-hub-- adapted from GW/YZ/BF code [DON'T CHANGE]
#******************************************************************************************************************
# Find superhub // hyperconnected cliques
# input is a 3-column sorted community size: 1- c.number, 2-size, 3- Ascending rank of size
find_SHs_from_size <- function(size_file, output_filename){
  
  # Take only the size column 
  inputVector = size_file[,2]
  
  # Calculate tangent point of curve as edge number cutoff for hyperconnected clique 
  cutoff_options <- calculate_cutoff(inputVector)

  superHubRows <- which(inputVector > cutoff_options$absolute)
  OBS <- size_file %>% mutate(hubClass = "norm")
  OBS[superHubRows, "hubClass"] = "hyper"
  
  write.table(OBS, paste(output_filename,'_hyperconnected_cliques_annotated.tsv',sep=''),
              row.names=F, quote=F, sep='\t', col.names=T)
  return(OBS)
  
}


#******************************************************************************************************************
# Code to annotate superhubs
#******************************************************************************************************************

#### If necessary, rename cols of input_df (makes it easier to merge later):
colnames(input_df) <- c("chr", "start", "end", "Hub_id", "Interaction_count")
input_dfs = subset(input_df, select=c(Hub_id, Interaction_count))

#### Mutate input dataframe to calculate hyperconnected cliques:
# col1 is clique number identifier ("c.number")
# col2 is edge_number of clique ("size")
# col3 is post-sort row number ("i")

colnames(input_dfs) <- c("c.number", "size")
S.DF= arrange(input_dfs, size)

#### Remove cliques that have edge_number < 6 
S.DF <- S.DF[S.DF$size > 5, ]

#### Add rank variable 'i' 
S.DF = mutate(S.DF, i = 1:dim(S.DF)[1])

#### Call hyperconnected clique function
S.DF_superHub = find_SHs_from_size(S.DF, output_filename)

#### Merge hyperconnected clique output back with original input file so that 
# hyperconnected cliques can be extracted for downstream analysis 
S.DF_superHub <- S.DF_superHub[S.DF_superHub$hubClass == "hyper",]
out = subset(S.DF_superHub, select=c(c.number, size))

colnames(out) <- c("Hub_id", "Interaction_count")
out_final = merge(out, input_df, by = "Hub_id")
out_final = subset(out_final, select=c(chr, start, end, Hub_id, Interaction_count.x))
colnames(out_final) <- c("chr", "start", "end", "Hub_id", "Interaction_count")
write.table(out_final, output2, row.names=F, quote=F, sep='\t', col.names=T)

#******************************************************************************************************************
# [OPTIONAL-- simplified version of code used to generate graph "A" in figures 2-4] 
# CODE TO GRAPH INTERACTION COUNT VS. HUB RANK GRAPHS WITH HYPERHUB LINES (USING INDIVIDUAL CASE/CONTROL HUB FILES)
#******************************************************************************************************************

#### Graphing style
p0 <- theme_bw() + theme( plot.title = element_text(hjust = 0.5),
                          panel.background = element_rect(fill = "white", colour = NA),
                          panel.border = element_rect( fill = NA, colour = "black", size = 1),
                          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                          text= element_text(family="sans", size=22) )

#### Load ggplot2
library("ggplot2")

#### Open file that has all cliques in 1 condition annotated by hyperconnected status
filename = paste(output_filename,'_hyperconnected_cliques_annotated.tsv',sep='')
clique_list <- read.table(file = filename, sep = '\t', header = TRUE)

#### Find edge_num and rank cutoffs for hyperconnected cliques
idx = match("hyper", clique_list$hubClass)
hyper_ycutoff = clique_list$size[idx]
hyper_xcutoff = idx

#### Graph & save as PDF
pdf(file = graph_out)

ggplot(clique_list, aes(y=size, x=i)) + geom_point() + p0 +
  geom_hline(yintercept=hyper_ycutoff, lty=2) +
  geom_vline(xintercept=hyper_xcutoff, lty=2)

dev.off()
