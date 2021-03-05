#!/bin/R

#########################
# Author: Molly Brothers
# Github: mollybrothers
# Date: 2021-03-04
#########################

#######################################################################################
# the per_read_base_calls.txt file itself is too large for RStudio's memory,
# so you'll need to use mawk on the command line to pick out the lines you want first
# EXAMPLE:
# cat per_read_modified_base_calls.txt | mawk '$2 ~ /^III$/ {print $0}' > chrIII.txt
#######################################################################################

library(data.table)
library(tidyverse)
library(wesanderson)

mega_directory <- "/Volumes/brothers_seq/Nanopore/201125_Turkey/megalodon_output_07/"
chr <- "III" #which chromosome?


probs <- fread(sprintf("%schr%s.txt", mega_directory, chr),
                    select = c(1, 2, 3, 4, 5),
                    col.names = c("read_id", "chrm", "strand", "pos", "mod_log_prob"))

#1. create a binary column; if 10^mod_log_prob is >0.8, set as methylated (m6A). if < 0.8, set as unmethylated (A)
#2. add start and end positions for each read_id
#3. order by position
probs <- probs[, methylation := ifelse(10^mod_log_prob > 0.8, TRUE, FALSE)][
  , list(read_id, pos, methylation, strand)][
    , start_pos := min(pos), by = read_id][
      , end_pos := max(pos), by = read_id][
        order(pos)]

#extract each unique read_id to set the order (in same order as in the data table to start with) for single read plots
read_names <- unique(probs$read_id)
probs$read_id = factor(probs$read_id, levels = read_names)

##################################
####PLOT AVERAGE PROBABILITIES####
##################################

plot_prob <- function(input, title) {
  ggplot(input, aes(x = pos, y = avg_prob)) +
    theme_classic() +
    geom_point(position = "jitter") +
    labs(
      title = sprintf("%s region", title),
      x = sprintf("chr %s", chr),
      y = "average log_mod_prob") +
    ylim(0, 1)
}

control_avg <- probs[pos %between% c(100e3, 105e3),
                     .(avg_prob = mean(mod_prob)),
                     by = pos
]
plot_prob(control_avg, "control")

HML_avg <- probs[pos %between% c(11e3, 16e3),
                 .(avg_prob = mean(mod_prob)),
                 by = pos
]
plot_prob(HML_avg, "HML")

HMR_avg <- probs[pos %between% c(291e3, 296e3),
                 .(avg_prob = mean(mod_prob)),
                 by = pos
]
plot_prob(HMR_avg, "HMR")

#########################
####PLOT SINGLE READS####
#########################

#plot function
plot_binary <- function(data) {
  ggplot(data, aes(x = pos, y = read_id, color = methylation)) +
    geom_point(shape = 15) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 15, color = "black", family = "Arial"),
          legend.position = "top",
          legend.title = element_blank(),
          legend.key = element_blank()) +
    scale_color_manual(values = c("grey90", "mediumpurple4"))
}

#plot unmethylated region
control_region <- c(95e3, 100e3)
control <- probs[start_pos <= control_region[1]][end_pos >= control_region[2]][
  pos %between% control_region]
plot_binary(control)

#plot HML
HML_region <- c(10.5e3, 15.5e3)
HML_silencers <- c(11146, 14849)
HML_linkers = c(#9407, 9587.5, 9067, 9747, 9923, 10166, 10331, 
  10585, 10748,
  10944, 11118, 11418, 11645, 12021, 12251, 12436, 12649, 12842,
  13017, 13396, 13558, 13829, 14011, 14221, 14883, 15229, 15406
  #15573, 15984, 16244
)
HML_segments <- readRDS("~/sequencing/sir_spreading/data/segmentsHML.rds")

HML <- probs[start_pos <= HML_region[1]][end_pos >= HML_region[2]][
  pos %between% HML_region]
hmlp <- plot_binary(HML)
hmlp
hmlp + geom_vline(xintercept = c(HML_segments[,1]), color = "red")
hmlp + geom_vline(xintercept = HML_silencers) #+
#geom_vline(xintercept = HML_linkers, size = 0.3, color = "black")

#plot HMR
HMR_region <- c(291e3, 296e3)
HMR_silencers = c(292388, 295034)
HMR_linkers = c(291100, 291312, 291644, 291863, 292129, 292322,
                292498, 292921, 293078, 293227, 293440, 293633, 293841, 294155,
                294515, 294699, 295239, 295555, 295743, 295906)
HMR_segments <- readRDS("~/sequencing/sir_spreading/data/segmentsHMR.rds")

HMR_all <- probs[start_pos <= HMR_region[1]][end_pos >= HMR_region[2]][
  pos %between% HMR_region]
HMR_plus <- probs[strand == "+"][start_pos <= HMR_region[1]][end_pos >= HMR_region[2]][
  pos %between% HMR_region]
HMR_minus <- probs[strand == "-"][start_pos <= HMR_region[1]][end_pos >= HMR_region[2]][
  pos %between% HMR_region]
hmrp <- plot_binary(HMR_all)
hmrp + geom_vline(xintercept = c(HMR_segments[,1]), color = "red")
hmrp + geom_vline(xintercept = HMR_silencers)

#plot tel3L
tel3L_region <- c(300, 5e3)
tel3L <- probs[start_pos <= tel3L_region[1]][end_pos >= tel3L_region[2]][
  pos %between% tel3L_region]
plot_binary(tel3L)

#plot tel6R
tel6R_region <- c(265e3, 270e3)
tel6R <- probs[start_pos <= tel6R_region[1]][end_pos >= tel6R_region[2]][pos %between% tel6R_region]
plot_binary(tel6R)

###############################
####HIERARCHICAL CLUSTERING####
###############################

#convert read_ids and positions to factors
HMR_plus$read_id <- factor(HMR_plus$read_id)
HMR_plus$pos <- factor(HMR_plus$pos)
HMRplus_readnames <- unique(HMR_plus$read_id)
HMRplus_positions <- unique(HMR_plus$pos)

#create a new dataframe with positions as column names and read_ids as row names
new_HMRplus <- data.frame()
for(i in HMRplus_positions){
  new_HMRplus[,i] <- NA
}
for(i in HMRplus_readnames){
  new_HMRplus[i,] <- logical()
}

#move the methylation data from the original dataframe to the new dataframe
for(x in HMRplus_readnames){
  for(y in HMRplus_positions){
    if(any(HMR_plus$read_id == x & HMR_plus$pos == y) == FALSE){
      next
    }
    new_HMRplus[as.character(x),as.character(y)] <- HMR_plus[which(HMR_plus$read_id == x & HMR_plus$pos == y),]$methylation
  }
}

#calculate distances and plot the resulting dendrogram
HMRplus_distances <- dist(new_HMRplus, method = "binary")
plot(hclust(HMRplus_distances, method = "average"))

# TO DO: figure out how to apply the clustering/distances to levels to plot the single reads by similarity/difference
