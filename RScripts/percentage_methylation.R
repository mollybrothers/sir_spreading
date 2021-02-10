#!/bin/R

#########################
# Author: Molly Brothers
# Github: mollybrothers
# Date: 2021-02-10
#########################

# With an input bedMethyl file:
# (https://www.encodeproject.org/data-standards/wgbs/),
# this script will plot the percentage of methylation
# at each position as a scatterplot

library(data.table)
library(tidyverse)
library(stringr)

dt_37C <- fread("/Volumes/brothers_seq/Nanopore/210205_Ocular/megalodon_output_11/modified_bases.aggregate11.6mA.bed")
colnames(dt_37C) <- c("chrom", "start", "end", "name", "score", 
                  "strand", "startCodon", "stopCodon", "color", 
                  "coverage", "percentage")
dt_25C <- fread("/Volumes/brothers_seq/Nanopore/210205_Ocular/megalodon_output_08/modified_bases.aggregate08.6mA.bed")
colnames(dt_25C) <- c("chrom", "start", "end", "name", "score", 
                      "strand", "startCodon", "stopCodon", "color", 
                      "coverage", "percentage")

#get a new data.table only containing chrom, start, coverage, and percentage
#filter out mitochondrial DNA (MT) and coverage < 10

select_cols <- c("chrom", "start", "coverage", "percentage")
relevant_37C <- dt_37C[chrom != "MT" & coverage > 10, ..select_cols]
relevant_25C <- dt_25C[chrom != "MT" & coverage > 10, ..select_cols]

#plot percentage of methylation in a particular region as a scatter plot
#opacity of dot corresponds to amount of coverage
plot_title = "sir3-8-EcoGII (JRY13114)"

plot_methylation_dot <- function(data, chr) {
  ggplot(data, aes(x = start, y = percentage)) +
    geom_point(color = "mediumpurple4", alpha = 1/4) +
    scale_x_continuous(limits = c(min(data$start), max(data$start)), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0,100), expand = c(0, 0)) +
    theme(panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.line = element_line(color = "black"),
          text = element_text(size=15, color = "black", family = "Arial")) +
    labs(title = plot_title,
         x = sprintf("position on chr %s", chr),
         y = "% methylated reads")
}

plot_methylation_bars <- function(data, data2, chr) {
  #break data into 200bp bins, find the mean and sd of % methylated reads in each bin
  numbreaks <- as.integer((max(data$start) - min(data$start)) / 200)
  data$window <- cut(data$start, breaks = numbreaks)
  data_mean <- aggregate(data$percentage, by = list(data$window), function(x){ mean(x)})
  data_sd <- aggregate(data$percentage, by = list(data$window), function(x){ sd(x)})
  
  numbreaks2 <- as.integer((max(data2$start) - min(data2$start)) / 200)
  data2$window <- cut(data2$start, breaks = numbreaks2)
  data2_mean <- aggregate(data2$percentage, by = list(data2$window), function(x){ mean(x)})
  data2_sd <- aggregate(data2$percentage, by = list(data2$window), function(x){ sd(x)})
  
  #make a data.frame of the stats found above and clean up for plotting
  data_stats <- data.frame(interval = data_mean[,1], mean = data_mean[,2], upper = data_mean[,2]+data_sd[, 2], lower = data_mean[,2]-data_sd[, 2])
  data_stats$interval <- gsub("[(]", "", data_stats$interval)
  data_stats$interval <- gsub("[,].*", "", data_stats$interval)
  data_stats$interval <- as.numeric(data_stats$interval)
  
  data2_stats <- data.frame(interval = data2_mean[,1], mean = data2_mean[,2], upper = data2_mean[,2]+data2_sd[, 2], lower = data2_mean[,2]-data2_sd[, 2])
  data2_stats$interval <- gsub("[(]", "", data2_stats$interval)
  data2_stats$interval <- gsub("[,].*", "", data2_stats$interval)
  data2_stats$interval <- as.numeric(data2_stats$interval)
  
  #get the points that have 100% reads methylated
  data_full <- data[data$percentage == 100, ]
  data2_full <- data2[data2$percentage == 100, ]
  
  #plot the bins and the points that are at 100%
  ggplot(mapping = aes(x = interval, y = mean)) +
    geom_point(data = data_stats, color = "mediumpurple4") +
    geom_point(data = data2_stats, color = "cyan3") +
    geom_point(data = data_full, mapping = aes(x = start, y = percentage), color = "mediumpurple4", inherit.aes = FALSE) +
    geom_point(data = data2_full, mapping = aes(x = start, y = percentage), color = "cyan3", inherit.aes = FALSE) +
    geom_errorbar(data = data_stats, mapping = aes(ymin = lower, ymax = upper), color = "mediumpurple4", inherit.aes = TRUE) + 
    geom_errorbar(data = data2_stats, mapping = aes(ymin = lower, ymax = upper), color = "cyan3", inherit.aes = TRUE) +
    theme(panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.line = element_line(color = "black"),
          text = element_text(size = 15, color = "black", family = "Arial")) +
    labs(title = plot_title,
         x = sprintf("position on chr %s", chr),
         y = "average % methylated reads")
}

#plot a negative control region
control_25C <- relevant_25C[chrom == "III" & start > 80e3 & start < 105e3]
control_37C <- relevant_37C[chrom == "III" & start > 80e3 & start < 105e3]
plot_methylation_bars(control_25C, control_37C, "III")
plot_methylation_dot(control, "III")

#plot HML
HML_25C <- relevant_25C[chrom == "III" & start > 0 & start < 25e3]
HML_37C <- relevant_37C[chrom == "III" & start > 0 & start < 25e3]
plot_methylation_bars(HML_25C, HML_37C, "III")
plot_methylation_dot(HML, "III")

#plot HMR
HMR_25C <- relevant_25C[chrom == "III" & start > 280e3 & start < 305e3]
HMR_37C <- relevant_37C[chrom == "III" & start > 280e3 & start < 305e3]
plot_methylation_bars(HMR_25C, HMR_37C, "III")
plot_methylation_dot(HMR, "III")

#plot telomere
tel6R_25C <- relevant_25C[chrom == "VI" & start > 250e3 & start < 270e3]
tel6R_37C <- relevant_37C[chrom == "VI" & start > 250e3 & start < 270e3]
plot_methylation_bars(tel6R_25C, tel6R_37C, "VI")
plot_methylation_dot(tel6R, "VI")

tel9L_25C <- relevant_25C[chrom == "IX" & start > 0 & start < 25e3]
tel9L_37C <- relevant_37C[chrom == "IX" & start > 0 & start < 25e3]
plot_methylation_bars(tel9L_25C, tel9L_37C, "IX")

#plot chromosome
chrXV <- relevant[chrom == "XV"]
plot_methylation_dot(chrXV, "XV")
plot_methylation_bars(chrXV, "XV")
