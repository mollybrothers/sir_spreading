#!/bin/R

#########################
# Author: Molly Brothers
# Github: mollybrothers
# Date: 2021-01-27
#########################

# With an input bedMethyl file:
# (https://www.encodeproject.org/data-standards/wgbs/),
# this script will plot the percentage of methylation
# at each position as a scatterplot

library(data.table)
library(tidyverse)
library(stringr)

dt <- fread("/Volumes/brothers_seq/Nanopore/201012_Doudna/megalodon_output_00/modified_bases.6mA.bed")
colnames(dt) <- c("chrom", "start", "end", "name", "score", 
                  "strand", "startCodon", "stopCodon", "color", 
                  "coverage", "percentage")

#get a new data.table only containing chrom, start, coverage, and percentage
#filter out mitochondrial DNA (MT) and coverage < 10

select_cols <- c("chrom", "start", "coverage", "percentage")
relevant <- dt[chrom != "MT" & coverage > 10, ..select_cols]

#plot percentage of methylation in a particular region as a scatter plot
#opacity of dot corresponds to amount of coverage
plot_title = "SIR3-EcoGII (JRY13027)"

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

plot_methylation_bars <- function(data, chr) {
  #break data into 200bp bins, find the mean and sd of % methylated reads in each bin
  numbreaks <- as.integer((max(data$start) - min(data$start)) / 200)
  data$window <- cut(data$start, breaks = numbreaks)
  data_mean <- aggregate(data$percentage, by = list(data$window), function(x){ mean(x)})
  data_sd <- aggregate(data$percentage, by = list(data$window), function(x){ sd(x)})
  
  #make a data.frame of the stats found above and clean up for plotting
  data_stats <- data.frame(interval = data_mean[,1], mean = data_mean[,2], upper = data_mean[,2]+data_sd[, 2], lower = data_mean[,2]-data_sd[, 2])
  data_stats$interval <- gsub("[(]", "", data_stats$interval)
  data_stats$interval <- gsub("[,].*", "", data_stats$interval)
  data_stats$interval <- as.numeric(data_stats$interval)
  
  #get the points that have 100% reads methylated
  data_full <- data[data$percentage == 100, ]
  
  #plot the bins and the points that are at 100%
  ggplot(data_stats, aes(x = interval, y = mean)) +
    geom_point(color = "mediumpurple4") +
    geom_point(data = data_full, aes(x = start, y = percentage), color = "mediumpurple4", inherit.aes = FALSE) +
    geom_errorbar(aes(ymin = lower, ymax = upper), color = "mediumpurple4", inherit.aes = TRUE) + 
    theme(panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.line = element_line(color = "black"),
          text = element_text(size = 15, color = "black", family = "Arial")) +
    labs(title = plot_title,
         x = sprintf("position on chr %s", chr),
         y = "% methylated reads")
}

#plot a negative control region
control <- relevant[chrom == "III" & start > 80e3 & start < 105e3]
plot_methylation_dot(control, "III")
plot_methylation_bars(control, "III")

#plot HML
HML <- relevant[chrom == "III" & start > 0 & start < 25e3]
plot_methylation_dot(HML, "III")
plot_methylation_bars(HML, "III")

#plot HMR
HMR <- relevant[chrom == "III" & start > 280e3 & start < 305e3]
plot_methylation_dot(HMR, "III")
plot_methylation_bars(HMR, "III")

#plot telomere
tel6R <- relevant[chrom == "VI" & start > 250e3 & start < 270e3]
plot_methylation_dot(tel6R, "VI")
plot_methylation_bars(tel6R, "VI")

#plot chromosome
chrXV <- relevant[chrom == "XV"]
plot_methylation_dot(chrXV, "XV")
plot_methylation_bars(chrXV, "XV")
