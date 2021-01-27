#!/bin/R

#########################
# Author: Molly Brothers
# Github: mollybrothers
# Date: 2021-01-26
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
plot_methylation <- function(data, chr) {
  ggplot(
    data,
    aes(x = start, y = percentage)
    ) +
    geom_point(position = "jitter", color = "mediumpurple4", alpha = 1/4) +
    ylim(0,100) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size=15)) +
    labs(title = "SIR3-EcoGII (JRY13027)",
         x = sprintf("position on chr %s", chr),
         y = "% methylated reads")
}

#plot a negative control region
control <- relevant[chrom == "X" & start > 50e3 & start < 75e3]
plot_methylation(control, "X")

#plot HML
HML <- relevant[chrom == "III" & start > 0 & start < 25e3]
plot_methylation(HML, "III")
dim(HML)

HML$window <- cut(HML$start, breaks = 150)
HMLmean <- aggregate(HML$percentage, by = list(HML$window), function(x){ mean(x[x != 0])})
HMLsd<- aggregate(HML$percentage, by = list(HML$window), function(x){ sd(x[x != 0])})

HMLfinal <- data.frame(interval = HMLmean[,1], mean = HMLmean[,2], upper = HMLmean[,2]+HMLsd[, 2], lower = HMLmean[,2]-HMLsd[, 2])
HMLfinal$interval <- gsub("[(]", "", HMLfinal$interval)
HMLfinal$interval <- gsub("]", "", HMLfinal$interval)
HMLfinal$start <- as.numeric(unlist(lapply(str_split(HMLfinal$interval, ","), "[", 1)))
HMLfull <- HML[HML$percentage == 100, ]

HML$window <- gsub("[(]", "", HML$window)
HML$window <- gsub("]", "", HML$window)
HML$window <- as.numeric(unlist(lapply(str_split(HML$window, ","), "[", 1)))


ggplot(HMLfinal, aes(x = start, y = mean)) +
  geom_point(color = "mediumpurple4") +
  geom_point(data = HMLfull, aes(x = start, y = percentage), color = "mediumpurple4", inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = lower, ymax = upper), color = "mediumpurple4", inherit.aes = TRUE) + 
  scale_x_continuous(breaks = seq(0, 25000, 2000)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, size = 10, vjust = 1, hjust = 1)) +
  labs(title = "SIR3-EcoGII (JRY13027)",
       x = "position on chr III",
       y = "% methylated reads")

#plot HMR
HMR <- relevant[chrom == "III" & start > 285e3 & start < 300e3]
plot_methylation(HMR, "III")

#plot telomere
tel6R <- relevant[chrom == "VI" & start > 250e3 & start < 270e3]
plot_methylation(tel6R, "VI")

#plot chromosome
chrXV <- relevant[chrom == "XV"]
plot_methylation(chrXV, "XV")
