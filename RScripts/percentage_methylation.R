#!/bin/R

#########################
# Author: Molly Brothers
# Github: mollybrothers
# Date: 2020-12-01
#########################

# With an input bedMethyl file:
# (https://www.encodeproject.org/data-standards/wgbs/),
# this script will plot the percentage of methylation
# at each position as a scatterplot

library(data.table)
library(tidyverse)

dt <- fread("/Volumes/brothers_seq/Nanopore/201125_Turkey/megalodon_output_05/modified_bases.aggregate05.6mA.bed")
colnames(dt) <- c("chrom", "start", "end", "name", "score", 
                  "strand", "startCodon", "stopCodon", "color", 
                  "coverage", "percentage")

#get a new data table only containing chrom, start, coverage, and percentage
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
    theme_minimal() +
    labs(title = "unfused EcoGII (JRY12838)",
         x = sprintf("position on chr %s", chr),
         y = "% methylated reads")
}

#plot a negative control region
control <- relevant[chrom == "X" & start > 50e3 & start < 75e3]
plot_methylation(control, "X")

#plot HML
HML <- relevant[chrom == "III" & start > 0 & start < 25e3]
plot_methylation(HML, "III")

#plot HMR
HMR <- relevant[chrom == "III" & start > 285e3 & start < 300e3]
plot_methylation(HMR, "III")

#plot telomere
tel6R <- relevant[chrom == "VI" & start > 250e3 & start < 270e3]
plot_methylation(tel6R, "VI")

#plot chromosome
chrXV <- relevant[chrom == "XV"]
plot_methylation(chrXV, "XV")
