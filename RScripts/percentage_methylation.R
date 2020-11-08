#!/bin/R

#######################
# Author: Molly Brothers
# Github: mollybrothers
# Date: 2020-11-08
#######################

# With an input bedMethyl file:
# (https://www.encodeproject.org/data-standards/wgbs/),
# this script will plot the percentage of methylation
# at each position as a scatterplot

library(data.table)
library(tidyverse)

dt <- fread("~/sequencing/201012_Doudna/modified_bases.6mA.bed")
colnames(dt) <- c("chrom", "start", "end", "name", "score", 
                  "strand", "startCodon", "stopCodon", "color", 
                  "coverage", "percentage")

#get a new data table only containing chrom, start, coverage, and percentage
#filter out mitochondrial DNA (MT), rDNA (XII), and coverage < 15

select_cols <- c("chrom", "start", "coverage", "percentage")
relevant <- dt[chrom != "MT" & chrom != "XII" & coverage > 15, ..select_cols]

#plot percentage of methylation in a particular region as a scatter plot
#opacity of dot corresponds to amount of coverage
plot_methylation <- function(data, chr) {
  ggplot(
    data,
    aes(x = start, y = percentage, alpha = coverage)
    ) +
      theme_minimal() +
      geom_point(position = "jitter", color = "mediumpurple4") +
      xlab(sprintf("position on %s", chr)) +
      ylab("% methylated reads")
}

#plot a negative control region
control <- relevant[chrom == "X" & start > 50e3 & start < 73e3]
plot_methylation(control, "X")

#plot HML
HML <- relevant[chrom == "III" & start > 0 & start < 25e3]
plot_methylation(HML, "III")

#plot HMR
HMR <- relevant[chrom == "III" & start > 285e3 & start < 300e3]
plot_methylation(HMR, "III")

