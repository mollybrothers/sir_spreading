library(data.table)
library(tidyverse)

dt <- fread("~/sequencing/201012_Doudna/modified_bases.6mA.bed")
colnames(dt) <- c("chrom", "start", "end", "name", "score", 
                  "strand", "startCodon", "stopCodon", "color", 
                  "coverage", "percentage")

#get a new data frame only containing chrom, start, end, coverage, and percentage
#filter out mitochondrial DNA (MT), rDNA (XII), and coverage < 15

select_cols <- c("chrom", "start", "end", "coverage", "percentage")
relevant <- dt[chrom != "MT" & chrom != "XII" & coverage > 15, ..select_cols]

#plot percentage of methylation in a particular region as a scatter plot
#size of dot corresponds to amount of coverage
control <- relevant[chrom == "X" & start > 50e3 & start < 73e3]
ggplot(
  control,
  aes(x = start, y = percentage, alpha = coverage)
  ) +
  theme_minimal() +
  geom_point(position = "jitter", color = "mediumpurple4")

HML <- relevant[chrom == "III" & start > 0 & start < 25e3]
ggplot(
  HML, 
  aes(x = start, y = percentage, alpha = coverage)
  ) +
  theme_minimal() +
  geom_point(position = "jitter", color = "mediumpurple4")

HMR <- relevant[chrom == "III" & start > 285e3 & start < 300e3]
ggplot(
  HMR,
  aes(x = start, y = percentage, alpha = coverage)
  ) +
  theme_minimal() +
  geom_point(position = "jitter", color = "mediumpurple4")

#bin into 10bp bins, take average % methylation of those 10bp
