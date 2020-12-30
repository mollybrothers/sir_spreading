#!/bin/R

#########################
# Author: Molly Brothers
# Github: mollybrothers
# Date: 2020-12-13
#########################

# IN PROGRESS
# goal of this script is to compare directly in the same plot nucleosome occupancy data and % methylation
# to reveal linker-region preference of SIR3-EcoGII methylation

library(data.table)
library(tidyverse)

#read in percentage methylation data
# Molly's path
# meth <- fread("/Volumes/brothers_seq/Nanopore/201012_Doudna/201012_Doudna_megalodon_output_00/modified_bases.6mA.bed")
# Savio path
# meth <- fread("/global/scratch/molly_brothers/201012_Doudna/megalodon_output_00/modified_bases.6mA.bed")
# Koen's path
meth <- fread("~/data/molly/modified_bases.6mA.bed")
colnames(meth) <- c("chrom", "start", "end", "name", "score", 
                    "strand", "startCodon", "stopCodon", "color", 
                    "coverage", "percentage")

meth_cols <- c("chrom", "start", "coverage", "percentage")
filtered_meth <- meth[coverage > 10, ..meth_cols]

#read in nucleosome occupancy data from Chereji et al 2018
# Molly's path
# nucs <- fread("/Volumes/brothers_seq/not_my_data/GSE97290_Henikoff_Chereji_2018/Chereji_Occupancy_rep1_singlebp.bedGraph")
# Savio path
# nucs <- fread("/global/scratch/molly_brothers/not_my_data/GSE97290_Henikoff_Chereji_2018/Chereji_Occupancy_rep1_singlebp.bedGraph")
# Koen's path
nucs <- fread("~/data/molly/Chereji_Occupancy_rep1_singlebp.bedGraph")
colnames(nucs) <- c("chrom", "start", "end", "occupancy")

#plot percentage methylation as a column chart
HMR <- filtered_meth[chrom == "III" & start > 291e3 & start < 294e3]
HMR_nucs <- nucs[chrom == "chrIII" & start > 291e3 & start < 294e3]

ggplot() +
  geom_col(HMR, mapping = aes(x = start, y = percentage), color = "mediumpurple4") +
  geom_col(HMR_nucs, mapping = aes(x = start, y = occupancy * 20), color = "black") +
  theme_minimal() +
  labs(title = "SIR3-EcoGII (JRY13027)",
       x = "position on chr III",
       y = "% methylated reads")


plot(x=HMR$start, y=rep(-3,3,length.out=nrow(HMR)), 
     type="n", ylim=c(-3,3), 
     ylab="Scaled average methylation / occupancy",
     xlab="Chromosome position")
# Methylation
loMethyl <- loess(percentage/100 ~ start, data=HMR, weights=HMR$coverage, enp.target = 50)
lines(x=loMethyl$x[order(loMethyl$x)], y=scale(loMethyl$fitted[order(loMethyl$x)]), 
      col="orange", lwd=3/2)
# Nucleosome occupancy
loOccup <- loess(occupancy/max(occupancy) ~ start, data=HMR_nucs, enp.target = 50)
lines(x=loOccup$x, y=scale(loOccup$fitted), col="darkseagreen3", lwd=3/2)
legend("topleft", c("Methylation", "Nucleosome occupancy"), 
       col=c("orange", "darkseagreen3"), lty=1, cex=1,
       bty='n', lwd=2)

# scatter plot of fitted values gives us a confused snake.
dfGrid <- data.frame(start = seq(291000, 294000), length.out=1000)
predMet <- predict(loMethyl, newdata=dfGrid)
predOccup <- predict(loOccup, newdata=dfGrid)
plot(x=predMet, y=predOccup)
