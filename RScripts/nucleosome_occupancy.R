#!/bin/R

#########################
# Author: Molly Brothers
# Github: mollybrothers
# Date: 2021-01-28
#########################

# IN PROGRESS
# goal of this script is to compare directly in the same plot nucleosome occupancy data and % methylation
# to reveal linker-region preference of SIR3-EcoGII methylation

library(data.table)
library(tidyverse)

####################
####READ IN DATA####
####################

#percentage methylation data
# Molly's path
meth <- fread("/Volumes/brothers_seq/Nanopore/201012_Doudna/megalodon_output_00/modified_bases.6mA.bed")
# Savio path
# meth <- fread("/global/scratch/molly_brothers/201012_Doudna/megalodon_output_00/modified_bases.6mA.bed")
# Koen's path
# meth <- fread("~/data/molly/modified_bases.6mA.bed")
colnames(meth) <- c("chrom", "start", "end", "name", "score", 
                    "strand", "startCodon", "stopCodon", "color", 
                    "coverage", "percentage")
meth_cols <- c("chrom", "start", "coverage", "percentage")
filtered_meth <- meth[coverage > 10, ..meth_cols]

#nucleosome occupancy data from Chereji et al 2018
# Molly's path
nucs <- fread("/Volumes/brothers_seq/not_my_data/GSE97290_Henikoff_Chereji_2018/Chereji_Occupancy_rep1_singlebp.bedGraph")
# Savio path
# nucs <- fread("/global/scratch/molly_brothers/not_my_data/GSE97290_Henikoff_Chereji_2018/Chereji_Occupancy_rep1_singlebp.bedGraph")
# Koen's path
# nucs <- fread("~/data/molly/Chereji_Occupancy_rep1_singlebp.bedGraph")
colnames(nucs) <- c("chrom", "start", "end", "occupancy")

#######################################################
####LOESS PLOT METHYLATION VS. NUCLEOSOME OCCUPANCY####
#######################################################

methyl_vs_nucs <- function (methyl_data, nuc_data) {
  plot(x=methyl_data$start, y=rep(-3,3,length.out=nrow(methyl_data)), 
       type="n", ylim=c(-3,4), 
       ylab="Scaled average methylation / occupancy",
       xlab="Chromosome position")
  #legend("topleft", c("Methylation", "Nucleosome occupancy"), 
         #col=c("mediumpurple4", "grey50"), lty=1, cex=1,
         #bty='n', lwd=3)
  
  # Methylation
  loMethyl <- loess(percentage/100 ~ start, 
                    data = methyl_data, 
                    weights = methyl_data$coverage, 
                    enp.target = 50)
  lines(x=loMethyl$x[order(loMethyl$x)], y=scale(loMethyl$fitted[order(loMethyl$x)]), 
        col="mediumpurple4", lwd=3)
  
  # Nucleosome occupancy
  loOccup <- loess(occupancy/max(occupancy) ~ start, 
                   data=nuc_data, 
                   enp.target = 50)
  lines(x=loOccup$x, y=scale(loOccup$fitted), col="grey70", lwd=3)
}

#control
#why does the control region show methylation?
ctrl_meth <- filtered_meth[chrom == "III" & start > 96e3 & start < 99e3]
ctrl_nucs <- nucs[chrom == "chrIII" & start > 96e3 & start < 99e3]
methyl_vs_nucs(ctrl_meth, ctrl_nucs)

#HMR
HMR_meth <- filtered_meth[chrom == "III" & start > 291e3 & start < 294e3]
HMR_nucs <- nucs[chrom == "chrIII" & start > 291e3 & start < 294e3]
methyl_vs_nucs(HMR_meth, HMR_nucs)

#HML
HML_meth <- filtered_meth[chrom == "III" & start > 11e3 & start < 14e3]
HML_nucs <- nucs[chrom == "chrIII" & start > 11e3 & start < 14e3]
methyl_vs_nucs(HML_meth, HML_nucs)

#tel6R
tel6R_meth <- filtered_meth[chrom == "VI" & start > 265e3 & start < 268e3]
tel6R_nucs <- nucs[chrom == "chrVI" & start > 265e3 & start < 268e3]
methyl_vs_nucs(tel6R_meth, tel6R_nucs)

#tel8L
tel10L_meth <- filtered_meth[chrom == "X" & start > 0 & start < 5e3]
tel10L_nucs <- nucs[chrom == "chrX" & start > 0 & start < 5e3]
methyl_vs_nucs(tel10L_meth, tel10L_nucs)

#####################################
####SCATTER PLOT OF FITTED VALUES####
#####################################

# scatter plot of fitted values gives us a confused snake.
dfGrid <- data.frame(start = seq(291000, 294000), length.out=1000)
predMet <- predict(loMethyl, newdata=dfGrid)
predOccup <- predict(loOccup, newdata=dfGrid)
plot(x=predMet, y=predOccup)

####################################################
####PLOT AVERAGE METHYLATION IN NUCS AND LINKERS####
####################################################
NtoL <- which(diff(HMR_nucs$occupancy < 0.5) == 1)
LtoN <- which(diff(HMR_nucs$occupancy < 0.5) == -1)

nuc_rows <- data.table(start = c(1,LtoN), end = c(NtoL, nrow(HMR_nucs)))
linker_rows <- data.table(start = NtoL, end = LtoN)

nuc_rows <- nuc_rows[(nuc_rows$end - nuc_rows$start) >=10,]
linker_rows <- linker_rows[(linker_rows$end - linker_rows$start) >=10,]
  