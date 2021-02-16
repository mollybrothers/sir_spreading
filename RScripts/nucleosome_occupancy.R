#!/bin/R

#################################################
# Author: Molly Brothers and Koen Van den Berge
# Github: mollybrothers
# Date: 2021-02-11:
#################################################

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

#############################################################
####FIND AND PLOT AVERAGE METHYLATION IN NUCS AND LINKERS####
#############################################################
average_methylation_nucs <- function(nucs_data, meth_data){  
  LtoN <- which(diff(nucs_data$occupancy < 0.4) == -1)
  NtoL <- which(diff(nucs_data$occupancy < 0.4) == 1)
  nuc_rows <- data.table(start = c(1,LtoN), end = c(NtoL, nrow(nucs_data)))
  nuc_rows <- nuc_rows[(nuc_rows$end - nuc_rows$start) >=10,]
  
  avg_meth_nucl <- c()
  for (i in 1:nrow(nuc_rows)) {
    sub_nucs <- nucs_data[nuc_rows$start[i]:nuc_rows$end[i]]
    begin <- min(sub_nucs$start)
    finish <- max(sub_nucs$end)
    
    sub_meth <- meth_data[start > begin & start < finish]
    avg_meth_nucl[i] <- mean(sub_meth$percentage)
  }
  return(avg_meth_nucl)
}


average_methylation_linkers <- function(nucs_data, meth_data) {
  LtoN <- which(diff(nucs_data$occupancy < 0.4) == -1)
  NtoL <- which(diff(nucs_data$occupancy < 0.4) == 1)
  linker_rows <- data.table(start = NtoL, end = LtoN)
  linker_rows <- linker_rows[(linker_rows$end - linker_rows$start) >=10,]
  
  avg_meth_linker <- c()
  for (i in 1:nrow(linker_rows)) {
    sub_linkers <- nucs_data[linker_rows$start[i]:linker_rows$end[i]]
    begin <- min(sub_linkers$start)
    finish <- max(sub_linkers$end)
    
    sub_meth <- meth_data[start > begin & start < finish]
    avg_meth_linker[i] <- mean(sub_meth$percentage)
  }
  return(avg_meth_linker)
}

methyl_boxplot <- function (nucs_methyl, linkers_methyl) {  
  combined <- data.table(region = c(rep("nucleosomes",length(nucs_methyl)),
                                    rep("linkers",length(linkers_methyl))),
                         avg_methylation = c(nucs_methyl, linkers_methyl))
  
  ggplot(combined, aes(x = region, y = avg_methylation)) +
    geom_boxplot(fill = "lavender") +
    geom_point() +
    ylim(c(0,80)) +
    theme(
      panel.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.line = element_line(color = "black"),
      text = element_text(size = 15, color = "black", family = "Arial")) +
    labs(x = "region",
         y = "average % methylated reads")
  
}

############
####DATA####
############

#control
ctrl_meth <- filtered_meth[chrom == "III" & start > 96e3 & start < 99e3]
ctrl_nucs <- nucs[chrom == "chrIII" & start > 96e3 & start < 99e3]
ctrl_linkers_meth <- average_methylation_linkers(ctrl_nucs, ctrl_meth)
ctrl_nucs_meth <- average_methylation_nucs(ctrl_nucs, ctrl_meth)
#methyl_vs_nucs(ctrl_meth, ctrl_nucs)
methyl_boxplot(ctrl_nucs_meth, ctrl_linkers_meth)

#HMR
HMR_nucs <- nucs[chrom == "chrIII" & start > 291e3 & start < 294e3]
HMR_meth <- filtered_meth[chrom == "III" & start > 291e3 & start < 294e3]
HMR_nucs_meth <- average_methylation_nucs(HMR_nucs, HMR_meth)
HMR_linkers_meth <- average_methylation_linkers(HMR_nucs, HMR_meth)
#methyl_vs_nucs(HMR_meth, HMR_nucs)
methyl_boxplot(HMR_nucs_meth, HMR_linkers_meth)

#HML
HML_nucs <- nucs[chrom == "chrIII" & start > 11e3 & start < 14e3]
HML_meth <- filtered_meth[chrom == "III" & start > 11e3 & start < 14e3]
HML_nucs_meth <- average_methylation_nucs(HML_nucs, HML_meth)
HML_linkers_meth <- average_methylation_linkers(HML_nucs, HML_meth)
#methyl_vs_nucs(HML_meth, HML_nucs)
methyl_boxplot(HML_nucs_meth, HML_linkers_meth)

all_linkers <- c(HMR_linkers_meth, HML_linkers_meth)
all_nucs <- c(HMR_nucs_meth, HML_nucs_meth)
methyl_boxplot(all_nucs, all_linkers)

#tel6R
tel6R_nucs <- nucs[chrom == "chrVI" & start > 265e3 & start < 268e3]
tel6R_meth <- filtered_meth[chrom == "VI" & start > 265e3 & start < 268e3]
#methyl_vs_nucs(tel6R_meth, tel6R_nucs)
methyl_boxplot(tel6R_nucs, tel6R_meth)

#tel10L
tel10L_nucs <- nucs[chrom == "chrX" & start > 0 & start < 5e3]
tel10L_meth <- filtered_meth[chrom == "X" & start > 0 & start < 5e3]
#methyl_vs_nucs(tel10L_meth, tel10L_nucs)

#####################################
####SCATTER PLOT OF FITTED VALUES####
#####################################

# scatter plot of fitted values gives us a confused snake.
dfGrid <- data.frame(start = seq(291000, 294000), length.out=1000)
predMet <- predict(loMethyl, newdata=dfGrid)
predOccup <- predict(loOccup, newdata=dfGrid)
plot(x=predMet, y=predOccup)