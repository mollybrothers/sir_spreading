#!/bin/R

###############################################
# Author: Molly Brothers / Koen Van den Berge
# Github: mollybrothers
# Date: 2021-04-24
###############################################

# Takes in aggregate bedMethyl file and plots loess-smoothed curve for different samples / time points

library(data.table)

columns <- c("chrom", "start", "end", "name", "score", 
             "strand", "startCodon", "stopCodon", "color", 
             "coverage", "percentage")
select_cols <- c("chrom", "start", "coverage", "percentage")
my_pal <- c("gray50", "forestgreen", "darkturquoise", "mediumpurple3", "deeppink", "black")

# change for different samples
methyl_0 <- fread("/home/mbrothers/nanopore/210531_Texas/modified_bases.aggregate01.6mA.bed",
                col.names = columns)
methyl_filtered_0 <- methyl_0[coverage > 10, ..select_cols]

methyl_30 <- fread("/home/mbrothers/nanopore/210531_Texas/modified_bases.aggregate02.6mA.bed",
                  col.names = columns)
methyl_filtered_30 <- methyl_30[coverage > 10, ..select_cols]

methyl_60 <- fread("/home/mbrothers/nanopore/210531_Texas/modified_bases.aggregate03.6mA.bed",
                  col.names = columns)
methyl_filtered_60 <- methyl_60[coverage > 10, ..select_cols]

methyl_90 <- fread("/home/mbrothers/nanopore/210531_Texas/modified_bases.aggregate04.6mA.bed",
                  col.names = columns)
methyl_filtered_90 <- methyl_90[coverage > 10, ..select_cols]

methyl_120 <- fread("/home/mbrothers/nanopore/210531_Texas/modified_bases.aggregate05.6mA.bed",
                   col.names = columns)
methyl_filtered_120 <- methyl_120[coverage > 10, ..select_cols]

methyl_150 <- fread("/home/mbrothers/nanopore/210531_Texas/modified_bases.aggregate06.6mA.bed",
                    col.names = columns)
methyl_filtered_150 <- methyl_150[coverage > 10, ..select_cols]

# uncomment for region of interest

# for HMR
chromo <- "III"
beg <- 292674 - 500 #285e3 #
end <- 294864 + 500 #305e3 #

# for HML
# chromo <- "III"
# beg <- 11237 - 500 #0  #
# end <- 14711 + 500 #20e3 #

# left telomeres
# chromo <- "XV"
# beg <- 0
# end <- 15000

# right telomeres
# chromo <- "VI"
# end <- 270161
# beg <- end - 15000

# other
# chromo <- "I"
# beg <- 124000
# end <- 125000


# filter to region of interest
region_methyl_0 <- methyl_filtered_0[chrom == chromo & start > beg & start < end]
region_methyl_30 <- methyl_filtered_30[chrom == chromo & start > beg & start < end]
region_methyl_60 <- methyl_filtered_60[chrom == chromo & start > beg & start < end]
region_methyl_90 <- methyl_filtered_90[chrom == chromo & start > beg & start < end]
region_methyl_120 <- methyl_filtered_120[chrom == chromo & start > beg & start < end]
region_methyl_150 <- methyl_filtered_150[chrom == chromo & start > beg & start < end]

#loess
lo_methyl_0 <- loess(percentage ~ start, data=region_methyl_0, weights=region_methyl_0$coverage, enp.target = 100)
lo_methyl_30 <- loess(percentage ~ start, data=region_methyl_30, weights=region_methyl_30$coverage, enp.target = 100)
lo_methyl_60 <- loess(percentage ~ start, data=region_methyl_60, weights=region_methyl_60$coverage, enp.target = 100)
lo_methyl_90 <- loess(percentage ~ start, data=region_methyl_90, weights=region_methyl_90$coverage, enp.target = 100)
lo_methyl_120 <- loess(percentage ~ start, data=region_methyl_120, weights=region_methyl_120$coverage, enp.target = 100)
lo_methyl_150 <- loess(percentage ~ start, data=region_methyl_150, weights=region_methyl_150$coverage, enp.target = 100)

#plot loess
plot(x=lo_methyl_0$x[order(lo_methyl_0$x)], y=lo_methyl_0$fitted[order(lo_methyl_0$x)], 
     type = 'l', lwd=2, ylim = c(0,100), col = my_pal[1], frame.plot = FALSE)
legend("topleft", legend = c("0min", "30min", "60min", "90min", "120min", "150min"), fill = my_pal)
lines(x=lo_methyl_30$x[order(lo_methyl_30$x)], y=lo_methyl_30$fitted[order(lo_methyl_30$x)], 
      type = 'l', lwd=2, ylim = c(0,100), col = my_pal[2])
# abline(v = c(68759, 68861))
lines(x=lo_methyl_60$x[order(lo_methyl_60$x)], y=lo_methyl_60$fitted[order(lo_methyl_60$x)], 
      type = 'l', lwd=2, ylim = c(0,100), col = my_pal[3])
lines(x=lo_methyl_90$x[order(lo_methyl_90$x)], y=lo_methyl_90$fitted[order(lo_methyl_90$x)], 
      type = 'l', lwd=2, ylim = c(0,100), col = my_pal[4])
lines(x=lo_methyl_120$x[order(lo_methyl_120$x)], y=lo_methyl_120$fitted[order(lo_methyl_120$x)], 
      type = 'l', lwd=2, ylim = c(0,100), col = my_pal[5])
lines(x=lo_methyl_150$x[order(lo_methyl_150$x)], y=lo_methyl_150$fitted[order(lo_methyl_150$x)], 
      type = 'l', lwd=2, ylim = c(0,100), col = my_pal[6])