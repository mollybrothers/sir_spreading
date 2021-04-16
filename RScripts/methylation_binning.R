#!/bin/R

###############################################
# Author: Molly Brothers / Koen Van den Berge
# Github: mollybrothers
# Date: 2021-04-15
###############################################

# The purpose of this script is to create high methylation bins and low methylation bins
# based on steady-state (grown at 25C) methylation data. Create a loess curve and draw a threshold.
# The threshold is slightly different based on the genomic region.
# The output file will contain the start and stop of each bin, its designation as high or low
# methylation, and the bin number. Use this output file in the aggregate-with-bins.R script

library(data.table)
library(GenomicRanges)

# change for different samples
methyl <- fread("/Volumes/brothers_seq/210205_Ocular/megalodon_output_08/modified_bases.aggregate08.6mA.bed")
columns <- c("chrom", "start", "end", "name", "score", 
             "strand", "startCodon", "stopCodon", "color", 
             "coverage", "percentage")
colnames(methyl) <- columns
select_cols <- c("chrom", "start", "coverage", "percentage")
methyl_filtered <- methyl[coverage > 10, ..select_cols]

# uncomment for region of interest

# # for HMR
# save <- "HMR"
# chromo <- "III"
# beg <- 292674 - 500
# end <- 294864 + 500
# threshold <- 40
# odd_bins <- "low"
# even_bins <- "high"

#for HML
save <- "HML"
chromo <- "III"
beg <- 11237 - 500
end <- 14711 + 500
threshold <- 40
odd_bins <- "low"
even_bins <- "high"

# #tel13L
# save <- "tel13L"
# chromo <- "XIII"
# beg <- 0
# end <- 20000
# threshold <- 40
# odd_bins <- "high"
# even_bins <- "low"

#tel14L
# save <- "tel14L"
# chromo <- "XIV"
# beg <- 0
# end <- 20000
# threshold <- 28
# odd_bins <- "low"
# even_bins <- "high"

#### GET METHYLATION BINS ####
# the same for every sample
region_methyl <- methyl_filtered[chrom == chromo & start > beg & start < end]

#loess
lo_methyl <- loess(percentage ~ start, data=region_methyl, weights=region_methyl$coverage, enp.target = 100)
plot(x=lo_methyl$x[order(lo_methyl$x)], y=lo_methyl$fitted[order(lo_methyl$x)], 
     type = 'l', lwd=3)
abline(h=threshold, col="red")

#binning
rle <- Rle(lo_methyl$fitted[order(lo_methyl$x)] > threshold)
starts <- c(1, cumsum(rle@lengths)[-length(rle@lengths)]+1)
ends <- cumsum(rle@lengths)
segments <- data.table(start = lo_methyl$x[order(lo_methyl$x)][starts], 
                       end = lo_methyl$x[order(lo_methyl$x)][ends])
meth_status <- rep_len(c(odd_bins, even_bins), nrow(segments))
segments$meth_status <- meth_status
segments <- segments[,size := (end - start)][size > 10,][,list(start, end, meth_status)]
segments$bin <- seq(1, nrow(segments), 1)

plot(x=lo_methyl$x[order(lo_methyl$x)], y=lo_methyl$fitted[order(lo_methyl$x)], 
     type = 'l', lwd=3)
abline(h=threshold, col="red")
abline(v=segments$start, lty=2, col="orange")

saveRDS(segments, file=paste0("~/sequencing/sir_spreading/data/",save,"_bins.rds"))
