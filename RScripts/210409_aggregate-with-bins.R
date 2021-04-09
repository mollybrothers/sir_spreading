dt <- fread("/Volumes/brothers_seq/201125_Turkey/megalodon_output_06/modified_bases.aggregate06.6mA.bed")
columns <- c("chrom", "start", "end", "name", "score", 
             "strand", "startCodon", "stopCodon", "color", 
             "coverage", "percentage")
colnames(dt) <- columns

select_cols <- c("chrom", "start", "coverage", "percentage")
relevant <- dt[chrom != "MT" & coverage > 10, ..select_cols]
HMR <- relevant[chrom == "III" & start > 292174 & start < 295364]
HMR_E = c(292674, 292769)
HMR_I = c(294805, 294864)

#read in the segments partitioned by methylation level from loess fitting to aggregate data
HMR_segments <- readRDS("~/sequencing/sir_spreading/data/segmentsHMR_extended.rds")

#add variables for methylation status and bin number to HMR_segments
meth_status <- rep_len(c("high", "low"), nrow(HMR_segments))
HMR_segments$meth_status <- meth_status
HMR_segments$bin <- seq(1, nrow(HMR_segments), 1)

#assign each of the rows in HMR to a bin
HMR$bin <- as.numeric()
for(y in 1:nrow(HMR_segments)){
  for(x in 1:nrow(HMR)){
    if (between(HMR[x,]$start, HMR_segments[y,]$start, HMR_segments[y,]$end)){
      HMR[x,]$bin <- HMR_segments[y,]$bin
    }
  }
}

HMR_bins <- HMR %>% group_by(bin) %>% summarize(start = start, average_meth = mean(percentage), sd = sd(percentage))

ggplot(HMR_bins, aes(x = start, y = average_meth)) +
  geom_point() +
  geom_line() +
  annotate("rect", xmin = c(HMR_E[1], HMR_I[1]), xmax = c(HMR_E[2], HMR_I[2]),
           fill = "black", ymin = 0, ymax = 100, alpha = 0.3) +
  annotate("rect", xmin = c(HMR_segments[which(HMR_segments$meth_status == "high"),]$start), xmax = c(HMR_segments[which(HMR_segments$meth_status == "high"),]$end),
           ymin = 0, ymax = 100, alpha = 0.3, fill = "mediumpurple4")
  #geom_errorbar(mapping = aes(ymin = average_meth-sd, ymax = average_meth+sd))
                