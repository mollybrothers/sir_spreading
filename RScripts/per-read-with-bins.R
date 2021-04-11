# IN PROGRESS: This seems like it won't be all that helpful

library(data.table)
library(tidyverse)
library(wesanderson)

bin <- function(data, segments){
  data$bin <- as.numeric()
  for(x in 1:nrow(data)){
    for(y in 1:nrow(segments)){
      if (between(data[x,]$pos, segments[y,]$start, segments[y,]$end)){
        data[x,]$bin <- segments[y,]$bin
      }
    }
  }
  return(data)
}

mega_directory <- "/Volumes/brothers_seq/201125_Turkey/megalodon_output_06/"
chr <- "III" #which chromosome?
barcode <- "06"

probs <- fread(sprintf("%schr%s_%s.txt", mega_directory, chr, barcode),
               select = c(1, 2, 3, 4, 5),
               col.names = c("read_id", "chrm", "strand", "pos", "mod_log_prob"))

#1. create a binary column; if 10^mod_log_prob is >0.8, set as methylated (m6A). if < 0.8, set as unmethylated (A)
#2. add start and end positions for each read_id
#3. find the average methylation of each read (for ordering on plot)
#4. order by strand, avg methyl and read_id
probs <- probs[, methylation := ifelse(10^mod_log_prob > 0.8, TRUE, FALSE)][
  , list(read_id, pos, methylation, strand)][
    , start_pos := min(pos), by = read_id][
      , end_pos := max(pos), by = read_id][
        , avg_methyl := mean(methylation == TRUE), by = read_id][
          order(strand, -avg_methyl, read_id)]

# INCLUDE NA's (methylation = NA if prob is between 0.2 and 0.8)
# probs <- probs[, methylation := ifelse(10^mod_log_prob > 0.8, TRUE,
#                                        ifelse(10^mod_log_prob < 0.2, FALSE, NA))][
#   , list(read_id, pos, methylation, strand)][
#     , start_pos := min(pos), by = read_id][
#       , end_pos := max(pos), by = read_id][
#         , avg_methyl := mean(methylation == TRUE, na.rm = TRUE), by = read_id][
#           order(strand, -avg_methyl, read_id)]

#extract each unique read_id to set the order (in same order as in the data table to start with) for single read plots
read_names <- unique(probs$read_id)
probs$read_id = factor(probs$read_id, levels = read_names)

#HMR
HMR_region <- c(291e3, 296e3)
HMR_E = c(292674, 292769)
HMR_I = c(294805, 294864)

HMR_all <- probs[start_pos <= HMR_region[1]][end_pos >= HMR_region[2]][
  pos %between% HMR_region]
HMR_all <- droplevels(HMR_all)

HMR_plus <- probs[strand == "+"][start_pos <= HMR_region[1]][end_pos >= HMR_region[2]][
  pos %between% HMR_region]
HMR_plus <- droplevels(HMR_plus)

HMR_minus <- probs[strand == "-"][start_pos <= HMR_region[1]][end_pos >= HMR_region[2]][
  pos %between% HMR_region]
HMR_minus <- droplevels(HMR_minus)

#segmenting HMR based on loess fit to aggregate methylation data
HMR_segments <- data.table(readRDS("~/sequencing/sir_spreading/data/segmentsHMR_extended.rds"))
meth_status_HMR <- rep_len(c("high", "low"), nrow(HMR_segments))
HMR_segments$meth_status <- meth_status_HMR
HMR_segments$bin <- seq(1, nrow(HMR_segments), 1)
HMR_segments_high <- HMR_segments[(bin %% 2) != 0]
HMR_segments_low <- HMR_segments[(bin %% 2) == 0]

HMR_binned <- bin(HMR_all, HMR_segments)
HMR_binned <- HMR_binned[!is.na(bin)]
HMR_binned_high_plus <- HMR_binned[(bin %% 2) != 0 & strand == "+"][,pmeth := mean(methylation), by = c("read_id", "bin")]

ggplot(test, aes(x = pos, y = read_id, color = thresh)) +
  geom_point(shape = 15) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 15, color = "black", family = "Arial"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key = element_blank()) +
  scale_color_manual(values = c("gray90", "mediumpurple4"))
