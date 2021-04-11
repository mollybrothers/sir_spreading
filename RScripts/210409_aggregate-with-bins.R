# IN PROGRESS

library(data.table)
library(tidyverse)
library(viridis)

columns <- c("chrom", "start", "end", "name", "score", 
             "strand", "startCodon", "stopCodon", "color", 
             "coverage", "percentage")

dt_1 <- fread("/Volumes/brothersseq/210304_Amanita/modified_bases.aggregate07.6mA.bed")
colnames(dt_1) <- columns

dt_2 <- fread("/Volumes/brothersseq/210304_Amanita/modified_bases.aggregate08.6mA.bed")
colnames(dt_2) <- columns

dt_3 <- fread("/Volumes/brothersseq/210304_Amanita/modified_bases.aggregate09.6mA.bed")
colnames(dt_3) <- columns

dt_4 <- fread("/Volumes/brothersseq/210304_Amanita/modified_bases.aggregate10.6mA.bed")
colnames(dt_4) <- columns

dt_5 <- fread("/Volumes/brothersseq/210304_Amanita/modified_bases.aggregate11.6mA.bed")
colnames(dt_5) <- columns

dt_6 <- fread("/Volumes/brothersseq/210304_Amanita/modified_bases.aggregate12.6mA.bed")
colnames(dt_6) <- columns

dt_ss <- fread("/Volumes/brothersseq/210403_Hello/modified_bases.aggregate02.6mA.bed")
colnames(dt_ss) <- columns

#get a new data.table only containing chrom, start, coverage, and percentage
#filter out mitochondrial DNA (MT) and coverage < 10
select_cols <- c("chrom", "start", "coverage", "percentage")
relevant_1 <- dt_1[chrom != "MT" & coverage > 10, ..select_cols]
relevant_2 <- dt_2[chrom != "MT" & coverage > 10, ..select_cols]
relevant_3 <- dt_3[chrom != "MT" & coverage > 10, ..select_cols]
relevant_4 <- dt_4[chrom != "MT" & coverage > 10, ..select_cols]
relevant_5 <- dt_5[chrom != "MT" & coverage > 10, ..select_cols]
relevant_6 <- dt_6[chrom != "MT" & coverage > 10, ..select_cols]
relevant_ss <- dt_ss[chrom != "MT" & coverage > 10, ..select_cols]

#assign each of the rows in methylation data frame to a bin
bin <- function(data, segments){
  data$bin <- as.numeric()
  for(y in 1:nrow(segments)){
    for(x in 1:nrow(data)){
      if (between(data[x,]$start, segments[y,]$start, segments[y,]$end)){
        data[x,]$bin <- segments[y,]$bin
      }
    }
  }
  return(data)
}

#### HMR ####
#read in the segments partitioned by methylation level from loess fitting to aggregate data
HMR_segments <- data.table(readRDS("~/sequencing/sir_spreading/data/segmentsHMR_extended.rds"))

#add variables for methylation status and bin number to HMR_segments
meth_status_HMR <- rep_len(c("high", "low"), nrow(HMR_segments))
HMR_segments$meth_status <- meth_status_HMR
HMR_segments$bin <- seq(1, nrow(HMR_segments), 1)
HMR_segments_high <- HMR_segments[(bin %% 2) != 0]
HMR_segments_low <- HMR_segments[(bin %% 2) == 0]

#HMR data
HMR_E = c(292674, 292769)
HMR_I = c(294805, 294864)

HMR_1 <- relevant_1[chrom == "III" & start > 292174 & start < 295364, list(start, percentage)]
HMR_2 <- relevant_2[chrom == "III" & start > 292174 & start < 295364, list(start, percentage)]
HMR_3 <- relevant_3[chrom == "III" & start > 292174 & start < 295364, list(start, percentage)]
HMR_4 <- relevant_4[chrom == "III" & start > 292174 & start < 295364, list(start, percentage)]
HMR_5 <- relevant_5[chrom == "III" & start > 292174 & start < 295364, list(start, percentage)]
HMR_6 <- relevant_6[chrom == "III" & start > 292174 & start < 295364, list(start, percentage)]
HMR_ss <- relevant_ss[chrom == "III" & start > 292174 & start < 295364, list(start, percentage)]

combined_HMR <- HMR_1 %>% full_join(HMR_2, by = "start", suffix = c("1", "2")) %>%
  full_join(HMR_3, by = "start", suffix = c("2", "3")) %>% full_join(HMR_4, by = "start", suffix = c("3", "4")) %>%
  full_join(HMR_5, by = "start", suffix = c("4", "5")) %>% full_join(HMR_6, by = "start", suffix = c("5", "6")) %>%
  full_join(HMR_ss, by = "start", suffix = c("6", "ss"))


combined_HMR <- bin(combined_HMR, HMR_segments_high)
combined_HMR <- combined_HMR[!is.na(bin)][order(bin)]
final_HMR <- combined_HMR %>% group_by(bin) %>% summarize(
  mean1 = mean(percentage1, na.rm = TRUE),
  mean2 = mean(percentage2, na.rm = TRUE),
  mean3 = mean(percentage3, na.rm = TRUE),
  mean4 = mean(percentage4, na.rm = TRUE),
  mean5 = mean(percentage5, na.rm = TRUE),
  mean6 = mean(percentage6, na.rm = TRUE),
  meanss = mean(percentage, na.rm = TRUE),
  start = start)
final_HMR$meanss = final_HMR$meanss - final_HMR$mean1
final_HMR$mean6 = final_HMR$mean6 - final_HMR$mean1
final_HMR$mean5 = final_HMR$mean5 - final_HMR$mean1
final_HMR$mean4 = final_HMR$mean4 - final_HMR$mean1
final_HMR$mean3 = final_HMR$mean3 - final_HMR$mean1
final_HMR$mean2 = final_HMR$mean2 - final_HMR$mean1
final_HMR$mean1 = final_HMR$mean1 - final_HMR$mean1

colors <- viridis(6)
ggplot(final_HMR, aes(x = start)) +
  geom_point(mapping = aes(y = mean1), color = colors[6]) +
  geom_point(mapping = aes(y = mean2), color = colors[5]) +
  geom_point(mapping = aes(y = mean3), color = colors[4]) +
  geom_point(mapping = aes(y = mean4), color = colors[3]) +
  geom_point(mapping = aes(y = mean5), color = colors[2]) +
  geom_point(mapping = aes(y = mean6), color = colors[1]) +
  geom_point(mapping = aes(y = meanss), color = "black") +
  geom_line(mapping = aes(y = mean1), color = colors[6]) +
  geom_line(mapping = aes(y = mean2), color = colors[5]) +
  geom_line(mapping = aes(y = mean3), color = colors[4]) +
  geom_line(mapping = aes(y = mean4), color = colors[3]) +
  geom_line(mapping = aes(y = mean5), color = colors[2]) +
  geom_line(mapping = aes(y = mean6), color = colors[1]) +
  geom_line(mapping = aes(y = meanss), color = "black") +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(size = 15, color = "black", family = "Arial")) +
  labs(x = "position on chr III",
       y = "average % methylated reads")+
  annotate("rect", xmin = c(HMR_E[1], HMR_I[1]), xmax = c(HMR_E[2], HMR_I[2]),
           fill = "black", ymin = 0, ymax = 70, alpha = 0.3)

#### HML ####
#read in the segments partitioned by methylation level from loess fitting to aggregate data
HML_segments <- data.table(readRDS("~/sequencing/sir_spreading/data/segmentsHML_extended.rds"))

#add variables for methylation status and bin number to HMR_segments
meth_status_HML <- rep_len(c("high", "low"), nrow(HML_segments))
HML_segments$meth_status <- meth_status_HML
HML_segments$bin <- seq(1, nrow(HML_segments), 1)
HML_segments_high <- HML_segments[(bin %% 2) != 0]
HML_segments_low <- HML_segments[(bin %% 2) == 0]

#HML data
HML_E <- c(11237, 11268)
HML_I <- c(14600, 14711)

HML_1 <- relevant_1[chrom == "III" & start > 10738 & start < 15211, list(start, percentage)]
HML_2 <- relevant_2[chrom == "III" & start > 10738 & start < 15211, list(start, percentage)]
HML_3 <- relevant_3[chrom == "III" & start > 10738 & start < 15211, list(start, percentage)]
HML_4 <- relevant_4[chrom == "III" & start > 10738 & start < 15211, list(start, percentage)]
HML_5 <- relevant_5[chrom == "III" & start > 10738 & start < 15211, list(start, percentage)]
HML_6 <- relevant_6[chrom == "III" & start > 10738 & start < 15211, list(start, percentage)]
HML_ss <- relevant_ss[chrom == "III" & start > 10738 & start < 15211, list(start, percentage)]

combined_HML <- HML_1 %>% full_join(HML_2, by = "start", suffix = c("1", "2")) %>%
  full_join(HML_3, by = "start", suffix = c("2", "3")) %>% full_join(HML_4, by = "start", suffix = c("3", "4")) %>%
  full_join(HML_5, by = "start", suffix = c("4", "5")) %>% full_join(HML_6, by = "start", suffix = c("5", "6")) %>%
  full_join(HML_ss, by = "start", suffix = c("6", "ss"))


combined_HML <- bin(combined_HML, HML_segments_high)
combined_HML <- combined_HML[!is.na(bin)][order(bin)]
final_HML <- combined_HML %>% group_by(bin) %>% summarize(
  mean1 = mean(percentage1, na.rm = TRUE),
  mean2 = mean(percentage2, na.rm = TRUE),
  mean3 = mean(percentage3, na.rm = TRUE),
  mean4 = mean(percentage4, na.rm = TRUE),
  mean5 = mean(percentage5, na.rm = TRUE),
  mean6 = mean(percentage6, na.rm = TRUE),
  meanss = mean(percentage, na.rm = TRUE),
  start = start)
final_HML$meanss = final_HML$meanss - final_HML$mean1
final_HML$mean6 = final_HML$mean6 - final_HML$mean1
final_HML$mean5 = final_HML$mean5 - final_HML$mean1
final_HML$mean4 = final_HML$mean4 - final_HML$mean1
final_HML$mean3 = final_HML$mean3 - final_HML$mean1
final_HML$mean2 = final_HML$mean2 - final_HML$mean1
final_HML$mean1 = final_HML$mean1 - final_HML$mean1

colors <- viridis(6)
ggplot(final_HML, aes(x = start)) +
  geom_point(mapping = aes(y = mean1), color = colors[6]) +
  geom_point(mapping = aes(y = mean2), color = colors[5]) +
  geom_point(mapping = aes(y = mean3), color = colors[4]) +
  geom_point(mapping = aes(y = mean4), color = colors[3]) +
  geom_point(mapping = aes(y = mean5), color = colors[2]) +
  geom_point(mapping = aes(y = mean6), color = colors[1]) +
  geom_point(mapping = aes(y = meanss), color = "black") +
  geom_line(mapping = aes(y = mean1), color = colors[6]) +
  geom_line(mapping = aes(y = mean2), color = colors[5]) +
  geom_line(mapping = aes(y = mean3), color = colors[4]) +
  geom_line(mapping = aes(y = mean4), color = colors[3]) +
  geom_line(mapping = aes(y = mean5), color = colors[2]) +
  geom_line(mapping = aes(y = mean6), color = colors[1]) +
  geom_line(mapping = aes(y = meanss), color = "black") +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(size = 15, color = "black", family = "Arial")) +
  labs(x = "position on chr III",
       y = "average % methylated reads")+
  annotate("rect", xmin = c(HML_E[1], HML_I[1]), xmax = c(HML_E[2], HML_I[2]),
           fill = "black", ymin = 0, ymax = 40, alpha = 0.3)
