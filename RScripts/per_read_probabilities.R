#!/bin/R

#########################
# Author: Molly Brothers
# Github: mollybrothers
# Date: 2021-02-14 <3
#########################

#######################################################################################
# the per_read_base_calls.txt file itself is too large for RStudio's memory,
# so you'll need to use mawk on the command line to pick out the lines you want first
# EXAMPLE:
# cat per_read_modified_base_calls.txt | mawk '$2 ~ /^III$/ {print $0}' > chrIII.txt
#######################################################################################

library(data.table)
library(tidyverse)
library(wesanderson)

mega_directory <- "/Volumes/brothers_seq/Nanopore/201125_Turkey/megalodon_output_07/"
chr <- "VI" #which chromosome?

# find all the reads with a qscore < 9 and filter out those reads from the input chr.txt file
summary <- fread(sprintf("%ssequencing_summary.txt", mega_directory),
                 select = c('read_id', 'mean_qscore_template'))

bad_reads <- summary[mean_qscore_template <= 9, read_id]

unfiltered <- fread(sprintf("%schr%s.txt", mega_directory, chr),
                    select = c(1, 2, 4, 5),
                    col.names = c("read_id", "chrm", "pos", "mod_log_prob"))

filtered <- unfiltered[!(read_id %in% bad_reads)]

#1. convert log probabilities into probabilities
#2. add a binary column; if mod_prob is >80, set as methylated (TRUE). if < 80, set as unmethylated (FALSE)
#3. add start and end positions for each read_id
#4. find the average methylation of each read
#5. order by average number of A's methylated then read_id
probs <- filtered[, list(read_id, pos, mod_prob = 10 ^ mod_log_prob)][
  , binary := ifelse(mod_prob > 0.8, "m6A", "A")][
    , start_pos := min(pos), by = read_id][
      , end_pos := max(pos), by = read_id][
        , avg_methyl := mean(binary == "m6A"), by = read_id][
          order(-avg_methyl, read_id)]

#extract each unique read_id to set the order (in same order as in the data table to start with) for single read plots
read_names <- unique(probs$read_id)
probs$read_id = factor(probs$read_id, levels = read_names)

##################################
####PLOT AVERAGE PROBABILITIES####
##################################

plot_prob <- function(input, title) {
  ggplot(input, aes(x = pos, y = avg_prob)) +
    theme_classic() +
    geom_point(position = "jitter") +
    labs(
      title = sprintf("%s region", title),
      x = sprintf("chr %s", chr),
      y = "average log_mod_prob") +
    ylim(0, 1)
}

control_avg <- probs[pos %between% c(100e3, 105e3),
                     .(avg_prob = mean(mod_prob)),
                     by = pos
]
plot_prob(control_avg, "control")

HML_avg <- probs[pos %between% c(11e3, 16e3),
                 .(avg_prob = mean(mod_prob)),
                 by = pos
]
plot_prob(HML_avg, "HML")

HMR_avg <- probs[pos %between% c(291e3, 296e3),
                 .(avg_prob = mean(mod_prob)),
                 by = pos
]
plot_prob(HMR_avg, "HMR")

#########################
####PLOT SINGLE READS####
#########################

#plot functions
probs_pal <- wes_palette("Zissou1", 2, type = "continuous")
binary_pal <- wes_palettes$Moonrise1

plot_probs <- function(data) {
  ggplot(data, aes(x = pos, y = read_id, color = mod_prob)) +
    geom_point(shape = 15) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.grid = element_blank()) +
    scale_color_gradientn(colors = probs_pal, name = "m6A probability")
}

plot_binary <- function(data) {
  ggplot(data, aes(x = pos, y = read_id, color = binary)) +
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
    scale_color_manual(values = c("grey90", "mediumpurple4"))
}

#plot unmethylated region
control_region <- c(95e3, 100e3)
control <- probs[start_pos <= control_region[1]][end_pos >= control_region[2]][
  pos %between% control_region]
plot_binary(control)
plot_probs(control)

#plot HML
HML_region <- c(10.5e3, 15.5e3)
HML_silencers <- c(11146, 14849)
HML_linkers = c(#9407, 9587.5, 9067, 9747, 9923, 10166, 10331, 
  10585, 10748,
  10944, 11118, 11418, 11645, 12021, 12251, 12436, 12649, 12842,
  13017, 13396, 13558, 13829, 14011, 14221, 14883, 15229, 15406
  #15573, 15984, 16244
)

HML <- probs[start_pos <= HML_region[1]][end_pos >= HML_region[2]][
  pos %between% HML_region]
hmlp <- plot_binary(HML)
hmlp
hmlp + geom_vline(xintercept = HML_silencers) #+
#geom_vline(xintercept = HML_linkers, size = 0.3, color = "black")
plot_probs(HML)

#plot HMR
HMR_region <- c(291e3, 296e3)
HMR_silencers = c(292388, 295034)
HMR_linkers = c(291100, 291312, 291644, 291863, 292129, 292322,
                292498, 292921, 293078, 293227, 293440, 293633, 293841, 294155,
                294515, 294699, 295239, 295555, 295743, 295906)

HMR <- probs[start_pos <= HMR_region[1]][end_pos >= HMR_region[2]][
  pos %between% HMR_region]
hmrp <- plot_binary(HMR)
hmrp
hmrp + geom_vline(xintercept = HMR_silencers)
plot_probs(HMR)

#plot tel3L
tel3L_region <- c(300, 5e3)
tel3L <- probs[start_pos <= tel3L_region[1]][end_pos >= tel3L_region[2]][
  pos %between% tel3L_region]
plot_binary(tel3L)

#plot tel6R
tel6R_region <- c(265e3, 270e3)
tel6R <- probs[start_pos <= tel6R_region[1]][end_pos >= tel6R_region[2]][pos %between% tel6R_region]
plot_binary(tel6R)