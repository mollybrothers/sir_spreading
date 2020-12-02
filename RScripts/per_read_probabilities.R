#!/bin/R

#########################
# Author: Molly Brothers
# Github: mollybrothers
# Date: 2020-12-01
#########################

#TO DO: figure out best way to plot 'binary' methylation calls for single reads...

#######################################################################################
# the per_read_base_calls.txt file itself is too large for RStudio's memory,
# so you'll need to use mawk on the command line to pick out the lines you want first
# EXAMPLE:
# cat per_read_modified_base_calls.txt | mawk '$2 ~ /^III$/ {print $0}' > chrIII.txt
#######################################################################################

library(data.table)
library(tidyverse)
library(wesanderson)

mega_directory <- "/Volumes/brothers_seq/Nanopore/201125_Turkey/megalodon_output_06/"
chr <- "III" #which chromosome?

# find all the reads with a qscore < 9 and filter out those reads from the input chr.txt file
dt <- fread(sprintf("%s/sequencing_summary.txt", mega_directory),
            select = c('read_id', 'mean_qscore_template'))

bad_reads <- dt[mean_qscore_template <= 9, read_id]

unfiltered <- fread(sprintf("%s/chr%s.txt", mega_directory, chr),
                    select = c(1, 2, 4, 5),
                    col.names = c("read_id", "chrm", "pos", "mod_log_prob"))

filtered <- unfiltered[!(read_id %in% bad_reads)]

#convert log probabilities into probabilities and order by position and read_id
probs <- filtered[mod_log_prob > -2, list(read_id, pos, mod_log_prob, 
                                    mod_prob = 10 ^ mod_log_prob)][
                                      order(pos, read_id)
                                    ]

#############################
#PLOT AVERAGE PROBABILITIES
#############################

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

####################
#PLOT SINGLE READS
####################

#extract each unique read_id to set the order for single read plots
read_names <- unique(probs$read_id)
probs$read_id = factor(probs$read_id, levels = read_names)

#plot functions
pal <- wes_palette("Zissou1", 2, type = "continuous")

plot_nucs <- function(data, bounds, linkers) {
  ggplot(data, aes(x = pos, y = read_id, color = mod_prob)) +
    geom_point(shape = 15, size = 0.5) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.grid = element_blank()) +
    scale_color_gradientn(colors = pal, name = "m6A probability") +
    geom_vline(xintercept = bounds) +
    geom_vline(xintercept = linkers, size = 0.3, color = "black")
}

plot_clean <- function(data, bounds) {
  ggplot(data, aes(x = pos, y = read_id, color = mod_prob)) +
    geom_point(shape = 15, size = 0.5) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.grid = element_blank()) +
    scale_color_gradientn(colors = pal, name = "m6A probability") +
    geom_vline(xintercept = bounds)
}

#plot unmethylated region
control <- probs[pos %between% c(180e3, 185e3), list(read_id, pos, mod_prob)]
ggplot(control, aes(x = pos, y = read_id, color = mod_prob)) +
  geom_point(shape = 15, size = 0.5) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank()) +
  scale_color_gradientn(colors = pal, name = "m6A probability")

#plot HML
HML <- probs[pos %between% c(10.5e3, 15.5e3), list(read_id, pos, mod_prob)]
HML_bounds = c(11146, 14849)
HML_linkers = c(#9407, 9587.5, 9067, 9747, 9923, 10166, 10331, 
                10585, 10748,
                10944, 11118, 11418, 11645, 12021, 12251, 12436, 12649, 12842,
                13017, 13396, 13558, 13829, 14011, 14221, 14883, 15229, 15406
                #15573, 15984, 16244
                )
plot_nucs(HML, HML_bounds, HML_linkers)
plot_clean(HML, HML_bounds)

#plot HMR
HMR <- probs[pos %between% c(291e3, 296e3), list(read_id, pos, mod_prob)]
HMR_bounds = c(292388, 295034)
HMR_linkers = c(291100, 291312, 291644, 291863, 292129, 292322,
                292498, 292921, 293078, 293227, 293440, 293633, 293841, 294155,
                294515, 294699, 295239, 295555, 295743, 295906)
plot_nucs(HMR, HMR_bounds, HMR_linkers)
plot_clean(HMR, HMR_bounds)

#plot tel3L
tel3L <- probs[pos %between% c(0, 5e3), list(read_id, pos, mod_prob)]
ggplot(tel3L, aes(x = pos, y = read_id, color = mod_prob)) +
  geom_point(shape = 15, size = 0.5) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c(0.1, 0.8),
        legend.direction = "vertical") +
  scale_color_gradientn(colors = pal, name = "m6A probability")                         

#plot tel3R
tel3R <- probs[pos %between% c(312e3, 317e3), list(read_id, pos, mod_prob)]
ggplot(tel3R, aes(x = pos, y = read_id, color = mod_prob)) +
  geom_point(shape = 15, size = 0.5) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c(0.1, 0.8),
        legend.direction = "vertical") +
  scale_color_gradientn(colors = pal, name = "m6A probability")

################
#PLOT BINARIES
################

#if mod_prob is > 80, set as methylated (1). if < 80, set as unmethylated (0).
binary <- probs[, binary := ifelse(mod_prob > 0.8, 1, 0)]

plot_clean_binaries <- function(data, bounds) {
  ggplot(data, aes(x = pos, y = read_id, color = binary)) +
    geom_point(shape = 15, size = 0.5) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.grid = element_blank()) +
    scale_color_gradientn(colors = pal, name = "m6A probability") +
    geom_vline(xintercept = bounds)
}

control <- binary[pos %between% c(180e3, 185e3), list(read_id, pos, binary)]
ggplot(control, aes(x = pos, y = read_id, color = binary)) +
  geom_point(shape = 15, size = 0.5) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank()) +
  scale_color_gradientn(colors = pal, name = "m6A probability")

#plot HML
HML <- binary[pos %between% c(10.5e3, 15.5e3), list(read_id, pos, binary)]
HML_bounds = c(11146, 14849)
HML_linkers = c(#9407, 9587.5, 9067, 9747, 9923, 10166, 10331, 
  10585, 10748,
  10944, 11118, 11418, 11645, 12021, 12251, 12436, 12649, 12842,
  13017, 13396, 13558, 13829, 14011, 14221, 14883, 15229, 15406
  #15573, 15984, 16244
)
plot_clean_binaries(HML, HML_bounds)

#plot HMR
HMR <- binary[pos %between% c(291e3, 296e3), list(read_id, pos, binary)]
HMR_bounds = c(292388, 295034)
HMR_linkers = c(291100, 291312, 291644, 291863, 292129, 292322,
                292498, 292921, 293078, 293227, 293440, 293633, 293841, 294155,
                294515, 294699, 295239, 295555, 295743, 295906)
plot_clean_binaries(HMR, HMR_bounds)
