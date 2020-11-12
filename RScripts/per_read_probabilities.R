#!/bin/R

#########################
# Author: Molly Brothers
# Github: mollybrothers
# Date: 2020-11-12
#########################

# IN PROGRESS

library(data.table)
library(tidyverse)
library(wesanderson)

chr = "III" #which chromosome?
strand = "+ strand" #which strand?

#the per_read_base_calls.txt file itself is too large for RStudio's memory,
#so you'll need to use mawk on the command line to pick out the lines you want first
#EXAMPLE:
#cat per_read_modified_base_calls.txt | mawk '$2 ~ /^III$/ {print $0}' > chrIII.txt
dt <- fread("/Volumes/brothers_seq/Nanopore/201012_Doudna/megalodon_outputv2/chrIII_per_read_modified_bases.txt",
            select = c(1, 2, 4, 5),
            col.names = c("read_id", "chrm", "pos", "mod_log_prob"))

#convert log probabilities into probabilities and order by position and read_id
probs <- dt[mod_log_prob > -2, list(read_id, pos, mod_log_prob, 
                                    mod_prob = 10 ^ mod_log_prob)][
                                      order(pos, read_id)
                                    ]

#extract each unique read_id to set the order for single read plots
read_names <- unique(probs$read_id)
probs$read_id = factor(probs$read_id, levels = read_names)

#plot average probabilities for each position
plot_prob <- function(input, title) {
  ggplot(input, aes(x = pos, y = avg_prob)) +
    theme_classic() +
    geom_point(position = "dodge") +
    labs(
      title = sprintf("%s region", title),
      x = sprintf("%s at chr %s position", strand, chr),
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

#plot functions for single reads
pal <- wes_palette("Zissou1", 100, type = "continuous")

plot_nucs <- function(data, bounds, linkers) {
  ggplot(data, aes(x = pos, y = read_id, color = mod_prob)) +
    geom_point(shape = 15, size = 0.5) +
    theme(axis.text.y = element_blank(), ) +
    scale_color_gradientn(colors = pal) +
    geom_vline(xintercept = bounds, size = 1) +
    geom_vline(xintercept = linkers, size = 0.3, color = "black")
}

plot_clean <- function(data, bounds) {
  ggplot(data, aes(x = pos, y = read_id, color = mod_prob)) +
    geom_point(shape = 15, size = 0.5) +
    theme(axis.text.y = element_blank(), ) +
    scale_color_gradientn(colors = pal) +
    geom_vline(xintercept = bounds, size = 1)
}

#plot unmethylated region
control <- probs[pos %between% c(180e3, 185e3), list(read_id, pos, mod_prob)]
ggplot(control, aes(x = pos, y = read_id, color = mod_prob)) +
  geom_point(shape = 15, size = 0.5) +
  theme(axis.text.y = element_blank())+
  scale_color_gradientn(colors = pal)

#plot HML
HML <- probs[pos %between% c(10e3, 15.5e3), list(read_id, pos, mod_prob)]
HML_bounds = c(11146, 14849)
HML_linkers = c(#9407, 9587.5, 9067, 9747, 9923, 10166, 10331, 10585, 10748,
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
  theme(axis.text.y = element_blank()) +
  scale_color_gradientn(colors = pal)                         

#plot tel3R
tel3R <- probs[pos %between% c(312e3, 317e3), list(read_id, pos, mod_prob)]
ggplot(tel3R, aes(x = pos, y = read_id, color = mod_prob)) +
  geom_point(shape = 15, size = 0.5) +
  theme(axis.text.y = element_blank()) +
  scale_color_gradientn(colors = pal)
