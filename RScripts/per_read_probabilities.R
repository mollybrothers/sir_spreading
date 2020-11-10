#!/bin/R

#########################
# Author: Molly Brothers
# Github: mollybrothers
# Date: 2020-11-09
#########################

# IN PROGRESS

library(data.table)
library(tidyverse)
library(wesanderson)

chr = "III" #which chromosome?
strand = "+ strand" #which strand?

#the per_read_base_calls.txt file itself is too large for RStudio's memory,
#so you'll need to use mawk on the command line to pick out the lines you want first
dt <- fread("/Volumes/brothers_seq/Nanopore/201012_Doudna/megalodon_outputv2/chrIII_plus.txt",
            select = c(1, 2, 4, 5),
            col.names = c("read_id", "chrm", "pos", "mod_log_prob"))

#convert log probabilities into probabilities
probs <- dt[mod_log_prob > -2, list(read_id, pos, mod_log_prob, 
                                    mod_prob = 10 ^ mod_log_prob)]

# plot average probabilities for each position
plot_prob <- function(input, title) {
  ggplot(input, aes(x = pos, y = avg_prob)) +
    geom_point() +
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
             .(read_id, pos, avg_prob = mean(mod_prob)),
             by = pos
             ]
plot_prob(HML_avg, "HML")

HMR_avg <- probs[pos %between% c(291e3, 296e3),
             .(read_id, pos, avg_prob = mean(mod_prob)),
             by = pos
             ]
plot_prob(HMR_avg, "HMR")

#plot single reads
control <- probs[pos %between% c(100e3, 105e3), list(read_id, pos, mod_prob)]

HML <- probs[pos %between% c(11e3, 16e3), list(read_id, pos, mod_prob)]
HML_bounds = c(11146, 14849)

HMR <- probs[pos %between% c(291e3, 296e3), list(read_id, pos, mod_prob)]
HMR_bounds = c(292388, 295034)

pal <- wes_palette("Zissou1", 100, type = "continuous")

ggplot(HMR, aes(x = pos, y = read_id, color = mod_prob)) +
  geom_point(shape = 15, size = 0.8) +
  theme(axis.text.y = element_blank()) +
  scale_color_gradientn(colors = pal) +
  geom_vline(xintercept = HMR_bounds)

                         