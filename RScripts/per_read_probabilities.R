#!/bin/R

#########################
# Author: Molly Brothers
# Github: mollybrothers
# Date: 2020-11-08
#########################

# IN PROGRESS


library(data.table)
library(tidyverse)

dt <- fread("~/sequencing/201012_Doudna/per_read_modified_base_calls.txt")
sorted <- dt[order(read_id, pos)]

III_positive_strand <- sorted[chrm == "III" & strand == 1, 
                      list(read_id, chrm, pos, mod_log_prob)]

III_negative_strand <- sorted[chrm == "III" & strand == -1,
                      list(read_id, chrm, pos, mod_log_prob)]


# plot regions, using the average mod_log_prob for each adenine
plot_prob <- function(input, title) {
  ggplot(input, aes(x = pos, y = avg_prob)) +
    geom_point() +
    labs(
      title = sprintf("%s region", title),
      x = "chr III position",
      y = "average log_mod_prob") +
    ylim(-2, 0)
}

control <- sorted[chrm == "III" & pos %between% c(100e3, 105e3),
                  .(avg_prob = mean(mod_log_prob)),
                  by = pos
                  ]
plot_prob(control, "control")

HML <- III_positive_strand[pos %between% c(11e3, 16e3),
                           .(avg_prob = mean(mod_log_prob)),
                           by = pos]
plot_prob(HML, "HML")

HMR <- III_positive_strand[pos %between% c(291e3, 296e3),
                           .(avg_prob = mean(mod_log_prob)),
                           by = pos]

plot_prob(HMR, "HMR")