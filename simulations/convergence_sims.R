library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

data <- read_tsv(args[1],
                 col_names=c("tmp","r1","r2","loci")) %>%
  mutate(r1 = r1 - 0.75,
         r2 = r2 - 0.75)


covariance <- cov(data$r1,data$r2)
variance_1 <- var(data$r1)
variance_2 <- var(data$r2)
total_variance_plus_cov = variance_1 + variance_2 + covariance
total_variance= variance_1 + variance_2

mean_variance = sqrt(variance_1 * variance_2)



G_t = covariance/total_variance


convergence_correlation=covariance/mean_variance

percent_convergent = covariance/total_variance
percent_convergent
