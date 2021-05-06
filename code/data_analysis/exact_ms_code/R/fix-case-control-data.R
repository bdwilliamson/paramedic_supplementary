## Goal -- turn relative abundance table into count table
library(tidyverse)
rsm <- readRDS("RSMp01.Rds")

# Which columns are Prevotella?
prev <- rsm %>% names %>% startsWith("Prev") %>% which
(rsm %>% names)[prev]
(rsm %>% names)
rsm 

# What is the format of the data?
(rsm %>% names)[6:183]
rsm[, 7:182] %>% apply(1, sum) 
# Percentages, not even relative abundances

# Goal: turn into counts

# Start by making the assumption that lowest abundance taxon observed once
inv_min_ex_zero <- function(vec) 1/min(vec[vec>0])
relative_abundances <- rsm[, 7:182]/100
relative_abundances %>% apply(1, inv_min_ex_zero) # mostly whole numbers, good

ms <- relative_abundances %>% apply(1, inv_min_ex_zero)
median(ms) # paper says median 4345 reads, 4176 is close

# matrix of counts
ms_matrix <- matrix(ms, nrow = nrow(rsm), ncol = ncol(rsm[, 7:182]), byrow=F)

# initial guess at counts
ws <- relative_abundances * ms_matrix
ws

problematic_samples <- ws %>%
  as_tibble %>%
  bind_cols("id" = rsm[,"id"], .) %>%
  pivot_longer(cols=-id) %>%
  mutate(rounded = value %% 1) %>%
  mutate(rounded = round(rounded, 5)) %>%
  filter(rounded != 0) %>%
  filter(rounded != 1) 

problematic_samples$id %>% unique

### 500019782501
d500019782501 <- ws[ rsm[,"id"] == "500019782501", ] 
d500019782501[d500019782501 > 0]
ws[which(rsm[,"id"] == "500019782501"), ] <- 6*ws[which(rsm[,"id"] == "500019782501"), ]

# 500032992501
d500032992501 <- ws[ rsm[,"id"] == "500032992501", ] 
d500032992501[d500032992501 > 0]
ws[which(rsm[,"id"] == "500032992501"), ] <- 2*ws[which(rsm[,"id"] == "500032992501"), ]

# 500071142501 
d500071142501 <- ws[ rsm[,"id"] == "500071142501", ] 
d500071142501[d500071142501 > 0]
ws[which(rsm[,"id"] == "500071142501"), ] <- 7*ws[which(rsm[,"id"] == "500071142501"), ]

# 510042472501 
d510042472501 <- ws[ rsm[,"id"] == "510042472501", ] 
d510042472501[d510042472501 > 0]
ws[which(rsm[,"id"] == "510042472501"), ] <- 3*ws[which(rsm[,"id"] == "510042472501"), ]

# 520034122501
d520034122501 <- ws[ rsm[,"id"] == "520034122501", ] 
d520034122501[d520034122501 > 0]
ws[which(rsm[,"id"] == "520034122501"), ] <- 2*ws[which(rsm[,"id"] == "520034122501"), ]

# 550004152501
d550004152501 <- ws[ rsm[,"id"] == "550004152501", ] 
d550004152501[d550004152501 > 0]
# yikes
li <- relative_abundances[which(rsm[,"id"] == "550004152501"), "Lactobacillusiners"] 
lv <- relative_abundances[which(rsm[,"id"] == "550004152501"), "Lactobacillusvaginalis"] 

#  unlikely that this was the most deeply sequenced sample
tibble("scaling" = 1:ceiling(max(ms)/(li/lv))) %>%
  mutate("Lactobacillusiners" = li/lv * scaling) %>%
  mutate(iners = Lactobacillusiners %% 1) %>%
  filter(iners == min(iners))
ws[which(rsm[,"id"] == "550004152501"), ] <- 19*ws[which(rsm[,"id"] == "550004152501"), ]


ws %>%
  as_tibble %>%
  bind_cols("id" = rsm[,"id"], .) %>%
  pivot_longer(cols=-id) %>%
  mutate(rounded = value %% 1) %>%
  mutate(rounded = round(rounded, 5)) %>%
  filter(rounded != 0) %>%
  filter(rounded != 1) 
# good 


w_counts <- ws %>%
  as_tibble %>%
  bind_cols("id" = rsm[,"id"], .) %>%
  pivot_longer(cols=-id) %>%
  mutate(value = round(value)) %>%
  pivot_wider 

w_counts %>% select(-id) %>% apply(1, sum) %>% hist(breaks=20) 
w_counts %>% select(-id) %>% apply(1, sum) %>% sort
w_counts %>% select(-id) %>% apply(1, sum) %>% median # because of tie



all_data <- w_counts %>%
  right_join(rsm[, -c(7:182)], .)
all_data


# write to file
all_data %>%
  write_csv("RSMp01_counts.csv")
