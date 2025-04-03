########################.
## TREE DATA - COUNTS ##
########################.

## PACKAGES ####

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

#######################.

setwd("data")

#######################.

## DATA ####

data <- read_xlsx("trees.xlsx")


data_counts <- data %>%
  mutate(total = rowSums(select(., ends_with("_pct")))) %>%
  pivot_longer(contains("pct"), names_to = "type", values_to = "proportion") %>%
  mutate(type = case_when(type == "Ash_pct"    ~ "ash",
                          type == "Beech_pct"  ~ "beech",
                          type == "Larch_pct"  ~ "larch",
                          type == "Oak_pct"    ~ "oak",
                          type == "ScPine_pct" ~ "scots_pine",
                          type == "SBirch_pct" ~ "silver_birch",
                          type == "SitSpr_pct" ~ "sitka_spruce",
                          type == "SwChes_pct" ~ "sweet_chesnut",
                          type == "Syca_pct"   ~ "sycamore",
                          type == "Shadow_pct" ~ "shadow")) %>%
  mutate(count = (proportion * 100) / total) %>%
  pivot_wider(., names_from = "type", values_from = "count")


data_count <- data %>%
  mutate(total = rowSums(select(., ends_with("_pct")))) %>%
  mutate(across(ends_with("_pct"), ~ round(.x * 100 / total))) %>%
  mutate(total = rowSums(select(., ends_with("_pct")))) %>%
  
  rename(ash = "Ash_pct",
         beech = "Beech_pct",
         larch = "Larch_pct",
         oak = "Oak_pct",
         scots_pine = "ScPine_pct",
         silver_birch = "SBirch_pct",
         sitka_spruce = "SitSpr_pct",
         sweet_chesnut = "SwChes_pct",
         sycamore = "Syca_pct",
         shadow = "Shadow_pct")


#######################.

## HEATMAP ####

pdf("output/tree_count_heatmap_grid.pdf", height = 9)
data_count %>%
  pivot_longer(cols = ash:sycamore, names_to = "type", values_to = "count") %>%
  mutate(type = str_to_title(str_replace_all(type, "_", " "))) %>%
  ggplot(., aes(x = X_coord, y = Y_coord, fill = count)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  facet_wrap(~type, ncol = 2) +
  labs(x = "x coordinate", y = "y coordinate") +
  theme(legend.position = "bottom") #+
  #ggtitle("Heatmap for Tree Species")
dev.off()


#######################.

## SAVE DATA ####

saveRDS(data_count, "tree_count_100_data.rds")

#######################.
