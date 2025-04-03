#################################################.
## COVID DATA - CLUSTERING SPLINE COEFFICIENTS ##
#################################################.


## PACKAGES ##

library(dplyr)
library(cluster) 
library(factoextra)
library(dendextend)
library(openxlsx)
library(tidyr)


#######################.

## SET WORKING DIRECTORY ##

setwd("clustering")

#######################.

## DATA ####

data <- readRDS(file = "") %>%
  select_if(~ !is.numeric(.) || sum(.) != 0) %>%
  mutate(variant = forcats::fct_relevel(variant, c("alpha", "beta", "delta", "gamma", "omicron", "VOI", "all")))


# SQRT DATA
sqrt_data <- data %>%
  mutate(across(starts_with("20"), ~ if_else(variant != "all", sqrt(.x), .x)))

cluster_data <- readRDS(file = "spline_coeff_data.RDS") %>%
  # replace NA with 0
  replace(is.na(.), 0)


#######################.

## CLUSTERING ####

set.seed(123)


## DISTANCE MATRIX
dist <- cluster_data %>%
  select(-country) %>%
  dist(., method = "euclidean")


## OPTIMAL CLUSTERS
# ELBOW
cluster_data %>%
  select(-country) %>%
  fviz_nbclust(., FUN = hcut, method = "wss")

png("output/elbow.png")
cluster_data %>%
  select(-country) %>%
  fviz_nbclust(., FUN = hcut, method = "wss")
dev.off()


# 3 CLUSTERS
k = 3


# WARD LINKAGE
hc <- hclust(dist, method = "ward.D" )


# PLOT
plot(hc)


# BOXES FOR CLUSTERS
hc %>%
  as.dendrogram() %>%
  set("labels_cex", 0.8) %>%
  set("branches_k_color", k = k) %>% 
  plot(main = "Dendrogram - Ward Linkage, Coloured Clusters") 


# SQRT
png("output/dend_colour.png")
hc %>%
  as.dendrogram() %>%
  set("labels_cex", 0.8) %>%
  set("branches_k_color", k = k) %>% 
  plot(main = "Dendrogram - Ward Linkage, Coloured Clusters") 
dev.off()


# CUT TREE
hc_cut <- cutree(hc, k = k)


# NUMBER IN EACH CLUSTER
table(hc_cut)


# CLUSTERED COUNTRIES
clusters <- cluster_data %>%
  select(country) %>%
  mutate(cluster = hc_cut)


data_clusters <- sqrt_data %>%
  left_join(clusters, by = "country") %>%
  select(country,
         variant,
         measure,
         cluster,
         starts_with("20"))


#######################.

## OUTPUT CLUSTERED COUNTRIES ##

clusters %>%
  arrange(cluster, country) %>%
  write.csv("output/clustered_countries.csv")


#######################.


# MAP CLUSTERS ####

clustered_map = map_data("world") %>% 
  
  filter(region != "Antarctica") %>%
  
  mutate(region = if_else(region == "UK" | region == "Isle of Man" | region == "Guernsey" | region == "Jersey",
                          "United Kingdom", region)) %>%
  mutate(region = if_else(region == "Antigua" | region == "Barbuda", "Antigua and Barbuda", region)) %>%
  mutate(region = if_else(region == "Bahamas", "The Bahamas", region)) %>%
  mutate(region = if_else(region == "Ivory Coast", "Cote d Ivoire", region)) %>%
  mutate(region = if_else(region == "Republic of Congo", "Republic of the Congo", region)) %>%
  mutate(region = if_else(region == "Cabo Verde", "Cape Verde", region)) %>%
  mutate(region = if_else(region == "Guinea Bissau", "Guinea-Bissau", region)) %>%
  mutate(region = if_else(region == "Saint Kitts" | region == "Nevis", "Saint Kitts and Nevis", region)) %>%
  mutate(region = if_else(region == "Madeira Islands", "Portugal", region)) %>%
  mutate(region = if_else(region == "Timor-Leste", "Timor Leste", region)) %>%
  mutate(region = if_else(region == "Trinidad" | region == "Tobago", "Trinidad and Tobago", region)) %>%
  mutate(region = if_else(region == "Saint Vincent" | region == "Grenadines", "Saint Vincent and the Grenadines", region)) %>%
  mutate(region = if_else(region == "Virgin Islands", "British Virgin Islands", region)) %>%
  mutate(region = if_else(region == "Wallis and Futuna", "Wallis and Futuna Islands", region)) %>%
  mutate(region = if_else(region == "Samoa", "American Samoa", region)) %>%
  
  rename(country = "region") %>%
  left_join(clusters, by = "country") %>%
  rename(region = "country") %>%
  mutate(cluster = if_else(is.na(cluster), 9, cluster))


# NO DATA
no_data <- clustered_map %>%
  filter(cluster == 9) %>%
  select(region) %>%
  group_by(region) %>%
  distinct() %>%
  ungroup()

# NO CLUSTER
no_cluster <- clustered_map %>%
  filter(is.na(cluster)) %>%
  select(region) %>%
  group_by(region) %>%
  distinct() %>%
  ungroup()


ggplot(data = clustered_map, 
       aes(x = long, y = lat, group = group)) +
  geom_polygon() +
  geom_map(data = clustered_map, map = clustered_map, aes(map_id = region, fill = as.factor(cluster))) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", #"#A6D854",
                               "gray"),
                    labels = c("Cluster 1", "Cluster 2", "Cluster 3", 
                               "No Data")) +
  theme(legend.title = element_blank(),
        axis.text = element_blank(), 
        axis.ticks = element_blank()) + 
  theme(legend.position = "bottom", 
        legend.direction = "horizontal") +
  guides(fill = guide_legend(nrow = 1)) +
  ggtitle("Clustered World Map") 


pdf("output/clustered_map.pdf")
ggplot(data = clustered_map, 
       aes(x = long, y = lat, group = group)) +
  geom_polygon() +
  geom_map(data = clustered_map, map = clustered_map, aes(map_id = region, fill = as.factor(cluster))) +
  scale_fill_manual(name = "",
                    values = c("#66C2A5", "#FC8D62", "#8DA0CB", 
                               "gray"),
                    labels = c("Cluster 1", "Cluster 2", "Cluster 3", 
                               "No Data")) +
  theme(legend.title = element_blank(),
        axis.text = element_blank(), 
        axis.ticks = element_blank()) +
  theme_void() +
  labs(x = "", y = "") +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(fill = guide_legend(nrow = 1))
dev.off()


#######################.

## SAVE CLUSTERED DATA ##

data_clusters %>%
  saveRDS(., file = "output/clustered_spline_data_counts.RDS")

#######################.
