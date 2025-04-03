##############################.
## PLOT CLUSTER TIME SERIES ##
##############################.


## PACKAGES ##

library(dplyr)
library(ggplot2)

#######################.

## SET WORKING DIRECTORY ####

setwd("clustering")


#######################.

## DATA ####

variant_data <- read_rds("")

# CLUSTERS
countries_cluster <- cbind(countries = c("United Kingdom",
                                         "France",
                                         "Italy",
                                         "Germany",
                                         "Spain",
                                         "Australia",
                                         "Brazil",
                                         "Mexico",
                                         "South Africa",
                                         "Canada"), 
                           cluster = 3) %>%
  rbind(cbind(countries = c("Armenia",
                            "Azerbaijan",
                            "China",
                            "Dominica",
                            "Estonia",
                            "Fiji",
                            "India",
                            "Kenya",
                            "Qatar",
                            "United Arab Emirates"),
              cluster = 1)) %>%
  rbind(cbind(countries = c("Argentina",
                            "Bulgaria",
                            "Dominican Republic",
                            "Greece",
                            "Jamaica",
                            "Morocco",
                            "New Zealand",
                            "Pakistan",
                            "Peru",
                            "Philippines"),
              cluster = 2)) %>%
  as.data.frame() %>%
  arrange(countries)


cluster_variant_data <- variant_data %>%
  filter(country %in% countries_cluster$countries) %>%
  left_join(countries_cluster %>% rename(country = 'countries'), by = "country")


#######################.

## PLOT ####

pdf("output/covid_time_series-cluster1.pdf")
cluster_variant_data %>%
  filter(cluster == 1) %>%
  mutate(variant = factor(variant, levels = c("alpha", "beta", "gamma", "delta", "omicron", "VOI"))) %>%
  ggplot(., aes(x = date, y = count, colour = variant, group = variant)) +
  geom_line() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  scale_x_date(date_breaks = "6 months", 
               date_labels = "%b-%y") +
  #scale_colour_brewer(palette = "Set2") + 
  scale_color_manual(name = "Variant", values = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494")) +
  facet_wrap(~country, scales = "free_y") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 1)) +
  xlab("Week Date") +
  ylab("Cases")
dev.off()

pdf("output/covid_time_series-cluster2.pdf")
cluster_variant_data %>%
  filter(cluster == 2) %>%
  mutate(variant = factor(variant, levels = c("alpha", "beta", "gamma", "delta", "omicron", "VOI"))) %>%
  ggplot(., aes(x = date, y = count, colour = variant, group = variant)) +
  geom_line() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  scale_x_date(date_breaks = "6 months", 
               date_labels = "%b-%y") +
  #scale_colour_brewer(palette = "Set2") + 
  scale_color_manual(name = "Variant", values = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494")) +
  facet_wrap(~country, scales = "free_y") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 1)) +
  xlab("Week Date") +
  ylab("Cases")
dev.off()

pdf("output/covid_time_series-cluster3.pdf")
cluster_variant_data %>%
  filter(cluster == 3) %>%
  mutate(variant = factor(variant, levels = c("alpha", "beta", "gamma", "delta", "omicron", "VOI"))) %>%
  ggplot(., aes(x = date, y = count, colour = variant, group = variant)) +
  geom_line() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  scale_x_date(date_breaks = "6 months", 
               date_labels = "%b-%y") +
  #scale_colour_brewer(palette = "Set2") + 
  scale_color_manual(name = "Variant", values = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494")) +
  facet_wrap(~country, scales = "free_y") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 1)) +
  xlab("Week Date") +
  ylab("Cases")
dev.off()

#######################.
