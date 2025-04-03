#######################.
## PLOTS FOR CHAPTER ##
#######################.


## PACKAGES ####

library(dplyr)
library(tidyr)
library(abind)
library(reshape2)
library(multidplyr)
library(ggplot2)
library(reshape2)
library(ggh4x)


#######################.


## LOAD FULL DATA ####

full_data <- readRDS(file = "~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Data/variant_data_extract-long.rds")


data_plot <- full_data %>%
  filter(country %in% c("United Kingdom",
                        "New Zealand",
                        "United Arab Emirates",
                        "Argentina",
                        "Kenya",
                        "Canada"
                        ))


pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD/R/COVID Variant/plots/continent_countries_time_series.pdf", 
    height = 5)
data_plot %>%
  mutate(variant = factor(variant, levels = c("alpha",
                                              "beta",
                                              "gamma",
                                              "delta",
                                              "omicron",
                                              "VOI"))) %>%
  
  ggplot(., aes(x = date, y = count, group = variant, col = variant)) +
  geom_line(size = 1) +
  scale_color_brewer(palette = "Set2") +
  scale_x_date(date_breaks = "6 month", date_labels = "%b %Y") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  facet_wrap(~country, scales = "free_y") +
  guides(color = guide_legend(nrow = 1)) + 
  theme(legend.position = "bottom",
        legend.box = "horizontal",   
        legend.direction = "horizontal",
        legend.margin = margin(t = -15)) +
  xlab("") +
  ylab("") +
  labs(color = "Variant") 
dev.off()


pdf("/Users/student/Documents/PhD/R/COVID Variant/plots/uk_count_time_series.pdf", width = 8, height = 6)
data_plot %>%
  filter(country == "United Kingdom") %>%
  mutate(variant = factor(variant, levels = c("alpha",
                                              "beta",
                                              "gamma",
                                              "delta",
                                              "omicron",
                                              "VOI"))) %>%
  
  ggplot(., aes(x = date, y = count, group = variant, col = variant)) +
  geom_line(size = 1) +
  scale_color_brewer(palette = "Set2") +
  scale_x_date(date_breaks = "6 month", date_labels = "%b %Y") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  #facet_wrap(~country, scales = "free_y") +
  guides(color = guide_legend(nrow = 1)) + 
  theme(legend.position = "bottom",
        legend.box = "horizontal",   
        legend.direction = "horizontal") +
  xlab("Week") +
  ylab("Cases") +
  labs(color = "Variant") 
dev.off()


pdf("/Users/student/Documents/PhD/R/COVID Variant/plots/uk_proportion_time_series.pdf", width = 8, height = 6)
data_plot %>%
  filter(country == "United Kingdom") %>%
  mutate(variant = factor(variant, levels = c("alpha",
                                              "beta",
                                              "gamma",
                                              "delta",
                                              "omicron",
                                              "VOI"))) %>%
  mutate(proportion = if_else(total == 0, 0, count / total)) %>%
  
  ggplot(., aes(x = date, y = proportion, group = variant, col = variant)) +
  geom_line(size = 1) +
  scale_color_brewer(palette = "Set2") +
  scale_x_date(date_breaks = "6 month", date_labels = "%b %Y") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  #facet_wrap(~country, scales = "free_y") +
  guides(color = guide_legend(nrow = 1)) + 
  theme(legend.position = "bottom",
        legend.box = "horizontal",   
        legend.direction = "horizontal") +
  scale_y_continuous(labels=scales::label_percent()) +
  xlab("Week") +
  ylab("Percentage of Cases") +
  labs(color = "Variant") 
dev.off()


pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD/R/COVID Variant/plots/uk_time_series.pdf", 
    height = 3)
data_plot %>%
  filter(country == "United Kingdom") %>%
  mutate(variant = factor(variant, levels = c("alpha",
                                              "beta",
                                              "gamma",
                                              "delta",
                                              "omicron",
                                              "VOI"))) %>%
  mutate(proportion = if_else(total == 0, 0, count / total)) %>%
  pivot_longer(., cols = c(count, proportion), names_to = 'measure', values_to = 'value') %>%
  mutate(value = if_else(measure == "proportion", value * 100, value)) %>%
  
  ggplot(., aes(x = date, y = value, group = variant, col = variant)) +
  geom_line(size = 1) +
  scale_color_brewer(palette = "Set2") +
  scale_x_date(date_breaks = "6 month", date_labels = "%b %Y") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  
  facet_wrap(~measure, scales = "free_y",  
             labeller = as_labeller(c(count = "Cases", proportion = "Percentage of Cases"))) +
  
  guides(color = guide_legend(nrow = 1)) + 
  theme(legend.position = "bottom",
        legend.box = "horizontal",   
        legend.direction = "horizontal",
        legend.margin = margin(t = -15)) +
  xlab("") +
  ylab("") +
  labs(color = "Variant") 
dev.off()



#######################.

## CLUSTER MAP ####

clusters <- read.csv("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD/R/COVID Variant/clustering/splines/data/clustered_countries_log.csv") %>%
  select(-X)


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


# pdf("/Users/student/Documents/PhD/R/COVID Variant/plots/spline_clustered_map-k3.pdf", width = 8, height = 6)
# ggplot(data = clustered_map,
#        aes(x = long, y = lat, group = group)) +
#   geom_polygon(color = "white") +
#   geom_map(data = clustered_map, map = clustered_map, aes(map_id = region, fill = as.factor(cluster))) +
#   scale_fill_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB", 
#                                "gray"),
#                     labels = c("Cluster 1", "Cluster 2", "Cluster 3", 
#                                "No Data")) +
#   #theme(legend.title = element_blank()) +
#   #coord_map() +
#   theme_void() +
#   theme(legend.position = "bottom", legend.direction = "horizontal") +
#   guides(fill = guide_legend(nrow = 1)) +
#   labs(fill = "Cluster") 
# dev.off()


pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD/R/COVID Variant/plots/spline_clustered_map-k3.pdf", 
    height = 5)
ggplot(data = clustered_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(aes(fill = as.factor(cluster)), color = "white") +  
  scale_fill_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB", "gray"),
                    labels = c("Cluster 1", "Cluster 2", "Cluster 3", "No Data")) +
  coord_quickmap() + 
  theme_void() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.title = element_blank()) + 
  guides(fill = guide_legend(nrow = 1)) 
dev.off()


#######################.

## LOAD IN DATA ####

# VARIANTS
variants <- c("alpha", 
              "beta",
              "gamma", 
              "delta",
              "omicron")

# COUNTRIES / CLUSTERS
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



variant_data <- readRDS(file = "~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Data/variant_data_extract-proportion_total.rds") %>%
  select_if(~ !is.numeric(.) || sum(.) != 0) 


data <- variant_data %>%
  filter(
    # COUNTRIES
    country %in% countries_cluster$countries,
    # VARIANTS
    variant %in% variants) %>%
  arrange(country, match(variant, variants))


data_long <- data %>%
  pivot_longer(starts_with("20"), names_to = "date", values_to = "n") %>%
  mutate(date = as.Date(date))

# DATES IN DATA
dates <- data_long %>%
  select(date) %>%
  distinct() %>%
  t() %>%
  as.Date()
dates


# NUMBER OF VARIANTS
N_variants = length(unique(variants))

# NUMBER OF COUNTRIES
N_countries = length(unique(countries_cluster$countries))

# NUMBER OF CLUSTERS
N_clusters = length(unique(countries_cluster$cluster))

# NUMBER OF TIME POINTS
N = length(dates)


#######################.

## CLUSTER TIME SERIES ####

data_clusters <- full_data %>%
  filter(country %in% countries_cluster$countries) %>%
  left_join(countries_cluster %>% rename(country = 'countries'), by = 'country')


pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD/R/COVID Variant/plots/countries_time_series_cluster_1.pdf", 
    height = 4.5)
data_clusters %>%
  filter(cluster == 1) %>%
  mutate(variant = factor(variant, levels = c("alpha",
                                              "beta",
                                              "gamma",
                                              "delta",
                                              "omicron",
                                              "VOI"))) %>%
  
  ggplot(., aes(x = date, y = count, group = variant, col = variant)) +
  geom_line(size = 1) +
  scale_color_brewer(palette = "Set2") +
  scale_x_date(date_breaks = "6 month", date_labels = "%b %Y") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  facet_wrap(~country, scales = "free_y", ncol = 5) +
  guides(color = guide_legend(nrow = 1)) + 
  theme(legend.position = "bottom",
        legend.box = "horizontal",   
        legend.direction = "horizontal",
        legend.margin = margin(t = -15),
        strip.text = element_text(size = 7)) +
  xlab("") +
  ylab("") +
  labs(color = "Variant") 
dev.off()

pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD/R/COVID Variant/plots/countries_time_series_cluster_2.pdf", 
    height = 4.5)
data_clusters %>%
  filter(cluster == 2) %>%
  mutate(variant = factor(variant, levels = c("alpha",
                                              "beta",
                                              "gamma",
                                              "delta",
                                              "omicron",
                                              "VOI"))) %>%
  
  ggplot(., aes(x = date, y = count, group = variant, col = variant)) +
  geom_line(size = 1) +
  scale_color_brewer(palette = "Set2") +
  scale_x_date(date_breaks = "6 month", date_labels = "%b %Y") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  facet_wrap(~country, scales = "free_y", ncol = 5) +
  guides(color = guide_legend(nrow = 1)) + 
  theme(legend.position = "bottom",
        legend.box = "horizontal",   
        legend.direction = "horizontal",
        legend.margin = margin(t = -15),
        strip.text = element_text(size = 7)) +
  xlab("") +
  ylab("") +
  labs(color = "Variant") 
dev.off()

pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD/R/COVID Variant/plots/countries_time_series_cluster_3.pdf", 
    height = 4.5)
data_clusters %>%
  filter(cluster == 3) %>%
  mutate(variant = factor(variant, levels = c("alpha",
                                              "beta",
                                              "gamma",
                                              "delta",
                                              "omicron",
                                              "VOI"))) %>%
  
  ggplot(., aes(x = date, y = count, group = variant, col = variant)) +
  geom_line(size = 1) +
  scale_color_brewer(palette = "Set2") +
  scale_x_date(date_breaks = "6 month", date_labels = "%b %Y") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  facet_wrap(~country, scales = "free_y", ncol = 5) +
  guides(color = guide_legend(nrow = 1)) + 
  theme(legend.position = "bottom",
        legend.box = "horizontal",   
        legend.direction = "horizontal",
        legend.margin = margin(t = -15),
        strip.text = element_text(size = 7)) +
  xlab("") +
  ylab("") +
  labs(color = "Variant") 
dev.off()


#######################.

## TRANSITION MATRIX ####

load("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD/R/COVID Variant/GDM/count model/mcmc_output.RData")

output <- as_tibble(do.call('rbind', samples))

p_df <- output%>%select(starts_with("p["))%>%
  unlist()%>%
  array(dim=c(nrow(output),5,5,5,3))%>%
  melt(varnames=c("sample","i","j","variant","cluster"))

kappa_df <- p_df%>%filter(i==j,i!=5,i!=1)%>%
  mutate(kappa=1/(1-value),state=i)%>%
  group_by(state,variant,cluster)%>%
  summarise(lower=quantile(kappa,0.025),
            median=median(kappa),
            upper=quantile(kappa,0.975))


pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD/R/COVID Variant/plots/transition_prob_plot.pdf", 
    height = 3)
kappa_df %>%
  mutate(variant = case_when(variant == 1 ~ "alpha",
                             variant == 2 ~ "beta",
                             variant == 3 ~ "gamma",
                             variant == 4 ~ "delta",
                             variant == 5 ~ "omicron"),
         cluster=as.factor(cluster),
         state = factor(case_match(state,
                                   1 ~ "Dormant State (State 1)",
                                   2 ~ "Increasing State (State 2)",
                                   3 ~ "Dominant State (State 3)",
                                   4 ~ "Decreasing State (State 4)"),
                        levels=c("Dormant State (State 1)",
                                 "Increasing State (State 2)",
                                 "Dominant State (State 3)",
                                 "Decreasing State (State 4)"))) %>%
  mutate(variant = factor(variant, levels = c("alpha",
                                              "beta",
                                              "gamma",
                                              "delta",
                                              "omicron",
                                              "VOI"))) %>%
  
  ggplot(.,aes(x=variant, col = cluster))+
  facet_wrap(~state)+
  #geom_errorbar(aes(ymin=lower,ymax=upper),alpha=0.5)+
  geom_point(aes(y=median,shape=cluster),size=2.5) +
  scale_color_brewer(palette = "Set2") +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))
  labs(col = "Cluster", 
       shape = "Cluster",
       x = NULL,
       #y = "Expected persistence length (weeks)") +
       y = "") +
  #ggtitle("Expected persistence lengths for the Active HMM States (States 2,3,4) by Variant and Cluster") +
  theme(legend.position = "bottom")
dev.off()


#######################.

## LOAD WINDOW SUMMARY ####

window_summary <- readRDS("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD/R/COVID Variant/comparison/hmm_rw_dlm/15_March/window_summary_hmm_rw_dlm.rds")


#######################.

## DATA ####

# X - VARIANT COUNTS
x <- data %>%
  filter(measure == "count") %>%
  split(.$country) %>%
  lapply(., function(df) {
    df %>%
      select(variant, starts_with("20")) %>%
      t() %>%
      as.data.frame() %>%
      janitor::row_to_names(row_number = 1) %>%
      mutate(across(everything(), as.numeric))
  })
x

add_missing_variants <- function(df) {
  missing_variants <- setdiff(variants, colnames(df))
  if(length(missing_variants) > 0) {
    df[, missing_variants] <- 0
  }
  df <- select(df, alpha, beta, gamma, delta, omicron)
  return(df)
}
x <- lapply(x, add_missing_variants)
x


# Y = TOTAL COUNT
y <- data %>%
  filter(measure == "total") %>%
  split(.$country) %>%
  lapply(., function(df) {
    df %>%
      select(variant, starts_with("20")) %>%
      t() %>%
      as.data.frame() %>%
      janitor::row_to_names(row_number = 1) %>%
      mutate(across(everything(), as.numeric)) %>%
      apply(., 1, max)
  })
y


y_total_count_array <- abind(y, along = 2) %>%
  melt(varnames = c("date", "m"), value.name = "y") %>%
  group_by(m) %>%
  mutate(t = row_number()) %>%
  ungroup()
y_total_count_array


### SUMMARISE DATA ###


cluster <- new_cluster(4)

# WINDOW WIDTH
L <- c(5, 10, 15, 20, 30) 


cluster_assign(cluster, window_stats=function(x,L,N,stat_fun){
  w_starts <- 1:(N-L+1)
  values <- sapply(w_starts,function(u)stat_fun(x[u:(u+L-1)]))
  return(values)
},L=L,N=N)



data_proportion <- abind(x,along=3)%>%
  
  melt(varnames=c("date","v","m"),value.name = "x") %>%
  
  # ADD Y - TOTAL COUNT
  left_join(y_total_count_array %>% select(date, m, y), by = c("date", "m")) %>%
  rename(t = "date") %>%
  
  group_by(t, m)%>%
  partition(cluster)%>%
  
  # PROPORTION
  mutate(proportion = dplyr::if_else(y == 0, 0, x/y)) %>%
  
  ungroup() %>%
  collect()


data_stats <- data_proportion %>%
  
  group_by(v,m)%>%
  partition(cluster)%>%
  
  summarise(total = sum(x),
            
            mean_count=mean(x),
            median_count=median(x),
            
            var_count=var(x),
            IQR_count=IQR(x),
            
            mean_SD_count=mean(window_stats(x,L[3],N,sd)),
            SD_upper_q_count=quantile(window_stats(x,L[3],N,sd), 0.75),
            SD_lower_q_count=quantile(window_stats(x,L[3],N,sd), 0.25),
            
            mean_prop=mean(proportion),
            median_prop=median(proportion),
            
            var_prop=var(proportion),
            IQR_prop=IQR(proportion),
            
            mean_SD_prop=mean(window_stats(proportion,L[3],N,sd)),
            SD_upper_q_prop=quantile(window_stats(proportion,L[3],N,sd), 0.75),
            SD_lower_q_prop=quantile(window_stats(proportion,L[3],N,sd), 0.25),)%>%
  collect()


#######################.

## PLOT ####

### DENSITY ####

#### MEAN SD PROPORTION ####
pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Comparison/output/density_log_mean_SD_proportion_UAE.pdf", width = 8, height = 6)
ggplot(window_summary %>%
         filter(window_length == 15,
                measure == "mean_SD_prop",
                m == "United Arab Emirates") %>% 
         mutate(method = factor(method, levels = c("rw",
                                                   "dlm",
                                                   "hmm"))) %>%
         mutate(v = factor(v, levels = c("alpha",
                                         "beta",
                                         "gamma",
                                         "delta",
                                         "omicron"))) %>%
         
         mutate(value = if_else(value == 0, 0.001, value)),
       
       aes(x = log(value), fill = method)) +
  
  geom_density(alpha = 0.4)+
  
  geom_vline(data = filter(data_stats, 
                           m == "United Arab Emirates"),
             aes(xintercept = log(mean_SD_prop)),
             col ="#FC8D62", linewidth = 0.8)+
  
  facet_wrap(~v,scale="free")+
  scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                "rw" = "#8DA0CB",
                                                "dlm" = "#E78AC3")) +
  theme(legend.position = "bottom") +
  labs(#title = "Density Plot of the Proportion Mean of the SD across Windows - United Arab Emirates",
       x = "Log Proportion Mean SD") 
dev.off()

pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Comparison/output/density_log_mean_SD_proportion_NZ.pdf", width = 8, height = 6)
ggplot(window_summary %>%
         filter(window_length == 15,
                measure == "mean_SD_prop",
                m == "New Zealand") %>% 
         mutate(method = factor(method, levels = c("rw",
                                                   "dlm",
                                                   "hmm"))) %>%
         mutate(v = factor(v, levels = c("alpha",
                                         "beta",
                                         "gamma",
                                         "delta",
                                         "omicron"))) %>%
         
         mutate(value = if_else(value == 0, 0.001, value)),
       
       aes(x = log(value), fill = method)) +
  
  geom_density(alpha = 0.4)+
  
  geom_vline(data = filter(data_stats, 
                           m == "New Zealand"),
             aes(xintercept = log(mean_SD_prop)),
             col ="#FC8D62", linewidth = 0.8)+
  
  facet_wrap(~v,scale="free")+
  scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                "rw" = "#8DA0CB",
                                                "dlm" = "#E78AC3")) +
  theme(legend.position = "bottom") +
  labs(#title = "Density Plot of the Proportion Mean of the SD across Windows - New Zealand",
       x = "Log Proportion Mean SD") 
dev.off()

pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Comparison/output/density_log_mean_SD_proportion_UK.pdf", width = 8, height = 6)
ggplot(window_summary %>%
         filter(window_length == 15,
                measure == "mean_SD_prop",
                m == "United Kingdom") %>% 
         mutate(method = factor(method, levels = c("rw",
                                                   "dlm",
                                                   "hmm"))) %>%
         mutate(v = factor(v, levels = c("alpha",
                                         "beta",
                                         "gamma",
                                         "delta",
                                         "omicron"))) %>%
         
         mutate(value = if_else(value == 0, 0.001, value)),
       
       aes(x = log(value), fill = method)) +
  
  geom_density(alpha = 0.4)+
  
  geom_vline(data = filter(data_stats, 
                           m == "United Kingdom"),
             aes(xintercept = log(mean_SD_prop)),
             col ="#FC8D62", linewidth = 0.8)+
  
  facet_wrap(~v,scale="free")+
  scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                "rw" = "#8DA0CB",
                                                "dlm" = "#E78AC3")) +
  theme(legend.position = "bottom") +
  labs(#title = "Density Plot of the Proportion Mean of the SD across Windows - United Kingdom",
       x = "Log Proportion Mean SD") 
dev.off()


pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Comparison/output/density_log_mean_SD_proportion-clusters.pdf",
    height = 4.5)
ggplot(window_summary %>%
         filter(window_length == 15,
                measure == "mean_SD_prop",
                m == "United Arab Emirates") %>%
         mutate(cluster = 'cluster 1 - UAE') %>%
         rbind(window_summary %>%
                 filter(window_length == 15,
                        measure == "mean_SD_prop",
                        m == "New Zealand") %>%
                 mutate(cluster = 'cluster 2 - NZ')) %>%
         rbind(window_summary %>%
         filter(window_length == 15,
                measure == "mean_SD_prop",
                m == "United Kingdom") %>%
         mutate(cluster = 'cluster 3 - UK')) %>%
         
         mutate(method = factor(method, levels = c("rw",
                                                   "dlm",
                                                   "hmm"))) %>%
         mutate(v = factor(v, levels = c("alpha",
                                         "beta",
                                         "gamma",
                                         "delta",
                                         "omicron"))) %>%
         
         mutate(value = if_else(value == 0, 0.001, value)),
       
       aes(x = log(value), fill = method)) +
  
  geom_density(alpha = 0.4)+
  
  geom_vline(data = filter(data_stats, 
                           m == "United Kingdom"),
             aes(xintercept = log(mean_SD_prop)),
             col ="#FC8D62", linewidth = 0.8)+
  
  facet_grid2(cluster~v, scales="free", independent='y')+
  scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                "rw" = "#8DA0CB",
                                                "dlm" = "#E78AC3")) +
  theme(legend.position = "bottom",
        legend.margin = margin(t = -15)) +
  labs(#title = "Density Plot of the Proportion Mean of the SD across Windows - United Kingdom",
    x = "",
    y = "") 
dev.off()



#### UPPER Q SD PROPORTION ####
pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Comparison/output/density_log_upper_q_SD_proportion_UAE.pdf", width = 8, height = 6)
ggplot(window_summary %>%
         filter(window_length == 15,
                measure == "SD_upper_q_prop",
                m == "United Arab Emirates") %>% 
         mutate(method = factor(method, levels = c("rw",
                                                   "dlm",
                                                   "hmm"))) %>%
         mutate(v = factor(v, levels = c("alpha",
                                         "beta",
                                         "gamma",
                                         "delta",
                                         "omicron"))) %>%
         
         mutate(value = if_else(value == 0, 0.001, value)),
       
       aes(x = log(value), fill = method)) +
  
  geom_density(alpha = 0.4)+
  
  geom_vline(data = filter(data_stats, 
                           m == "United Arab Emirates"),
             aes(xintercept = log(SD_upper_q_prop)),
             col ="#FC8D62", linewidth = 0.8)+
  
  facet_wrap(~v,scale="free")+
  scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                "rw" = "#8DA0CB",
                                                "dlm" = "#E78AC3")) +
  theme(legend.position = "bottom") +
  labs(#title = "Density Plot of the Proportion Upper Quartile of the SD across Windows - United Arab Emirates",
       x = "Log Proportion Upper Quartile SD") 
dev.off()

pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Comparison/output/density_log_upper_q_SD_proportion_NZ.pdf", width = 8, height = 6)
ggplot(window_summary %>%
         filter(window_length == 15,
                measure == "SD_upper_q_prop",
                m == "New Zealand") %>% 
         mutate(method = factor(method, levels = c("rw",
                                                   "dlm",
                                                   "hmm"))) %>%
         mutate(v = factor(v, levels = c("alpha",
                                         "beta",
                                         "gamma",
                                         "delta",
                                         "omicron"))) %>%
         
         mutate(value = if_else(value == 0, 0.001, value)),
       
       aes(x = log(value), fill = method)) +
  
  geom_density(alpha = 0.4)+
  
  geom_vline(data = filter(data_stats, 
                           m == "New Zealand"),
             aes(xintercept = log(SD_upper_q_prop)),
             col ="#FC8D62", linewidth = 0.8)+
  
  facet_wrap(~v,scale="free")+
  scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                "rw" = "#8DA0CB",
                                                "dlm" = "#E78AC3")) +
  theme(legend.position = "bottom") +
  labs(#title = "Density Plot of the Proportion Upper Quartile of the SD across Windows - New Zealand",
       x = "Log Proportion Upper Quartile SD")  
dev.off()

pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Comparison/output/density_log_upper_q_SD_proportion_UK.pdf", width = 8, height = 6)
ggplot(window_summary %>%
         filter(window_length == 15,
                measure == "SD_upper_q_prop",
                m == "United Kingdom") %>% 
         mutate(method = factor(method, levels = c("rw",
                                                   "dlm",
                                                   "hmm"))) %>%
         mutate(v = factor(v, levels = c("alpha",
                                         "beta",
                                         "gamma",
                                         "delta",
                                         "omicron"))) %>%
         
         mutate(value = if_else(value == 0, 0.001, value)),
       
       aes(x = log(value), fill = method)) +
  
  geom_density(alpha = 0.4)+
  
  geom_vline(data = filter(data_stats, 
                           m == "United Kingdom"),
             aes(xintercept = log(SD_upper_q_prop)),
             col ="#FC8D62", linewidth = 0.8)+
  
  facet_wrap(~v,scale="free")+
  scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                "rw" = "#8DA0CB",
                                                "dlm" = "#E78AC3")) +
  theme(legend.position = "bottom") +
  labs(#title = "Density Plot of the Proportion Upper Quartile of the SD across Windows - United Kingdom",
       x = "Log Proportion Upper Quartile SD") 
dev.off()


pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Comparison/output/density_log_upper_q_SD_proportion-clusters.pdf",
    height = 4.5)
ggplot(window_summary %>%
         filter(window_length == 15,
                measure == "SD_upper_q_prop",
                m == "United Arab Emirates") %>%
         mutate(cluster = 'cluster 1 - UAE') %>%
         rbind(window_summary %>%
                 filter(window_length == 15,
                        measure == "SD_upper_q_prop",
                        m == "New Zealand") %>%
                 mutate(cluster = 'cluster 2 - NZ')) %>%
         rbind(window_summary %>%
                 filter(window_length == 15,
                        measure == "SD_upper_q_prop",
                        m == "United Kingdom") %>%
                 mutate(cluster = 'cluster 3 - UK')) %>%
         
         mutate(method = factor(method, levels = c("rw",
                                                   "dlm",
                                                   "hmm"))) %>%
         mutate(v = factor(v, levels = c("alpha",
                                         "beta",
                                         "gamma",
                                         "delta",
                                         "omicron"))) %>%
         
         mutate(value = if_else(value == 0, 0.001, value)),
       
       aes(x = log(value), fill = method)) +
  
  geom_density(alpha = 0.4)+
  
  geom_vline(data = filter(data_stats, 
                           m == "United Kingdom"),
             aes(xintercept = log(SD_upper_q_prop)),
             col ="#FC8D62", linewidth = 0.8)+
  
  facet_grid2(cluster~v, scales="free", independent='y')+
  scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                "rw" = "#8DA0CB",
                                                "dlm" = "#E78AC3")) +
  theme(axis.text.x = element_text(angle = 45)) +
  theme(legend.position = "bottom",
        legend.margin = margin(t = -15)) +
  labs(
    x = "",
    y = "") 
dev.off()


### MAE ####

window_summary_mae <- window_summary %>%
  filter(window_length == 15,
         measure %in% c("mean_SD_prop", "SD_upper_q_prop")) %>%
  left_join(data_stats %>%
              select(v, m, mean_SD_prop, SD_upper_q_prop) %>%
              pivot_longer(., cols = c(mean_SD_prop, SD_upper_q_prop), names_to = "measure", values_to = "original_data_value"), 
            by = c("v", "m", "measure")) %>%
  
  group_by(v, m, method, measure, cluster) %>%
  summarise(mae = mean(abs(value - original_data_value))) %>%
  ungroup() %>%
  mutate(v = factor(v, levels = c("alpha",
                                  "beta",
                                  "gamma",
                                  "delta",
                                  "omicron"))) %>%
  mutate(m = case_when(m == 'United Arab Emirates' ~ 'UAE',
                       m == 'United Kingdom' ~ 'UK',
                       m == 'New Zealand' ~ 'NZ',
                       m == 'Dominican Republic' ~ 'DR',
                       m == 'South Africa' ~ 'SA',
                       TRUE ~ m))


#### MEAN SD PROPORTION ####
pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Comparison/output/mae_mean_SD_proportion_cluster_1.pdf", 
    height = 2.25)
ggplot(data = window_summary_mae %>%
         filter(cluster == 1,
                measure == "mean_SD_prop") %>%
         mutate(cluster = as.character('cluster 1')) %>%
         mutate(method = factor(method, levels = c("rw",
                                                   "dlm",
                                                   "hmm"))),
       
       aes(x = m,
           y = mae,
           fill = method)) +
  
  geom_col(position=position_dodge()) +
  
  geom_hline(data = window_summary_mae %>% 
               filter(cluster == 1,
                      measure == "mean_SD_prop") %>%
               group_by(v, method) %>% 
               summarise(mae = mean(mae)) %>%
               mutate(method = factor(method, levels = c("rw",
                                                         "dlm",
                                                         "hmm"))),
             
             aes(yintercept = mae, 
                 col = method,
                 linetype = method),
             size = 1) +
  
  facet_grid2(cluster~v, scales = "free", independent='y') +
  scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                "rw" = "#8DA0CB",
                                                "dlm" = "#E78AC3")) +
  scale_color_manual(name = "Method", values = c("hmm" = "#FC8D62",
                                                 "rw" = "#FFD92F",
                                                 "dlm" = "#E5C494")) +
  
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  #theme(legend.position = "bottom") +
  #labs(color  = "Method", fill = "Method", linetype = "Method") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 7)) +
  ylab("") +
  xlab("")
dev.off()

pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Comparison/output/mae_mean_SD_proportion_cluster_2.pdf", 
    height = 2)
ggplot(data = window_summary_mae %>%
         filter(cluster == 2,
                measure == "mean_SD_prop") %>%
         mutate(cluster = as.character('cluster 2')) %>%
         mutate(method = factor(method, levels = c("rw",
                                                   "dlm",
                                                   "hmm"))),
       
       aes(x = m,
           y = mae,
           fill = method)) +
  
  geom_col(position=position_dodge()) +
  
  geom_hline(data = window_summary_mae %>% 
               filter(cluster == 2,
                      measure == "mean_SD_prop") %>%
               group_by(v, method) %>% 
               summarise(mae = mean(mae)) %>%
               mutate(method = factor(method, levels = c("rw",
                                                         "dlm",
                                                         "hmm"))),
             
             aes(yintercept = mae, 
                 col = method,
                 linetype = method),
             size = 1) +
  
  facet_grid2(cluster~v, scales = "free", independent='y') +
  theme(strip.text.x = element_blank()) +
  scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                "rw" = "#8DA0CB",
                                                "dlm" = "#E78AC3")) +
  scale_color_manual(name = "Method", values = c("hmm" = "#FC8D62",
                                                 "rw" = "#FFD92F",
                                                 "dlm" = "#E5C494")) +
  
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  #theme(legend.position = "bottom") +
  #labs(color  = "Method", fill = "Method", linetype = "Method") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 7)) +
  ylab("") +
  xlab("")
dev.off()

pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Comparison/output/mae_mean_SD_proportion_cluster_3.pdf", 
    height = 2)
ggplot(data = window_summary_mae %>%
         filter(cluster == 3,
                measure == "mean_SD_prop") %>%
         mutate(cluster = as.character('cluster 3')) %>%
         mutate(method = factor(method, levels = c("rw",
                                                   "dlm",
                                                   "hmm"))),
       
       aes(x = m,
           y = mae,
           fill = method)) +
  
  geom_col(position=position_dodge()) +
  
  geom_hline(data = window_summary_mae %>% 
               filter(cluster == 3,
                      measure == "mean_SD_prop") %>%
               group_by(v, method) %>% 
               summarise(mae = mean(mae)) %>%
               mutate(method = factor(method, levels = c("rw",
                                                         "dlm",
                                                         "hmm"))),
             
             aes(yintercept = mae, 
                 col = method,
                 linetype = method),
             size = 1) +
  
  facet_grid2(cluster~v, scales = "free", independent='y') +
  theme(strip.text.x = element_blank()) +
  scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                "rw" = "#8DA0CB",
                                                "dlm" = "#E78AC3")) +
  scale_color_manual(name = "Method", values = c("hmm" = "#FC8D62",
                                                 "rw" = "#FFD92F",
                                                 "dlm" = "#E5C494")) +
  
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  theme(
    legend.position = "bottom",         
    legend.justification = "center",     
    legend.margin = margin(t = -15)      
  ) +
  labs(color  = "Method", fill = "Method", linetype = "Method") +
  theme(axis.text.x = element_text(size = 7)) +
  ylab("") +
  xlab("") 
dev.off()


#### UPPER Q SD PROPORTION ####
pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Comparison/output/mae_upper_q_SD_proportion_cluster_1.pdf", 
    height = 2.25)
ggplot(data = window_summary_mae %>%
         filter(cluster == 1,
                measure == "SD_upper_q_prop") %>%
         mutate(cluster = as.character('cluster 1')) %>%
         mutate(method = factor(method, levels = c("rw",
                                                   "dlm",
                                                   "hmm"))),
       
       aes(x = m,
           y = mae,
           fill = method)) +
  
  geom_col(position=position_dodge()) +
  
  geom_hline(data = window_summary_mae %>% 
               filter(cluster == 1,
                      measure == "SD_upper_q_prop") %>%
               group_by(v, method) %>% 
               summarise(mae = mean(mae)) %>%
               mutate(method = factor(method, levels = c("rw",
                                                         "dlm",
                                                         "hmm"))),
             
             aes(yintercept = mae, 
                 col = method,
                 linetype = method),
             size = 1) +
  
  facet_grid2(cluster~v, scales = "free", independent='y') +
  scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                "rw" = "#8DA0CB",
                                                "dlm" = "#E78AC3")) +
  scale_color_manual(name = "Method", values = c("hmm" = "#FC8D62",
                                                 "rw" = "#FFD92F",
                                                 "dlm" = "#E5C494")) +
  
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  #theme(legend.position = "bottom") +
  #labs(color  = "Method", fill = "Method", linetype = "Method") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 7)) +
  ylab("") +
  xlab("")
dev.off()

pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Comparison/output/mae_upper_q_SD_proportion_cluster_2.pdf", 
    height = 2)
ggplot(data = window_summary_mae %>%
         filter(cluster == 2,
                measure == "SD_upper_q_prop") %>%
         mutate(cluster = as.character('cluster 2')) %>%
         mutate(method = factor(method, levels = c("rw",
                                                   "dlm",
                                                   "hmm"))),
       
       aes(x = m,
           y = mae,
           fill = method)) +
  
  geom_col(position=position_dodge()) +
  
  geom_hline(data = window_summary_mae %>% 
               filter(cluster == 2,
                      measure == "SD_upper_q_prop") %>%
               group_by(v, method) %>% 
               summarise(mae = mean(mae)) %>%
               mutate(method = factor(method, levels = c("rw",
                                                         "dlm",
                                                         "hmm"))),
             
             aes(yintercept = mae, 
                 col = method,
                 linetype = method),
             size = 1) +
  
  facet_grid2(cluster~v, scales = "free", independent='y') +
  theme(strip.text.x = element_blank()) +
  scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                "rw" = "#8DA0CB",
                                                "dlm" = "#E78AC3")) +
  scale_color_manual(name = "Method", values = c("hmm" = "#FC8D62",
                                                 "rw" = "#FFD92F",
                                                 "dlm" = "#E5C494")) +
  
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  #theme(legend.position = "bottom") +
  #labs(color  = "Method", fill = "Method", linetype = "Method") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 7)) +
  ylab("") +
  xlab("")
dev.off()

pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Comparison/output/mae_upper_q_SD_proportion_cluster_3.pdf", 
    height = 2)
ggplot(data = window_summary_mae %>%
         filter(cluster == 3,
                measure == "SD_upper_q_prop") %>%
         mutate(cluster = as.character('cluster 3')) %>%
         mutate(method = factor(method, levels = c("rw",
                                                   "dlm",
                                                   "hmm"))),
       
       aes(x = m,
           y = mae,
           fill = method)) +
  
  geom_col(position=position_dodge()) +
  
  geom_hline(data = window_summary_mae %>% 
               filter(cluster == 3,
                      measure == "SD_upper_q_prop") %>%
               group_by(v, method) %>% 
               summarise(mae = mean(mae)) %>%
               mutate(method = factor(method, levels = c("rw",
                                                         "dlm",
                                                         "hmm"))),
             
             aes(yintercept = mae, 
                 col = method,
                 linetype = method),
             size = 1) +
  
  facet_grid2(cluster~v, scales = "free", independent='y') +
  theme(strip.text.x = element_blank()) +
  scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                "rw" = "#8DA0CB",
                                                "dlm" = "#E78AC3")) +
  scale_color_manual(name = "Method", values = c("hmm" = "#FC8D62",
                                                 "rw" = "#FFD92F",
                                                 "dlm" = "#E5C494")) +
  
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  theme(
    legend.position = "bottom",         
    legend.justification = "center",     
    legend.margin = margin(t = -15)      
  ) +
  labs(color  = "Method", fill = "Method", linetype = "Method") +
  theme(axis.text.x = element_text(size = 7)) +
  ylab("") +
  xlab("") 
dev.off()


### WINDOW COMARISON ####

data_stats_window <- data_proportion %>%
  
  group_by(v,m)%>%
  partition(cluster)%>%
  
  summarise(mean_SD_count_5    = mean(window_stats(x,L[1],N,sd)),
            SD_upper_q_count_5 = quantile(window_stats(x,L[1],N,sd), 0.75),
            SD_lower_q_count_5 = quantile(window_stats(x,L[1],N,sd), 0.25),
            
            mean_SD_prop_5    = mean(window_stats(proportion,L[1],N,sd)),
            SD_upper_q_prop_5 = quantile(window_stats(proportion,L[1],N,sd), 0.75),
            SD_lower_q_prop_5 = quantile(window_stats(proportion,L[1],N,sd), 0.25),
            
            mean_SD_count_10    = mean(window_stats(x,L[2],N,sd)),
            SD_upper_q_count_10 = quantile(window_stats(x,L[2],N,sd), 0.75),
            SD_lower_q_count_10 = quantile(window_stats(x,L[2],N,sd), 0.25),
            
            mean_SD_prop_10    = mean(window_stats(proportion,L[2],N,sd)),
            SD_upper_q_prop_10 = quantile(window_stats(proportion,L[2],N,sd), 0.75),
            SD_lower_q_prop_10 = quantile(window_stats(proportion,L[2],N,sd), 0.25),
            
            mean_SD_count_15    = mean(window_stats(x,L[3],N,sd)),
            SD_upper_q_count_15 = quantile(window_stats(x,L[3],N,sd), 0.75),
            SD_lower_q_count_15 = quantile(window_stats(x,L[3],N,sd), 0.25),
            
            mean_SD_prop_15    = mean(window_stats(proportion,L[3],N,sd)),
            SD_upper_q_prop_15 = quantile(window_stats(proportion,L[3],N,sd), 0.75),
            SD_lower_q_prop_15 = quantile(window_stats(proportion,L[3],N,sd), 0.25),
            
            mean_SD_count_20    = mean(window_stats(x,L[4],N,sd)),
            SD_upper_q_count_20 = quantile(window_stats(x,L[4],N,sd), 0.75),
            SD_lower_q_count_20 = quantile(window_stats(x,L[4],N,sd), 0.25),
            
            mean_SD_prop_20    = mean(window_stats(proportion,L[4],N,sd)),
            SD_upper_q_prop_20 = quantile(window_stats(proportion,L[4],N,sd), 0.75),
            SD_lower_q_prop_20 = quantile(window_stats(proportion,L[4],N,sd), 0.25))%>%
  collect()

window_summary_comparison <- window_summary %>%
  left_join(data_stats_window %>%
              pivot_longer(., cols = mean_SD_count_5:SD_lower_q_prop_20, names_to = "window_length", values_to = "original_data_value") %>%
              mutate(measure = window_length, .before = "original_data_value") %>%
              mutate(window_length = gsub(".*?(\\d+)", "\\1", window_length),
                     window_length = as.numeric(window_length),
                     measure = gsub("_\\d+", "", measure)), 
            by = c("v", "m", "window_length", "measure")) %>%
  
  group_by(v, m, method, measure, window_length, cluster) %>%
  summarise(mae = mean(abs(value - original_data_value))) %>%
  ungroup() %>%
  mutate(v = factor(v, levels = c("alpha",
                                  "beta",
                                  "gamma",
                                  "delta",
                                  "omicron")))


#### MEAN SD PROPORTION ####
pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Comparison/output/line_mae_mean_SD_proportion_cluster_1.pdf", width = 8, height = 6)
ggplot(data = window_summary_comparison %>%
         filter(measure == "mean_SD_prop",
                cluster == 1) %>%
         mutate(method = factor(method, levels = c("rw",
                                                   "dlm",
                                                   "hmm"))) %>%
         group_by(v, window_length, cluster, method) %>%
         summarise(mae = median(mae)) %>%
         ungroup(),
       
       aes(x = window_length,
           y = mae,
           col = method)) +
  geom_line(aes(linetype = method),
            size = 1) +
  facet_wrap(~v, scale = "free") +
  
  scale_color_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
  theme(legend.position = "bottom") +
  
  labs(color  = "Method", linetype = "Method") +
  ylab("MAE") +
  xlab("Window Length") 
dev.off()

pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Comparison/output/line_mae_mean_SD_proportion_cluster_2.pdf", width = 8, height = 6)
ggplot(data = window_summary_comparison %>%
         filter(measure == "mean_SD_prop",
                cluster == 2) %>%
         mutate(method = factor(method, levels = c("rw",
                                                   "dlm",
                                                   "hmm"))) %>%
         group_by(v, window_length, cluster, method) %>%
         summarise(mae = median(mae)) %>%
         ungroup(),
       
       aes(x = window_length,
           y = mae,
           col = method)) +
  geom_line(aes(linetype = method),
            size = 1) +
  facet_wrap(~v, scale = "free") +
  
  scale_color_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
  theme(legend.position = "bottom") +
  
  labs(color  = "Method", linetype = "Method") +
  ylab("MAE") +
  xlab("Window Length") 
dev.off()

pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Comparison/output/line_mae_mean_SD_proportion_cluster_3.pdf", width = 8, height = 6)
ggplot(data = window_summary_comparison %>%
         filter(measure == "mean_SD_prop",
                cluster == 3) %>%
         mutate(method = factor(method, levels = c("rw",
                                                   "dlm",
                                                   "hmm"))) %>%
         group_by(v, window_length, cluster, method) %>%
         summarise(mae = median(mae)) %>%
         ungroup(),
       
       aes(x = window_length,
           y = mae,
           col = method)) +
  geom_line(aes(linetype = method),
            size = 1) +
  facet_wrap(~v, scale = "free") +
  
  scale_color_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
  theme(legend.position = "bottom") +
  
  labs(color  = "Method", linetype = "Method") +
  ylab("MAE") +
  xlab("Window Length")
dev.off()


pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Comparison/output/line_mae_mean_SD_proportion_cluster_all.pdf", 
    height = 4)
ggplot(data = window_summary_comparison %>%
         filter(measure == "mean_SD_prop") %>%
         mutate(method = factor(method, levels = c("rw",
                                                   "dlm",
                                                   "hmm"))) %>%
         mutate(cluster = case_when(cluster == 1 ~ "cluster 1",
                                    cluster == 2 ~ "cluster 2",
                                    cluster == 3 ~ "cluster 3")) %>%
         group_by(v, window_length, cluster, method) %>%
         summarise(mae = median(mae)) %>%
         ungroup(),
       
       aes(x = window_length,
           y = mae,
           col = method)) +
  geom_line(aes(linetype = method),
            size = 1) +
  facet_grid(cluster~v, 
             scale = "free") +
  
  scale_color_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
  theme(legend.position = "bottom",
        legend.margin = margin(t = -15)) +
  
  labs(color  = "Method", linetype = "Method") +
  ylab("") +
  xlab("") 
dev.off()


#### UPPER Q SD PROPORTION ####
pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Comparison/output/line_mae_upper_q_SD_proportion_cluster_1.pdf", width = 8, height = 6)
ggplot(data = window_summary_comparison %>%
         filter(measure == "SD_upper_q_prop",
                cluster == 1) %>%
         mutate(method = factor(method, levels = c("rw",
                                                   "dlm",
                                                   "hmm"))) %>%
         group_by(v, window_length, cluster, method) %>%
         summarise(mae = median(mae)) %>%
         ungroup(),
       
       aes(x = window_length,
           y = mae,
           col = method)) +
  geom_line(aes(linetype = method),
            size = 1) +
  facet_wrap(~v, scale = "free") +
  
  scale_color_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
  theme(legend.position = "bottom") +
  
  labs(color  = "Method", linetype = "Method") +
  ylab("MAE") +
  xlab("Window Length") 
dev.off()

pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Comparison/output/line_mae_upper_q_SD_proportion_cluster_2.pdf", width = 8, height = 6)
ggplot(data = window_summary_comparison %>%
         filter(measure == "SD_upper_q_prop",
                cluster == 2) %>%
         mutate(method = factor(method, levels = c("rw",
                                                   "dlm",
                                                   "hmm"))) %>%
         group_by(v, window_length, cluster, method) %>%
         summarise(mae = median(mae)) %>%
         ungroup(),
       
       aes(x = window_length,
           y = mae,
           col = method)) +
  geom_line(aes(linetype = method),
            size = 1) +
  facet_wrap(~v, scale = "free") +
  
  scale_color_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
  theme(legend.position = "bottom") +
  
  labs(color  = "Method", linetype = "Method") +
  ylab("MAE") +
  xlab("Window Length") 
dev.off()

pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Comparison/output/line_mae_upper_q_SD_proportion_cluster_3.pdf", width = 8, height = 6)
ggplot(data = window_summary_comparison %>%
         filter(measure == "SD_upper_q_prop",
                cluster == 3) %>%
         mutate(method = factor(method, levels = c("rw",
                                                   "dlm",
                                                   "hmm"))) %>%
         group_by(v, window_length, cluster, method) %>%
         summarise(mae = median(mae)) %>%
         ungroup(),
       
       aes(x = window_length,
           y = mae,
           col = method)) +
  geom_line(aes(linetype = method),
            size = 1) +
  facet_wrap(~v, scale = "free") +
  
  scale_color_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
  theme(legend.position = "bottom") +
  
  labs(color  = "Method", linetype = "Method") +
  ylab("MAE") +
  xlab("Window Length") 
dev.off()



pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/Comparison/output/line_mae_upper_Q_SD_proportion_cluster_all.pdf", 
    height = 4)
ggplot(data = window_summary_comparison %>%
         filter(measure == "SD_upper_q_prop") %>%
         mutate(method = factor(method, levels = c("rw",
                                                   "dlm",
                                                   "hmm"))) %>%
         mutate(cluster = case_when(cluster == 1 ~ "cluster 1",
                                    cluster == 2 ~ "cluster 2",
                                    cluster == 3 ~ "cluster 3")) %>%
         group_by(v, window_length, cluster, method) %>%
         summarise(mae = median(mae)) %>%
         ungroup(),
       
       aes(x = window_length,
           y = mae,
           col = method)) +
  geom_line(aes(linetype = method),
            size = 1) +
  facet_grid(cluster~v, 
             scale = "free") +
  
  scale_color_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
  theme(legend.position = "bottom",
        legend.margin = margin(t = -15)) +
  
  labs(color  = "Method", linetype = "Method") +
  ylab("") +
  xlab("") 
dev.off()


#######################.


#######################.