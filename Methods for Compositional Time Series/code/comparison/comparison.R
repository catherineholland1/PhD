#################################################.
## WINDOW DENSITY COMPARISON - HMM / RW  / DLM ##
#################################################.


## PACKAGES ####

library(dplyr)
library(tidyr)
library(abind)
library(reshape2)
library(multidplyr)
library(ggplot2)


#######################.

## LOAD IN OUTPUT ####

## SET WORKING DIRECTORY ##

setwd("comparison")


## LOAD IN DATA ##

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



variant_data <- readRDS(file = "") %>%
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


#######################.

## SIMULATION OUTPUT ####


## LOAD IN HMM ##

load('HMM/output/replicates_GDM.RData')
hmm_simulation_replicates <- simulation_replicates

rm(simulation_replicates)


## LOAD IN RW ##

load('RW/output/replicates_RW.RData')
rw_simulation_replicates <- simulation_replicates

rm(simulation_replicates)


## LOAD IN DLM ##

load('replicates_DLM.RData')
dlm_simulation_replicates <- simulation_replicates

rm(simulation_replicates)



### NUMBER OF SAMPLES ####

N_samples <- c(hmm = dim(hmm_simulation_replicates)[4], 
               rw = dim(rw_simulation_replicates)[4],
               dlm = dim(dlm_simulation_replicates)[4])


#######################.

## WINDOWS ####

### SETUP ####

cluster <- new_cluster(4)

# WINDOW WIDTH
L <- c(5, 10, 15, 20, 30) 


cluster_assign(cluster, window_stats=function(x,L,N,stat_fun){
  w_starts <- 1:(N-L+1)
  values <- sapply(w_starts,function(u)stat_fun(x[u:(u+L-1)]))
  return(values)
},L=L,N=N)


#######################.

## SUMMARISE DATA ####

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
            
            mean_SD_count=mean(window_stats(x,L[2],N,sd)),
            SD_upper_q_count=quantile(window_stats(x,L[2],N,sd), 0.75),
            SD_lower_q_count=quantile(window_stats(x,L[2],N,sd), 0.25),
            
            mean_prop=mean(proportion),
            median_prop=median(proportion),
            
            var_prop=var(proportion),
            IQR_prop=IQR(proportion),
            
            mean_SD_prop=mean(window_stats(proportion,L[2],N,sd)),
            SD_upper_q_prop=quantile(window_stats(proportion,L[2],N,sd), 0.75),
            SD_lower_q_prop=quantile(window_stats(proportion,L[2],N,sd), 0.25),)%>%
  collect()


#######################.

## SUMMARISE HMM ####

hmm_replicates <- hmm_simulation_replicates %>%
  
  melt(varnames=c("t","v","m", "sample"), value.name = "x") %>%
  
  # ADD Y - TOTAL COUNT
  left_join(y_total_count_array %>% 
              select(-date) %>% 
              mutate(m_numeric = as.numeric(factor(m, levels = unique(m)))) %>%
              select(-m) %>%
              rename(m = "m_numeric"), by = c("t", "m")) %>%
  
  group_by(t, m)%>%
  partition(cluster)%>%
  
  # PROPORTION
  mutate(proportion = dplyr::if_else(y == 0, 0, x/y)) %>%
  
  ungroup() %>%
  collect()


hmm_replicates_stats <- hmm_replicates %>%
  
  group_by(v, m, sample)%>%
  partition(cluster)%>%
  
  summarise(total = sum(x),
            
            mean_count=mean(x),
            median_count=median(x),
            
            var_count=var(x),
            IQR_count=IQR(x),
            
            mean_SD_count=mean(window_stats(x,L[2],N,sd)),
            SD_upper_q_count=quantile(window_stats(x,L[2],N,sd), 0.75),
            SD_lower_q_count=quantile(window_stats(x,L[2],N,sd), 0.25),
            
            mean_prop=mean(proportion),
            median_prop=median(proportion),
            
            var_prop=var(proportion),
            IQR_prop=IQR(proportion),
            
            mean_SD_prop=mean(window_stats(proportion,L[2],N,sd)),
            SD_upper_q_prop=quantile(window_stats(proportion,L[2],N,sd), 0.75),
            SD_lower_q_prop=quantile(window_stats(proportion,L[2],N,sd), 0.25))%>%
  collect()


#######################.

## SUMMARISE RW ####

rw_replicates <- rw_simulation_replicates %>%
  
  melt(varnames=c("t","v","m", "sample"), value.name = "x") %>%
  
  # ADD Y - TOTAL COUNT
  left_join(y_total_count_array %>% 
              select(-date) %>% 
              mutate(m_numeric = as.numeric(factor(m, levels = unique(m)))) %>%
              select(-m) %>%
              rename(m = "m_numeric"), by = c("t", "m")) %>%
  
  group_by(t, m)%>%
  partition(cluster)%>%
  
  # PROPORTION
  mutate(proportion = dplyr::if_else(y == 0, 0, x/y)) %>%
  
  ungroup() %>%
  collect()


rw_replicates_stats <- rw_replicates %>%
  
  group_by(v, m, sample)%>%
  partition(cluster)%>%
  
  summarise(total = sum(x),
            
            mean_count=mean(x),
            median_count=median(x),
            
            var_count=var(x),
            IQR_count=IQR(x),
            
            mean_SD_count=mean(window_stats(x,L[2],N,sd)),
            SD_upper_q_count=quantile(window_stats(x,L[2],N,sd), 0.75),
            SD_lower_q_count=quantile(window_stats(x,L[2],N,sd), 0.25),
            
            mean_prop=mean(proportion),
            median_prop=median(proportion),
            
            var_prop=var(proportion),
            IQR_prop=IQR(proportion),
            
            mean_SD_prop=mean(window_stats(proportion,L[2],N,sd)),
            SD_upper_q_prop=quantile(window_stats(proportion,L[2],N,sd), 0.75),
            SD_lower_q_prop=quantile(window_stats(proportion,L[2],N,sd), 0.25))%>%
  collect()


#######################.

## SUMMARISE DLM ####

dlm_replicates <- dlm_simulation_replicates %>%
  
  melt(varnames=c("t","v","m", "sample"), value.name = "x") %>%
  
  # ADD Y - TOTAL COUNT
  left_join(y_total_count_array %>% 
              select(-date) %>% 
              mutate(m_numeric = as.numeric(factor(m, levels = unique(m)))) %>%
              select(-m) %>%
              rename(m = "m_numeric"), by = c("t", "m")) %>%
  
  group_by(t, m)%>%
  partition(cluster)%>%
  
  # PROPORTION
  mutate(proportion = dplyr::if_else(y == 0, 0, x/y)) %>%
  
  ungroup() %>%
  collect()


dlm_replicates_stats <- dlm_replicates %>%
  
  group_by(v, m, sample)%>%
  partition(cluster)%>%
  
  summarise(total = sum(x),
            
            mean_count=mean(x),
            median_count=median(x),
            
            var_count=var(x),
            IQR_count=IQR(x),
            
            mean_SD_count=mean(window_stats(x,L[1],N,sd)),
            SD_upper_q_count=quantile(window_stats(x,L[1],N,sd), 0.75),
            SD_lower_q_count=quantile(window_stats(x,L[1],N,sd), 0.25),
            
            mean_prop=mean(proportion),
            median_prop=median(proportion),
            
            var_prop=var(proportion),
            IQR_prop=IQR(proportion),
            
            mean_SD_prop=mean(window_stats(proportion,L[1],N,sd)),
            SD_upper_q_prop=quantile(window_stats(proportion,L[1],N,sd), 0.75),
            SD_lower_q_prop=quantile(window_stats(proportion,L[1],N,sd), 0.25))%>%
  collect()


#######################.

## COMPARE REPLICATES ####

replicates_compare <- rbind(mutate(hmm_replicates_stats, method = "hmm"),
                            mutate(rw_replicates_stats, method = "rw"),
                            mutate(dlm_replicates_stats, method = "dlm"))%>%
  mutate(v=variants[v], m=countries_cluster$countries[m])


### PLOT ####


#### MEAN SD - COUNT ####

pdf("output/mean_SD/countries_mean_SD_count_log.pdf",width=8,height=4.5)
for(i in 1:N_countries){
  plot(ggplot(replicates_compare %>%
                filter(m == countries_cluster$countries[i]) %>% 
                mutate(across(total:SD_lower_q_prop, ~if_else(.x == 0, 0.001, .x))),
              aes(x = log(mean_SD_count), fill = method))+
         
         geom_density(alpha=0.5)+
         
         geom_vline(data = filter(data_stats, 
                                  m == countries_cluster$countries[i]),
                    aes(xintercept = log(mean_SD_count)),
                    col ="#FC8D62", linewidth = 0.8)+
         
         facet_wrap(~v,scale="free")+
         scale_fill_manual(name = "", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
         labs(title=countries_cluster$countries[i])) 
}
dev.off()


#### MEAN SD - PROPORTION ####

pdf("output/mean_SD/countries_mean_SD_proportion_log.pdf",width=8,height=4.5)
for(i in 1:N_countries){
  plot(ggplot(replicates_compare %>%
                filter(m == countries_cluster$countries[i]) %>% 
                mutate(across(total:SD_lower_q_prop, ~if_else(.x == 0, 0.001, .x))),
              aes(x = log(mean_SD_prop), fill = method))+
         
         geom_density(alpha=0.5)+
         
         geom_vline(data = filter(data_stats, 
                                  m == countries_cluster$countries[i]),
                    aes(xintercept = log(mean_SD_prop)),
                    col ="#FC8D62", linewidth = 0.8)+
         
         facet_wrap(~v,scale="free")+
         scale_fill_manual(name = "", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
         labs(title=countries_cluster$countries[i])) 
}
dev.off()


#### SD UPPER Q - COUNT ####

pdf("output/SD_upper_q/countries_upper_q_SD_count_log.pdf",width=8,height=4.5)
for(i in 1:N_countries){
  plot(ggplot(replicates_compare %>%
                filter(m==countries_cluster$countries[i]) %>% 
                mutate(across(total:SD_lower_q_prop, ~if_else(.x == 0, 0.001, .x))),
              aes(x=log(SD_upper_q_count), fill=method))+
         
         geom_density(alpha=0.5)+
         
         geom_vline(data=filter(data_stats,
                                m==countries_cluster$countries[i]),
                    aes(xintercept = log(SD_upper_q_count)),
                    col ="#FC8D62", linewidth = 0.8)+
         
         facet_wrap(~v,scale="free")+
         scale_fill_manual(name = "", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
         labs(title=countries_cluster$countries[i])) 
}
dev.off()


#### SD UPPER Q - PROPORTION ####

pdf("output/SD_upper_q/countries_upper_q_SD_proportion_log.pdf",width=8,height=4.5)
for(i in 1:N_countries){
  plot(ggplot(replicates_compare %>%
                filter(m==countries_cluster$countries[i]) %>% 
                mutate(across(total:SD_lower_q_prop, ~if_else(.x == 0, 0.001, .x))),
              aes(x=log(SD_upper_q_prop), fill=method))+
         geom_density(alpha=0.5)+
         geom_vline(data=filter(data_stats,
                                m==countries_cluster$countries[i]),
                    aes(xintercept = log(SD_upper_q_prop)),
                    col ="#FC8D62", linewidth = 0.8)+
         facet_wrap(~v,scale="free")+
         scale_fill_manual(name = "", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
         labs(title=countries_cluster$countries[i])) 
}
dev.off()


#### SD LOWER Q - COUNT ####

pdf("output/SD_lower_q/countries_lower_q_SD_count_log.pdf",width=8,height=4.5)
for(i in 1:N_countries){
  plot(ggplot(replicates_compare %>%
                filter(m==countries_cluster$countries[i]) %>% 
                mutate(across(total:SD_lower_q_prop, ~if_else(.x == 0, 0.001, .x))),
              aes(x=log(SD_lower_q_count), fill=method))+
         
         geom_density(alpha=0.5)+
         
         geom_vline(data=filter(data_stats,
                                m==countries_cluster$countries[i]),
                    aes(xintercept = log(SD_lower_q_count)),
                    col ="#FC8D62", linewidth = 0.8)+
         
         facet_wrap(~v,scale="free")+
         scale_fill_manual(name = "", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
         labs(title=countries_cluster$countries[i])) 
}
dev.off()


#### SD LOWER Q - PROPORTION ####

pdf("output/SD_lower_q/countries_lower_q_SD_proportion_log.pdf",width=8,height=4.5)
for(i in 1:N_countries){
  plot(ggplot(replicates_compare %>%
                filter(m==countries_cluster$countries[i]) %>% 
                mutate(across(total:SD_lower_q_prop, ~if_else(.x == 0, 0.001, .x))),
              aes(x=log(SD_lower_q_prop), fill=method))+
         geom_density(alpha=0.5)+
         geom_vline(data=filter(data_stats,
                                m==countries_cluster$countries[i]),
                    aes(xintercept = log(SD_lower_q_prop)),
                    col ="#FC8D62", linewidth = 0.8)+
         facet_wrap(~v,scale="free")+
         scale_fill_manual(name = "", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
         labs(title=countries_cluster$countries[i])) 
}
dev.off()


#### TOTAL ####

pdf("output/total/countries_total_cases_log.pdf",width=8,height=4.5)
for(i in 1:N_countries){
  plot(ggplot(replicates_compare %>%
                filter(m==countries_cluster$countries[i]) %>% 
                mutate(across(total:SD_lower_q_prop, ~if_else(.x == 0, 0.001, .x))),
              aes(x=log(total), fill=method))+
         
         geom_density(alpha=0.5)+
         
         geom_vline(data=filter(data_stats,
                                m==countries_cluster$countries[i]),
                    aes(xintercept = log(total)),
                    col ="#FC8D62", linewidth = 0.8)+
         
         facet_wrap(~v,scale="free")+
         scale_fill_manual(name = "", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
         labs(title=countries_cluster$countries[i])) 
}
dev.off()



### PLOT - BY CLUSTER ####

replicates_compare_cluster <- replicates_compare %>%
  left_join(countries_cluster %>%
              rename(m = "countries"), by = "m") %>%
  mutate(cluster = as.numeric(cluster))


countries_cluster_order <- countries_cluster %>%
  arrange(cluster)


#### MEAN SD - PROPORTION ####

##### CLUSTER 1 ####

cluster_value = 1

pdf("output/mean_SD/countries_mean_SD_proportion_log-cluster1.pdf",width=8,height=4.5)
for(i in 1:(N_countries/N_clusters)){
  plot(ggplot(replicates_compare_cluster %>%
                filter(cluster == cluster_value) %>%
                filter(m == (countries_cluster_order %>% 
                               filter(cluster == cluster_value) %>% 
                               select(countries) %>% 
                               unlist())[i]) %>% 
                mutate(across(total:SD_lower_q_prop, ~if_else(.x == 0, 0.001, .x))),
              
              aes(x = log(mean_SD_prop), fill = method))+
         
         geom_density(alpha=0.5)+
         
         geom_vline(data = filter(data_stats, 
                                  m == (countries_cluster_order %>% 
                                          filter(cluster == cluster_value) %>% 
                                          select(countries) %>% 
                                          unlist())[i]),
                    aes(xintercept = log(mean_SD_prop)),
                    col ="#FC8D62", linewidth = 0.8)+
         
         facet_wrap(~v,scale="free")+
         scale_fill_manual(name = "", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
         labs(title=paste0((countries_cluster_order %>% 
                              filter(cluster == cluster_value) %>% 
                              select(countries) %>% 
                              unlist())[i], " - Cluster ", cluster_value))
  )
}
dev.off()


##### CLUSTER 2 ####

cluster_value = 2

pdf("output/mean_SD/countries_mean_SD_proportion_log-cluster2.pdf",width=8,height=4.5)
for(i in 1:(N_countries/N_clusters)){
  plot(ggplot(replicates_compare_cluster %>%
                filter(cluster == cluster_value) %>%
                filter(m == (countries_cluster_order %>% 
                               filter(cluster == cluster_value) %>% 
                               select(countries) %>% 
                               unlist())[i]) %>% 
                mutate(across(total:SD_lower_q_prop, ~if_else(.x == 0, 0.001, .x))),
              
              aes(x = log(mean_SD_prop), fill = method))+
         
         geom_density(alpha=0.5)+
         
         geom_vline(data = filter(data_stats, 
                                  m == (countries_cluster_order %>% 
                                          filter(cluster == cluster_value) %>% 
                                          select(countries) %>% 
                                          unlist())[i]),
                    aes(xintercept = log(mean_SD_prop)),
                    col ="#FC8D62", linewidth = 0.8)+
         
         facet_wrap(~v,scale="free")+
         scale_fill_manual(name = "", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
         labs(title=paste0((countries_cluster_order %>% 
                              filter(cluster == cluster_value) %>% 
                              select(countries) %>% 
                              unlist())[i], " - Cluster ", cluster_value))
  )
}
dev.off()


##### CLUSTER 3 ####

cluster_value = 3

pdf("output/mean_SD/countries_mean_SD_proportion_log-cluster3.pdf",width=8,height=4.5)
for(i in 1:(N_countries/N_clusters)){
  plot(ggplot(replicates_compare_cluster %>%
                filter(cluster == cluster_value) %>%
                filter(m == (countries_cluster_order %>% 
                               filter(cluster == cluster_value) %>% 
                               select(countries) %>% 
                               unlist())[i]) %>% 
                mutate(across(total:SD_lower_q_prop, ~if_else(.x == 0, 0.001, .x))),
              
              aes(x = log(mean_SD_prop), fill = method))+
         
         geom_density(alpha=0.5)+
         
         geom_vline(data = filter(data_stats, 
                                  m == (countries_cluster_order %>% 
                                          filter(cluster == cluster_value) %>% 
                                          select(countries) %>% 
                                          unlist())[i]),
                    aes(xintercept = log(mean_SD_prop)),
                    col ="#FC8D62", linewidth = 0.8)+
         
         facet_wrap(~v,scale="free")+
         scale_fill_manual(name = "", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
         labs(title=paste0((countries_cluster_order %>% 
                              filter(cluster == cluster_value) %>% 
                              select(countries) %>% 
                              unlist())[i], " - Cluster ", cluster_value))
  )
}
dev.off()



#### SD UPPER Q - PROPORTION ####

##### CLUSTER 1 ####

cluster_value = 1

pdf("output/SD_upper_q/countries_upper_q_SD_proportion_log-cluster1.pdf",width=8,height=4.5)
for(i in 1:(N_countries/N_clusters)){
  plot(ggplot(replicates_compare_cluster %>%
                filter(cluster == cluster_value) %>%
                filter(m == (countries_cluster_order %>% 
                               filter(cluster == cluster_value) %>% 
                               select(countries) %>% 
                               unlist())[i]) %>% 
                mutate(across(total:SD_lower_q_prop, ~if_else(.x == 0, 0.001, .x))),
              
              aes(x = log(SD_upper_q_prop), fill = method))+
         
         geom_density(alpha=0.5)+
         
         geom_vline(data = filter(data_stats, 
                                  m == (countries_cluster_order %>% 
                                          filter(cluster == cluster_value) %>% 
                                          select(countries) %>% 
                                          unlist())[i]),
                    aes(xintercept = log(SD_upper_q_prop)),
                    col ="#FC8D62", linewidth = 0.8)+
         
         facet_wrap(~v,scale="free")+
         scale_fill_manual(name = "", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
         labs(title=paste0((countries_cluster_order %>% 
                              filter(cluster == cluster_value) %>% 
                              select(countries) %>% 
                              unlist())[i], " - Cluster ", cluster_value))
  )
}
dev.off()


##### CLUSTER 2 ####

cluster_value = 2

pdf("output/SD_upper_q/countries_upper_q_SD_proportion_log-cluster2.pdf",width=8,height=4.5)
for(i in 1:(N_countries/N_clusters)){
  plot(ggplot(replicates_compare_cluster %>%
                filter(cluster == cluster_value) %>%
                filter(m == (countries_cluster_order %>% 
                               filter(cluster == cluster_value) %>% 
                               select(countries) %>% 
                               unlist())[i]) %>% 
                mutate(across(total:SD_lower_q_prop, ~if_else(.x == 0, 0.001, .x))),
              
              aes(x = log(SD_upper_q_prop), fill = method))+
         
         geom_density(alpha=0.5)+
         
         geom_vline(data = filter(data_stats, 
                                  m == (countries_cluster_order %>% 
                                          filter(cluster == cluster_value) %>% 
                                          select(countries) %>% 
                                          unlist())[i]),
                    aes(xintercept = log(SD_upper_q_prop)),
                    col ="#FC8D62", linewidth = 0.8)+
         
         facet_wrap(~v,scale="free")+
         scale_fill_manual(name = "", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
         labs(title=paste0((countries_cluster_order %>% 
                              filter(cluster == cluster_value) %>% 
                              select(countries) %>% 
                              unlist())[i], " - Cluster ", cluster_value))
  )
}
dev.off()


##### CLUSTER 3 ####

cluster_value = 3

pdf("output/SD_upper_q/countries_upper_q_SD_proportion_log-cluster3.pdf",width=8,height=4.5)
for(i in 1:(N_countries/N_clusters)){
  plot(ggplot(replicates_compare_cluster %>%
                filter(cluster == cluster_value) %>%
                filter(m == (countries_cluster_order %>% 
                               filter(cluster == cluster_value) %>% 
                               select(countries) %>% 
                               unlist())[i]) %>% 
                mutate(across(total:SD_lower_q_prop, ~if_else(.x == 0, 0.001, .x))),
              
              aes(x = log(SD_upper_q_prop), fill = method))+
         
         geom_density(alpha=0.5)+
         
         geom_vline(data = filter(data_stats, 
                                  m == (countries_cluster_order %>% 
                                          filter(cluster == cluster_value) %>% 
                                          select(countries) %>% 
                                          unlist())[i]),
                    aes(xintercept = log(SD_upper_q_prop)),
                    col ="#FC8D62", linewidth = 0.8)+
         
         facet_wrap(~v,scale="free")+
         scale_fill_manual(name = "", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
         labs(title=paste0((countries_cluster_order %>% 
                              filter(cluster == cluster_value) %>% 
                              select(countries) %>% 
                              unlist())[i], " - Cluster ", cluster_value))
  )
}
dev.off()



#### SD LOWER Q - PROPORTION ####

##### CLUSTER 1 ####

cluster_value = 1

pdf("output/SD_lower_q/countries_lower_q_SD_proportion_log-cluster1.pdf",width=8,height=4.5)
for(i in 1:(N_countries/N_clusters)){
  plot(ggplot(replicates_compare_cluster %>%
                filter(cluster == cluster_value) %>%
                filter(m == (countries_cluster_order %>% 
                               filter(cluster == cluster_value) %>% 
                               select(countries) %>% 
                               unlist())[i]) %>% 
                mutate(across(total:SD_lower_q_prop, ~if_else(.x == 0, 0.001, .x))),
              
              aes(x = log(SD_lower_q_prop), fill = method))+
         
         geom_density(alpha=0.5)+
         
         geom_vline(data = filter(data_stats, 
                                  m == (countries_cluster_order %>% 
                                          filter(cluster == cluster_value) %>% 
                                          select(countries) %>% 
                                          unlist())[i]),
                    aes(xintercept = log(SD_lower_q_prop)),
                    col ="#FC8D62", linewidth = 0.8)+
         
         facet_wrap(~v,scale="free")+
         scale_fill_manual(name = "", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
         labs(title=paste0((countries_cluster_order %>% 
                              filter(cluster == cluster_value) %>% 
                              select(countries) %>% 
                              unlist())[i], " - Cluster ", cluster_value))
  )
}
dev.off()


##### CLUSTER 2 ####

cluster_value = 2

pdf("output/SD_lower_q/countries_lower_q_SD_proportion_log-cluster2.pdf",width=8,height=4.5)
for(i in 1:(N_countries/N_clusters)){
  plot(ggplot(replicates_compare_cluster %>%
                filter(cluster == cluster_value) %>%
                filter(m == (countries_cluster_order %>% 
                               filter(cluster == cluster_value) %>% 
                               select(countries) %>% 
                               unlist())[i]) %>% 
                mutate(across(total:SD_lower_q_prop, ~if_else(.x == 0, 0.001, .x))),
              
              aes(x = log(SD_lower_q_prop), fill = method))+
         
         geom_density(alpha=0.5)+
         
         geom_vline(data = filter(data_stats, 
                                  m == (countries_cluster_order %>% 
                                          filter(cluster == cluster_value) %>% 
                                          select(countries) %>% 
                                          unlist())[i]),
                    aes(xintercept = log(SD_lower_q_prop)),
                    col ="#FC8D62", linewidth = 0.8)+
         
         facet_wrap(~v,scale="free")+
         scale_fill_manual(name = "", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
         labs(title=paste0((countries_cluster_order %>% 
                              filter(cluster == cluster_value) %>% 
                              select(countries) %>% 
                              unlist())[i], " - Cluster ", cluster_value))
  )
}
dev.off()


##### CLUSTER 3 ####

cluster_value = 3

pdf("output/SD_lower_q/countries_lower_q_SD_proportion_log-cluster3.pdf",width=8,height=4.5)
for(i in 1:(N_countries/N_clusters)){
  plot(ggplot(replicates_compare_cluster %>%
                filter(cluster == cluster_value) %>%
                filter(m == (countries_cluster_order %>% 
                               filter(cluster == cluster_value) %>% 
                               select(countries) %>% 
                               unlist())[i]) %>% 
                mutate(across(total:SD_lower_q_prop, ~if_else(.x == 0, 0.001, .x))),
              
              aes(x = log(SD_lower_q_prop), fill = method))+
         
         geom_density(alpha=0.5)+
         
         geom_vline(data = filter(data_stats, 
                                  m == (countries_cluster_order %>% 
                                          filter(cluster == cluster_value) %>% 
                                          select(countries) %>% 
                                          unlist())[i]),
                    aes(xintercept = log(SD_lower_q_prop)),
                    col ="#FC8D62", linewidth = 0.8)+
         
         facet_wrap(~v,scale="free")+
         scale_fill_manual(name = "", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
         labs(title=paste0((countries_cluster_order %>% 
                              filter(cluster == cluster_value) %>% 
                              select(countries) %>% 
                              unlist())[i], " - Cluster ", cluster_value))
  )
}
dev.off()


#######################.

## MEAN ABSOLUTE ERROR (MAE) ####

mae <- replicates_compare %>%
  left_join(data_stats %>% 
              as.data.frame(.) %>% 
              rename_with(~ paste("data", .x, sep = "_"), total:SD_lower_q_prop),
            by = c("v", "m")) %>%
  
  group_by(v, m, method) %>%
  summarise(
    mae_mean_SD_count = mean(abs(mean_SD_count - data_mean_SD_count)),
    mae_mean_SD_prop = mean(abs(mean_SD_prop - data_mean_SD_prop)),
    mae_SD_upper_q_count = mean(abs(SD_upper_q_count - data_SD_upper_q_count)),
    mae_SD_upper_q_prop = mean(abs(SD_upper_q_prop - data_SD_upper_q_prop)),
    mae_SD_lower_q_count = mean(abs(SD_lower_q_count - data_SD_lower_q_count)),
    mae_SD_lower_q_prop = mean(abs(SD_lower_q_prop - data_SD_lower_q_prop))) %>%
  ungroup() %>%
  
  left_join(countries_cluster %>%
              rename(m = "countries"), by = "m") %>%
  mutate(cluster = as.numeric(cluster))


### PLOT - BY CLUSTER ####

#### MEAN SD - COUNT ####

pdf("output/mean_SD/mae_mean_SD_count-cluster.pdf",width=12,height=6)
for (i in 1:N_clusters) {
  plot(ggplot(data = mae %>%
                filter(cluster == i),
              
              aes(x = m,
                  y = mae_mean_SD_count,
                  fill = method)) +
         geom_col(position=position_dodge()) +
         
         geom_hline(data = mae %>% 
                      filter(cluster == i) %>%
                      group_by(v, method) %>% 
                      summarise(across(starts_with("mae"), mean)),
                    
                    aes(yintercept = mae_mean_SD_count, 
                        col = method,
                        linetype = method),
                    size = 1) +
         
         facet_wrap(~v, scale = "free") +
         scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                       "rw" = "#8DA0CB",
                                                       "dlm" = "#E78AC3")) +
         scale_color_manual(name = "Method", values = c("hmm" = "#FC8D62",
                                                        "rw" = "#FFD92F",
                                                        "dlm" = "#E5C494")) +
         
         theme(axis.text.x=element_text(angle=60, hjust=1)) +
         labs(color  = "Method", fill = "Method", linetype = "Method") +
         ylab("MAE") +
         xlab("Country") +
         ggtitle(paste0("MAE - Mean SD Count - Cluster ", i))
  )
}
dev.off()


#### MEAN SD - PROPORTION ####

pdf("output/mean_SD/mae_mean_SD_proportion-cluster.pdf",width=12,height=6)
for (i in 1:N_clusters) {
  plot(
    ggplot(data = mae %>%
             filter(cluster == i),
           
           aes(x = m,
               y = mae_mean_SD_prop,
               fill = method)) +
      geom_col(position=position_dodge()) +
      geom_hline(data = mae %>% 
                   filter(cluster == i) %>%
                   group_by(v, method) %>% 
                   summarise(across(starts_with("mae"), mean)),
                 aes(yintercept = mae_mean_SD_prop, 
                     col = method,
                     linetype = method),
                 size = 1) +
      facet_wrap(~v, scale = "free") +
      scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                    "rw" = "#8DA0CB",
                                                    "dlm" = "#E78AC3")) +
      scale_color_manual(name = "Method", values = c("hmm" = "#FC8D62",
                                                     "rw" = "#FFD92F",
                                                     "dlm" = "#E5C494")) +
      
      theme(axis.text.x=element_text(angle=60, hjust=1)) +
      labs(color  = "Method", fill = "Method", linetype = "Method") +
      ylab("MAE") +
      xlab("Country") +
      ggtitle(paste0("MAE - Mean SD Proportion - Cluster ", i))
  )
}
dev.off()


#### SD UPPER Q - COUNT ####

pdf("output/SD_upper_q/mae_upper_q_SD_count-cluster.pdf",width=12,height=6)
for (i in 1:N_clusters) {
  plot(
    ggplot(data = mae %>%
             filter(cluster == i),
           aes(x = m,
               y = mae_SD_upper_q_count,
               fill = method)) +
      geom_col(position=position_dodge()) +
      geom_hline(data = mae %>% 
                   filter(cluster == i) %>%
                   group_by(v, method) %>% 
                   summarise(across(starts_with("mae"), mean)),
                 aes(yintercept = mae_SD_upper_q_count, 
                     col = method,
                     linetype = method),
                 size = 1) +
      facet_wrap(~v, scale = "free") +
      scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                    "rw" = "#8DA0CB",
                                                    "dlm" = "#E78AC3")) +
      scale_color_manual(name = "Method", values = c("hmm" = "#FC8D62",
                                                     "rw" = "#FFD92F",
                                                     "dlm" = "#E5C494")) +
      
      theme(axis.text.x=element_text(angle=60, hjust=1)) +
      labs(color  = "Method", fill = "Method", linetype = "Method") +
      ylab("MAE") +
      xlab("Country") +
      ggtitle(paste0("MAE - SD Upper Q Count - Cluster ", i))
  )
}
dev.off()


#### SD UPPER Q - PROPORTION ####

pdf("output/SD_upper_q/mae_upper_q_SD_proportion-cluster.pdf",width=12,height=6)
for (i in 1:N_clusters) {
  plot(
    ggplot(data = mae %>%
             filter(cluster == i),
           aes(x = m,
               y = mae_SD_upper_q_prop,
               fill = method)) +
      
      geom_col(position=position_dodge()) +
      
      geom_hline(data = mae %>% 
                   filter(cluster == i) %>%
                   group_by(v, method) %>% 
                   summarise(across(starts_with("mae"), mean)),
                 aes(yintercept = mae_SD_upper_q_prop, 
                     col = method,
                     linetype = method),
                 size = 1) +
      
      facet_wrap(~v, scale = "free") +
      scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                    "rw" = "#8DA0CB",
                                                    "dlm" = "#E78AC3")) +
      scale_color_manual(name = "Method", values = c("hmm" = "#FC8D62",
                                                     "rw" = "#FFD92F",
                                                     "dlm" = "#E5C494")) +
      
      theme(axis.text.x=element_text(angle=60, hjust=1)) +
      labs(color  = "Method", fill = "Method", linetype = "Method") +
      ylab("MAE") +
      xlab("Country") +
      ggtitle(paste0("MAE - SD Upper Q Proportion - Cluster ", i))
  )
}
dev.off()


#### SD LOWER Q - COUNT ####

pdf("output/SD_lower_q/mae_lower_q_SD_count-cluster.pdf",width=12,height=6)
for (i in 1:N_clusters) {
  plot(
    ggplot(data = mae %>%
             filter(cluster == i),
           aes(x = m,
               y = mae_SD_lower_q_count,
               fill = method)) +
      
      geom_col(position=position_dodge()) +
      
      geom_hline(data = mae %>% 
                   filter(cluster == i) %>%
                   group_by(v, method) %>% 
                   summarise(across(starts_with("mae"), mean)),
                 aes(yintercept = mae_SD_lower_q_count, 
                     col = method,
                     linetype = method),
                 size = 1) +
      
      facet_wrap(~v, scale = "free") +
      scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                    "rw" = "#8DA0CB",
                                                    "dlm" = "#E78AC3")) +
      scale_color_manual(name = "Method", values = c("hmm" = "#FC8D62",
                                                     "rw" = "#FFD92F",
                                                     "dlm" = "#E5C494")) +
      
      theme(axis.text.x=element_text(angle=60, hjust=1)) +
      labs(color  = "Method", fill = "Method", linetype = "Method") +
      ylab("MAE") +
      xlab("Country") +
      ggtitle(paste0("MAE - SD Lower Q Count - Cluster ", i))
  )
}
dev.off()


#### SD LOWER Q - PROPORTION ####

pdf("output/SD_lower_q/mae_lower_q_SD_proportion-cluster.pdf",width=12,height=6)
for (i in 1:N_clusters) {
  plot(
    ggplot(data = mae %>%
             filter(cluster == i),
           aes(x = m,
               y = mae_SD_lower_q_prop,
               fill = method)) +
      
      geom_col(position=position_dodge()) +
      
      geom_hline(data = mae %>% 
                   filter(cluster == i) %>%
                   group_by(v, method) %>% 
                   summarise(across(starts_with("mae"), mean)),
                 aes(yintercept = mae_SD_lower_q_prop, 
                     col = method,
                     linetype = method),
                 size = 1) +
      
      facet_wrap(~v, scale = "free") +
      scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                    "rw" = "#8DA0CB",
                                                    "dlm" = "#E78AC3")) +
      scale_color_manual(name = "Method", values = c("hmm" = "#FC8D62",
                                                     "rw" = "#FFD92F",
                                                     "dlm" = "#E5C494")) +
      
      theme(axis.text.x=element_text(angle=60, hjust=1)) +
      labs(color  = "Method", fill = "Method", linetype = "Method") +
      ylab("MAE") +
      xlab("Country") +
      ggtitle(paste0("MAE - SD Lower Q Proportion - Cluster ", i))
  )
}
dev.off()


#######################.

## TAIL PROBABILITIES ####

tail_data <- replicates_compare_cluster %>%
  left_join(data_stats %>% 
              rename(data_mean_SD_count = "mean_SD_count",
                     data_SD_upper_q_count = "SD_upper_q_count",
                     data_SD_lower_q_count = "SD_lower_q_count",
                     data_mean_SD_prop = "mean_SD_prop",
                     data_SD_upper_q_prop = "SD_upper_q_prop",
                     data_SD_lower_q_prop = "SD_lower_q_prop") %>%
              select(v, m, starts_with("data")), 
            by = c("v", "m")) 


tail_prob <- tail_data %>%
  
  group_by(v, m, method) %>%
  
  summarise(proportion_mean_SD_count = (sum(mean_SD_count <= data_mean_SD_count)),
            proportion_mean_SD_prop = (sum(mean_SD_prop <= data_mean_SD_prop)),
            proportion_SD_upper_q_count = (sum(SD_upper_q_count <= data_SD_upper_q_count)),
            proportion_SD_upper_q_prop = (sum(SD_upper_q_prop <= data_SD_upper_q_prop)),
            proportion_SD_lower_q_count = (sum(SD_lower_q_count <= data_SD_lower_q_count)),
            proportion_SD_lower_q_prop = (sum(SD_lower_q_prop <= data_SD_lower_q_prop))) %>%
  
  mutate(N_samples = if_else(method == "hmm", N_samples[1], N_samples[2])) %>%
  
  group_by(v, m, method) %>%
  summarise(tail_prob_mean_SD_count = (proportion_mean_SD_count / N_samples) * 100,
            tail_prob_mean_SD_prop = (proportion_mean_SD_prop / N_samples) * 100,
            tail_prob_SD_upper_q_count = (proportion_SD_upper_q_count / N_samples) * 100,
            tail_prob_SD_upper_q_prop = (proportion_SD_upper_q_prop / N_samples) * 100,
            tail_prob_SD_lower_q_count = (proportion_SD_lower_q_count / N_samples) * 100,
            tail_prob_SD_lower_q_prop = (proportion_SD_lower_q_prop / N_samples) * 100) %>%
  
  ungroup() %>%
  mutate(across(starts_with("tail"), ~ round(., 2)))


### EXTREME TAIL PROBABILITES ####

extreme_tail_prob <- tail_prob %>%
  pivot_longer(., cols = starts_with("tail"),
               names_to = "tail_prob",
               values_to = "value") %>%
  #filter(value > 97.5 | value < 2.5)
  filter(value > 97.5)



### PLOT ####

plot_tail_prob <- replicates_compare_cluster %>%
  left_join(tail_prob, by = c("v", "m", "method"))


#### MEAN SD - COUNT ####

pdf("output/mean_SD/countries_mean_SD_count_log-tail.pdf",width=8,height=4.5)
for(i in 1:N_countries){
  plot(ggplot(plot_tail_prob %>%
                filter(m == countries_cluster$countries[i]) %>% 
                mutate(across(total:SD_lower_q_prop, ~if_else(.x == 0, 0.001, .x))),
              aes(x = log(mean_SD_count), fill = method))+
         
         geom_density(alpha=0.5)+
         
         geom_vline(data = filter(data_stats, 
                                  m == countries_cluster$countries[i]),
                    aes(xintercept = log(mean_SD_count)),
                    col ="#FC8D62", linewidth = 0.8)+
         
         geom_text(data = tail_prob %>%
                     filter(m == countries_cluster$countries[i],
                            method == "hmm"),
                   aes(label = paste0("hmm: ", tail_prob_mean_SD_count, "%"),
                       x = Inf, y = Inf),
                   hjust = 1, vjust = 2,
                   size = 3,
                   fontface = "bold",
                   col = "#66C2A5") +
         
         geom_text(data = tail_prob %>%
                     filter(m == countries_cluster$countries[i],
                            method == "rw"),
                   aes(label = paste0("rw: ", tail_prob_mean_SD_count, "%"),
                       x = Inf, y = Inf),
                   hjust = 1, vjust = 4,
                   size = 3,
                   fontface = "bold",
                   col = "#8DA0CB") +
         
         geom_text(data = tail_prob %>%
                     filter(m == countries_cluster$countries[i],
                            method == "dlm"),
                   aes(label = paste0("dlm: ", tail_prob_mean_SD_count, "%"),
                       x = Inf, y = Inf),
                   hjust = 1, vjust = 6,
                   size = 3,
                   fontface = "bold",
                   col = "#E78AC3") +
         
         facet_wrap(~v,scale="free")+
         scale_fill_manual(name = "", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
         labs(title=paste0(countries_cluster$countries[i], " - Cluster ", countries_cluster$cluster[i]))
  )
}
dev.off()


#### MEAN SD - PROPORTION ####

pdf("output/mean_SD/countries_mean_SD_proportion_log-tail.pdf",width=8,height=4.5)
for(i in 1:N_countries){
  plot(ggplot(plot_tail_prob %>%
                filter(m == countries_cluster$countries[i]) %>% 
                mutate(across(total:SD_lower_q_prop, ~if_else(.x == 0, 0.001, .x))),
              aes(x = log(mean_SD_prop), fill = method))+
         
         geom_density(alpha=0.5)+
         
         geom_vline(data = filter(data_stats, 
                                  m == countries_cluster$countries[i]),
                    aes(xintercept = log(mean_SD_prop)),
                    col ="#FC8D62", linewidth = 0.8)+
         
         geom_text(data = tail_prob %>%
                     filter(m == countries_cluster$countries[i],
                            method == "hmm"),
                   aes(label = paste0("hmm: ", tail_prob_mean_SD_prop, "%"),
                       x = Inf, y = Inf),
                   hjust = 1, vjust = 2,
                   size = 3,
                   fontface = "bold",
                   col = "#66C2A5") +
         
         geom_text(data = tail_prob %>%
                     filter(m == countries_cluster$countries[i],
                            method == "rw"),
                   aes(label = paste0("rw: ", tail_prob_mean_SD_prop, "%"),
                       x = Inf, y = Inf),
                   hjust = 1, vjust = 4,
                   size = 3,
                   fontface = "bold",
                   col = "#8DA0CB") +
         
         geom_text(data = tail_prob %>%
                     filter(m == countries_cluster$countries[i],
                            method == "dlm"),
                   aes(label = paste0("dlm: ", tail_prob_mean_SD_prop, "%"),
                       x = Inf, y = Inf),
                   hjust = 1, vjust = 6,
                   size = 3,
                   fontface = "bold",
                   col =  "#E78AC3") +
         
         facet_wrap(~v,scale="free")+
         scale_fill_manual(name = "", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
         labs(title=paste0(countries_cluster$countries[i], " - Cluster ", countries_cluster$cluster[i]))
  )
}
dev.off()


#### SD UPPER Q - COUNT ####

pdf("output/SD_upper_q/countries_upper_q_SD_count_log-tail.pdf",width=8,height=4.5)
for(i in 1:N_countries){
  plot(ggplot(plot_tail_prob %>%
                filter(m == countries_cluster$countries[i]) %>% 
                mutate(across(total:SD_lower_q_prop, ~if_else(.x == 0, 0.001, .x))),
              aes(x = log(SD_upper_q_count), fill = method))+
         
         geom_density(alpha=0.5)+
         
         geom_vline(data = filter(data_stats, 
                                  m == countries_cluster$countries[i]),
                    aes(xintercept = log(SD_upper_q_count)),
                    col ="#FC8D62", linewidth = 0.8)+
         
         geom_text(data = tail_prob %>%
                     filter(m == countries_cluster$countries[i],
                            method == "hmm"),
                   aes(label = paste0("hmm: ", tail_prob_SD_upper_q_count, "%"),
                       x = 0, y = Inf),
                   hjust = 1, vjust = 2,
                   size = 3,
                   fontface = "bold",
                   col = "#66C2A5") +
         
         geom_text(data = tail_prob %>%
                     filter(m == countries_cluster$countries[i],
                            method == "rw"),
                   aes(label = paste0("rw: ", tail_prob_SD_upper_q_count, "%"),
                       x = 0, y = Inf),
                   hjust = 1, vjust = 4,
                   size = 3,
                   fontface = "bold",
                   col = "#8DA0CB") +
         
         geom_text(data = tail_prob %>%
                     filter(m == countries_cluster$countries[i],
                            method == "dlm"),
                   aes(label = paste0("dlm: ", tail_prob_SD_upper_q_count, "%"),
                       x = 0, y = Inf),
                   hjust = 1, vjust = 6,
                   size = 3,
                   fontface = "bold",
                   col = "#E78AC3") +
         
         facet_wrap(~v,scale="free")+
         scale_fill_manual(name = "", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
         labs(title=paste0(countries_cluster$countries[i], " - Cluster ", countries_cluster$cluster[i]))
  )
}
dev.off()


#### SD UPPER Q - PROPORTION ####

pdf("output/SD_upper_q/countries_upper_q_SD_proportion_log-tail.pdf",width=8,height=4.5)
for(i in 1:N_countries){
  plot(ggplot(plot_tail_prob %>%
                filter(m == countries_cluster$countries[i]) %>% 
                mutate(across(total:SD_lower_q_prop, ~if_else(.x == 0, 0.001, .x))),
              aes(x = log(SD_upper_q_prop), fill = method))+
         
         geom_density(alpha=0.5)+
         
         geom_vline(data = filter(data_stats, 
                                  m == countries_cluster$countries[i]),
                    aes(xintercept = log(SD_upper_q_prop)),
                    col ="#FC8D62", linewidth = 0.8)+
         
         geom_text(data = tail_prob %>%
                     filter(m == countries_cluster$countries[i],
                            method == "hmm"),
                   aes(label = paste0("hmm: ", tail_prob_SD_upper_q_prop, "%"),
                       x = Inf, y = Inf),
                   hjust = 1, vjust = 2,
                   size = 3,
                   fontface = "bold",
                   col = "#66C2A5") +
         
         geom_text(data = tail_prob %>%
                     filter(m == countries_cluster$countries[i],
                            method == "rw"),
                   aes(label = paste0("rw: ", tail_prob_SD_upper_q_prop, "%"),
                       x = Inf, y = Inf),
                   hjust = 1, vjust = 4,
                   size = 3,
                   fontface = "bold",
                   col = "#8DA0CB") +
         
         geom_text(data = tail_prob %>%
                     filter(m == countries_cluster$countries[i],
                            method == "dlm"),
                   aes(label = paste0("dlm: ", tail_prob_SD_upper_q_prop, "%"),
                       x = Inf, y = Inf),
                   hjust = 1, vjust = 6,
                   size = 3,
                   fontface = "bold",
                   col = "#E78AC3") +
         
         facet_wrap(~v,scale="free")+
         scale_fill_manual(name = "", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
         labs(title=paste0(countries_cluster$countries[i], " - Cluster ", countries_cluster$cluster[i]))
  ) 
}
dev.off()


#### SD LOWER Q - COUNT ####

pdf("output/SD_lower_q/countries_lower_q_SD_count_log-tail.pdf",width=8,height=4.5)
for(i in 1:N_countries){
  plot(ggplot(plot_tail_prob %>%
                filter(m == countries_cluster$countries[i]) %>% 
                mutate(across(total:SD_lower_q_prop, ~if_else(.x == 0, 0.001, .x))),
              aes(x = log(SD_lower_q_count), fill = method))+
         
         geom_density(alpha=0.5)+
         
         geom_vline(data = filter(data_stats, 
                                  m == countries_cluster$countries[i]),
                    aes(xintercept = log(SD_lower_q_count)),
                    col ="#FC8D62", linewidth = 0.8)+
         
         geom_text(data = tail_prob %>%
                     filter(m == countries_cluster$countries[i],
                            method == "hmm"),
                   aes(label = paste0("hmm: ", tail_prob_SD_lower_q_count, "%"),
                       x = 0, y = Inf),
                   hjust = 1, vjust = 2,
                   size = 3,
                   fontface = "bold",
                   col = "#66C2A5") +
         
         geom_text(data = tail_prob %>%
                     filter(m == countries_cluster$countries[i],
                            method == "rw"),
                   aes(label = paste0("rw: ", tail_prob_SD_lower_q_count, "%"),
                       x = 0, y = Inf),
                   hjust = 1, vjust = 4,
                   size = 3,
                   fontface = "bold",
                   col = "#8DA0CB") +
         
         geom_text(data = tail_prob %>%
                     filter(m == countries_cluster$countries[i],
                            method == "dlm"),
                   aes(label = paste0("dlm: ", tail_prob_SD_lower_q_count, "%"),
                       x = 0, y = Inf),
                   hjust = 1, vjust = 6,
                   size = 3,
                   fontface = "bold",
                   col = "#E78AC3") +
         
         facet_wrap(~v,scale="free")+
         scale_fill_manual(name = "", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
         labs(title=paste0(countries_cluster$countries[i], " - Cluster ", countries_cluster$cluster[i]))
  )
}
dev.off()


#### SD LOWER Q - PROPORTION ####

pdf("output/SD_lower_q/countries_lower_q_SD_proportion_log-tail.pdf",width=8,height=4.5)
for(i in 1:N_countries){
  plot(ggplot(plot_tail_prob %>%
                filter(m == countries_cluster$countries[i]) %>% 
                mutate(across(total:SD_lower_q_prop, ~if_else(.x == 0, 0.001, .x))),
              aes(x = log(SD_lower_q_prop), fill = method))+
         
         geom_density(alpha=0.5)+
         
         geom_vline(data = filter(data_stats, 
                                  m == countries_cluster$countries[i]),
                    aes(xintercept = log(SD_lower_q_prop)),
                    col ="#FC8D62", linewidth = 0.8)+
         
         geom_text(data = tail_prob %>%
                     filter(m == countries_cluster$countries[i],
                            method == "hmm"),
                   aes(label = paste0("hmm: ", tail_prob_SD_lower_q_prop, "%"),
                       x = Inf, y = Inf),
                   hjust = 1, vjust = 2,
                   size = 3,
                   fontface = "bold",
                   col = "#66C2A5") +
         
         geom_text(data = tail_prob %>%
                     filter(m == countries_cluster$countries[i],
                            method == "rw"),
                   aes(label = paste0("rw: ", tail_prob_SD_lower_q_prop, "%"),
                       x = Inf, y = Inf),
                   hjust = 1, vjust = 4,
                   size = 3,
                   fontface = "bold",
                   col = "#8DA0CB") +
         
         geom_text(data = tail_prob %>%
                     filter(m == countries_cluster$countries[i],
                            method == "dlm"),
                   aes(label = paste0("dlm: ", tail_prob_SD_lower_q_prop, "%"),
                       x = Inf, y = Inf),
                   hjust = 1, vjust = 6,
                   size = 3,
                   fontface = "bold",
                   col = "#E78AC3") +
         
         facet_wrap(~v,scale="free")+
         scale_fill_manual(name = "", values = c("hmm" = "#66C2A5",
                                                 "rw" = "#8DA0CB",
                                                 "dlm" = "#E78AC3")) +
         labs(title=paste0(countries_cluster$countries[i], " - Cluster ", countries_cluster$cluster[i]))
  )
}
dev.off()


#######################.

## WINDOW LENGTHS ####

### DATA ####

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


### HMM ####

hmm_replicates_stats_window_5 <- hmm_replicates %>%
  
  group_by(v, m, sample)%>%
  partition(cluster)%>%
  
  summarise(mean_SD_count    = mean(window_stats(x,L[1],N,sd)),
            SD_upper_q_count = quantile(window_stats(x,L[1],N,sd), 0.75),
            SD_lower_q_count = quantile(window_stats(x,L[1],N,sd), 0.25),
            
            mean_SD_prop    = mean(window_stats(proportion,L[1],N,sd)),
            SD_upper_q_prop = quantile(window_stats(proportion,L[1],N,sd), 0.75),
            SD_lower_q_prop = quantile(window_stats(proportion,L[1],N,sd), 0.25)) %>%
  collect()

hmm_replicates_stats_window_10 <- hmm_replicates %>%
  
  group_by(v, m, sample)%>%
  partition(cluster)%>%
  
  summarise(mean_SD_count    = mean(window_stats(x,L[2],N,sd)),
            SD_upper_q_count = quantile(window_stats(x,L[2],N,sd), 0.75),
            SD_lower_q_count = quantile(window_stats(x,L[2],N,sd), 0.25),
            
            mean_SD_prop    = mean(window_stats(proportion,L[2],N,sd)),
            SD_upper_q_prop = quantile(window_stats(proportion,L[2],N,sd), 0.75),
            SD_lower_q_prop = quantile(window_stats(proportion,L[2],N,sd), 0.25)) %>%
  collect()

hmm_replicates_stats_window_15 <- hmm_replicates %>%
  
  group_by(v, m, sample)%>%
  partition(cluster)%>%
  
  summarise(mean_SD_count    = mean(window_stats(x,L[3],N,sd)),
            SD_upper_q_count = quantile(window_stats(x,L[3],N,sd), 0.75),
            SD_lower_q_count = quantile(window_stats(x,L[3],N,sd), 0.25),
            
            mean_SD_prop    = mean(window_stats(proportion,L[3],N,sd)),
            SD_upper_q_prop = quantile(window_stats(proportion,L[3],N,sd), 0.75),
            SD_lower_q_prop = quantile(window_stats(proportion,L[3],N,sd), 0.25)) %>%
  collect()

hmm_replicates_stats_window_20 <- hmm_replicates %>%
  
  group_by(v, m, sample)%>%
  partition(cluster)%>%
  
  summarise(mean_SD_count    = mean(window_stats(x,L[4],N,sd)),
            SD_upper_q_count = quantile(window_stats(x,L[4],N,sd), 0.75),
            SD_lower_q_count = quantile(window_stats(x,L[4],N,sd), 0.25),
            
            mean_SD_prop    = mean(window_stats(proportion,L[4],N,sd)),
            SD_upper_q_prop = quantile(window_stats(proportion,L[4],N,sd), 0.75),
            SD_lower_q_prop = quantile(window_stats(proportion,L[4],N,sd), 0.25)) %>%
  collect()


#### COMBINE WINDOW STATS #####

hmm_replicates_stats_window <- hmm_replicates_stats_window_5 %>%
  pivot_longer(., cols = mean_SD_count_5:SD_lower_q_prop_5, names_to = "window_length", values_to = "value") %>%
  mutate(measure = window_length, .before = "value") %>%
  mutate(window_length = gsub(".*?(\\d+)", "\\1", window_length),
         measure = gsub("_\\d+", "", measure)) %>%
  rbind(hmm_replicates_stats_window_10 %>%
          pivot_longer(., cols = mean_SD_count_10:SD_lower_q_prop_10, names_to = "window_length", values_to = "value") %>%
          mutate(measure = window_length, .before = "value") %>%
          mutate(window_length = gsub(".*?(\\d+)", "\\1", window_length),
                 measure = gsub("_\\d+", "", measure))) %>%
  rbind(hmm_replicates_stats_window_15 %>%
          pivot_longer(., cols = mean_SD_count_15:SD_lower_q_prop_15, names_to = "window_length", values_to = "value") %>%
          mutate(measure = window_length, .before = "value") %>%
          mutate(window_length = gsub(".*?(\\d+)", "\\1", window_length),
                 measure = gsub("_\\d+", "", measure))) %>%
  rbind(hmm_replicates_stats_window_20 %>%
          pivot_longer(., cols = mean_SD_count_20:SD_lower_q_prop_20, names_to = "window_length", values_to = "value") %>%
          mutate(measure = window_length, .before = "value") %>%
          mutate(window_length = gsub(".*?(\\d+)", "\\1", window_length),
                 measure = gsub("_\\d+", "", measure)))


### RW ####

rw_replicates_stats_window_5 <- rw_replicates %>%
  
  group_by(v, m, sample)%>%
  partition(cluster)%>%
  
  summarise(mean_SD_count    = mean(window_stats(x,L[1],N,sd)),
            SD_upper_q_count = quantile(window_stats(x,L[1],N,sd), 0.75),
            SD_lower_q_count = quantile(window_stats(x,L[1],N,sd), 0.25),
            
            mean_SD_prop    = mean(window_stats(proportion,L[1],N,sd)),
            SD_upper_q_prop = quantile(window_stats(proportion,L[1],N,sd), 0.75),
            SD_lower_q_prop = quantile(window_stats(proportion,L[1],N,sd), 0.25)) %>%
  collect()

rw_replicates_stats_window_10 <- rw_replicates %>%
  
  group_by(v, m, sample)%>%
  partition(cluster)%>%
  
  summarise(mean_SD_count    = mean(window_stats(x,L[2],N,sd)),
            SD_upper_q_count = quantile(window_stats(x,L[2],N,sd), 0.75),
            SD_lower_q_count = quantile(window_stats(x,L[2],N,sd), 0.25),
            
            mean_SD_prop    = mean(window_stats(proportion,L[2],N,sd)),
            SD_upper_q_prop = quantile(window_stats(proportion,L[2],N,sd), 0.75),
            SD_lower_q_prop = quantile(window_stats(proportion,L[2],N,sd), 0.25)) %>%
  collect()

rw_replicates_stats_window_15 <- rw_replicates %>%
  
  group_by(v, m, sample)%>%
  partition(cluster)%>%
  
  summarise(mean_SD_count    = mean(window_stats(x,L[3],N,sd)),
            SD_upper_q_count = quantile(window_stats(x,L[3],N,sd), 0.75),
            SD_lower_q_count = quantile(window_stats(x,L[3],N,sd), 0.25),
            
            mean_SD_prop    = mean(window_stats(proportion,L[3],N,sd)),
            SD_upper_q_prop = quantile(window_stats(proportion,L[3],N,sd), 0.75),
            SD_lower_q_prop = quantile(window_stats(proportion,L[3],N,sd), 0.25)) %>%
  collect()

rw_replicates_stats_window_20 <- rw_replicates %>%
  
  group_by(v, m, sample)%>%
  partition(cluster)%>%
  
  summarise(mean_SD_count    = mean(window_stats(x,L[4],N,sd)),
            SD_upper_q_count = quantile(window_stats(x,L[4],N,sd), 0.75),
            SD_lower_q_count = quantile(window_stats(x,L[4],N,sd), 0.25),
            
            mean_SD_prop    = mean(window_stats(proportion,L[4],N,sd)),
            SD_upper_q_prop = quantile(window_stats(proportion,L[4],N,sd), 0.75),
            SD_lower_q_prop = quantile(window_stats(proportion,L[4],N,sd), 0.25)) %>%
  collect()


#### COMBINE WINDOW STATS #####

rw_replicates_stats_window <- rw_replicates_stats_window_5 %>%
  pivot_longer(., cols = mean_SD_count_5:SD_lower_q_prop_5, names_to = "window_length", values_to = "value") %>%
  mutate(measure = window_length, .before = "value") %>%
  mutate(window_length = gsub(".*?(\\d+)", "\\1", window_length),
         measure = gsub("_\\d+", "", measure)) %>%
  rbind(rw_replicates_stats_window_10 %>%
          pivot_longer(., cols = mean_SD_count_10:SD_lower_q_prop_10, names_to = "window_length", values_to = "value") %>%
          mutate(measure = window_length, .before = "value") %>%
          mutate(window_length = gsub(".*?(\\d+)", "\\1", window_length),
                 measure = gsub("_\\d+", "", measure))) %>%
  rbind(rw_replicates_stats_window_15 %>%
          pivot_longer(., cols = mean_SD_count_15:SD_lower_q_prop_15, names_to = "window_length", values_to = "value") %>%
          mutate(measure = window_length, .before = "value") %>%
          mutate(window_length = gsub(".*?(\\d+)", "\\1", window_length),
                 measure = gsub("_\\d+", "", measure))) %>%
  rbind(rw_replicates_stats_window_20 %>%
          pivot_longer(., cols = mean_SD_count_20:SD_lower_q_prop_20, names_to = "window_length", values_to = "value") %>%
          mutate(measure = window_length, .before = "value") %>%
          mutate(window_length = gsub(".*?(\\d+)", "\\1", window_length),
                 measure = gsub("_\\d+", "", measure)))


### DLM ####

dlm_replicates_stats_window_5 <- dlm_replicates %>%
  
  group_by(v, m, sample)%>%
  partition(cluster)%>%
  
  summarise(mean_SD_count    = mean(window_stats(x,L[1],N,sd)),
            SD_upper_q_count = quantile(window_stats(x,L[1],N,sd), 0.75),
            SD_lower_q_count = quantile(window_stats(x,L[1],N,sd), 0.25),
            
            mean_SD_prop    = mean(window_stats(proportion,L[1],N,sd)),
            SD_upper_q_prop = quantile(window_stats(proportion,L[1],N,sd), 0.75),
            SD_lower_q_prop = quantile(window_stats(proportion,L[1],N,sd), 0.25)) %>%
  collect()

dlm_replicates_stats_window_10 <- dlm_replicates %>%
  
  group_by(v, m, sample)%>%
  partition(cluster)%>%
  
  summarise(mean_SD_count    = mean(window_stats(x,L[2],N,sd)),
            SD_upper_q_count = quantile(window_stats(x,L[2],N,sd), 0.75),
            SD_lower_q_count = quantile(window_stats(x,L[2],N,sd), 0.25),
            
            mean_SD_prop    = mean(window_stats(proportion,L[2],N,sd)),
            SD_upper_q_prop = quantile(window_stats(proportion,L[2],N,sd), 0.75),
            SD_lower_q_prop = quantile(window_stats(proportion,L[2],N,sd), 0.25)) %>%
  collect()

dlm_replicates_stats_window_15 <- dlm_replicates %>%
  
  group_by(v, m, sample)%>%
  partition(cluster)%>%
  
  summarise(mean_SD_count    = mean(window_stats(x,L[3],N,sd)),
            SD_upper_q_count = quantile(window_stats(x,L[3],N,sd), 0.75),
            SD_lower_q_count = quantile(window_stats(x,L[3],N,sd), 0.25),
            
            mean_SD_prop    = mean(window_stats(proportion,L[3],N,sd)),
            SD_upper_q_prop = quantile(window_stats(proportion,L[3],N,sd), 0.75),
            SD_lower_q_prop = quantile(window_stats(proportion,L[3],N,sd), 0.25)) %>%
  collect()

dlm_replicates_stats_window_20 <- dlm_replicates %>%
  
  group_by(v, m, sample)%>%
  partition(cluster)%>%
  
  summarise(mean_SD_count    = mean(window_stats(x,L[4],N,sd)),
            SD_upper_q_count = quantile(window_stats(x,L[4],N,sd), 0.75),
            SD_lower_q_count = quantile(window_stats(x,L[4],N,sd), 0.25),
            
            mean_SD_prop    = mean(window_stats(proportion,L[4],N,sd)),
            SD_upper_q_prop = quantile(window_stats(proportion,L[4],N,sd), 0.75),
            SD_lower_q_prop = quantile(window_stats(proportion,L[4],N,sd), 0.25)) %>%
  collect()


#### COMBINE WINDOW STATS #####

dlm_replicates_stats_window <- dlm_replicates_stats_window_5 %>%
  pivot_longer(., cols = mean_SD_count_5:SD_lower_q_prop_5, names_to = "window_length", values_to = "value") %>%
  mutate(measure = window_length, .before = "value") %>%
  mutate(window_length = gsub(".*?(\\d+)", "\\1", window_length),
         measure = gsub("_\\d+", "", measure)) %>%
  rbind(dlm_replicates_stats_window_10 %>%
          pivot_longer(., cols = mean_SD_count_10:SD_lower_q_prop_10, names_to = "window_length", values_to = "value") %>%
          mutate(measure = window_length, .before = "value") %>%
          mutate(window_length = gsub(".*?(\\d+)", "\\1", window_length),
                 measure = gsub("_\\d+", "", measure))) %>%
  rbind(dlm_replicates_stats_window_15 %>%
          pivot_longer(., cols = mean_SD_count_15:SD_lower_q_prop_15, names_to = "window_length", values_to = "value") %>%
          mutate(measure = window_length, .before = "value") %>%
          mutate(window_length = gsub(".*?(\\d+)", "\\1", window_length),
                 measure = gsub("_\\d+", "", measure))) %>%
  rbind(dlm_replicates_stats_window_20 %>%
          pivot_longer(., cols = mean_SD_count_20:SD_lower_q_prop_20, names_to = "window_length", values_to = "value") %>%
          mutate(measure = window_length, .before = "value") %>%
          mutate(window_length = gsub(".*?(\\d+)", "\\1", window_length),
                 measure = gsub("_\\d+", "", measure)))


### COMBINE METHODS ####

window_summary <- hmm_replicates_stats_window %>%
  mutate(method = "hmm") %>%
  rbind(rw_replicates_stats_window %>% 
          mutate(method = "rw")) %>%
  rbind(dlm_replicates_stats_window %>%
          mutate(method = "dlm")) %>%
  left_join(countries_cluster %>%
              rename(m = "countries") %>%
              mutate(m_numeric = as.numeric(factor(m, levels = unique(m)))) %>%
              select(-m) %>%
              rename(m = "m_numeric"), by = "m") %>%
  mutate(window_length = as.numeric(window_length)) %>%
  mutate(v=variants[v], m=countries_cluster$countries[m])
  


saveRDS(window_summary, "output/window_summary.rds")
#window_summary <- readRDS("output/window_summary.rds")


#######################.

### WINDOW LENGTH MAE ####

window_summary_mae <- window_summary %>%
  left_join(data_stats_window %>%
              pivot_longer(., cols = mean_SD_count_5:SD_lower_q_prop_20, names_to = "window_length", values_to = "original_data_value") %>%
              mutate(measure = window_length, .before = "original_data_value") %>%
              mutate(window_length = gsub(".*?(\\d+)", "\\1", window_length),
                     window_length = as.numeric(window_length),
                     measure = gsub("_\\d+", "", measure)), 
            by = c("v", "m", "window_length", "measure")) %>%
  
  group_by(v, m, method, measure, window_length, cluster) %>%
  summarise(mae = mean(abs(value - original_data_value))) %>%
  ungroup() 


#### PLOT ####

##### MEAN SD - COUNT ####

pdf("output/mae_mean_SD_count-line.pdf", width = 12, height = 6)
for (i in 1:N_clusters) {
  
  plot(ggplot(data = window_summary_mae %>%
                filter(measure == "mean_SD_count",
                       cluster == i) %>%
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
         
         labs(color  = "Method", linetype = "Method") +
         ylab("MAE Mean SD Count") +
         xlab("Window Length") +
         ggtitle(paste0("MAE - Mean SD Count - Cluster ", i, " by Window Length"))
  )
}
dev.off()


pdf("output/mae_mean_SD_count.pdf",width=12,height=6)
for (i in 1:N_clusters) {
  for (j in 1:length(unique(window_summary_mae$window_length))) {
    
    plot(ggplot(data = window_summary_mae %>%
                  filter(cluster == i,
                         window_length == unique(window_summary_mae$window_length)[j],
                         measure == "mean_SD_count"),
                
                aes(x = m,
                    y = mae,
                    fill = method)) +
           
           geom_col(position=position_dodge()) +
           
           geom_hline(data = window_summary_mae %>% 
                        filter(cluster == i,
                               window_length == unique(window_summary_mae$window_length)[j],
                               measure == "mean_SD_count") %>%
                        group_by(v, method, window_length) %>% 
                        summarise(mae = mean(mae)),
                      
                      aes(yintercept = mae, 
                          col = method,
                          linetype = method),
                      size = 1) +
           
           facet_wrap(~v, scale = "free") +
           scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                         "rw" = "#8DA0CB",
                                                         "dlm" = "#E78AC3")) +
           scale_color_manual(name = "Method", values = c("hmm" = "#FC8D62",
                                                          "rw" = "#FFD92F",
                                                          "dlm" = "#E5C494")) +
           
           theme(axis.text.x=element_text(angle=60, hjust=1)) +
           labs(color  = "Method", fill = "Method", linetype = "Method") +
           ylab("MAE") +
           xlab("Country") +
           ggtitle(paste0("MAE - Mean SD Count - Cluster ", i, " - Window Length ", unique(window_summary_mae$window_length)[j]))
    )
  }
}
dev.off()



##### MEAN SD - PROPORTION ####

pdf("output/mae_mean_SD_proportion-line.pdf", width = 12, height = 6)
for (i in 1:N_clusters) {
  
  plot(ggplot(data = window_summary_mae %>%
                filter(measure == "mean_SD_prop",
                       cluster == i) %>%
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
         
         labs(color  = "Method", linetype = "Method") +
         ylab("MAE Mean SD Proportion") +
         xlab("Country") +
         ggtitle(paste0("MAE - Mean SD Proportion - Cluster ", i, " by Window Length"))
  )
}
dev.off()


pdf("output/mae_mean_SD_proportion.pdf",width=12,height=6)
for (i in 1:N_clusters) {
  for (j in 1:length(unique(window_summary_mae$window_length))) {
    
    plot(ggplot(data = window_summary_mae %>%
                  filter(cluster == i,
                         window_length == unique(window_summary_mae$window_length)[j],
                         measure == "mean_SD_prop"),
                
                aes(x = m,
                    y = mae,
                    fill = method)) +
           
           geom_col(position=position_dodge()) +
           
           geom_hline(data = window_summary_mae %>% 
                        filter(cluster == i,
                               window_length == unique(window_summary_mae$window_length)[j],
                               measure == "mean_SD_prop") %>%
                        group_by(v, method, window_length) %>% 
                        summarise(mae = mean(mae)),
                      
                      aes(yintercept = mae, 
                          col = method,
                          linetype = method),
                      size = 1) +
           
           facet_wrap(~v, scale = "free") +
           scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                         "rw" = "#8DA0CB",
                                                         "dlm" = "#E78AC3")) +
           scale_color_manual(name = "Method", values = c("hmm" = "#FC8D62",
                                                          "rw" = "#FFD92F",
                                                          "dlm" = "#E5C494")) +
           
           theme(axis.text.x=element_text(angle=60, hjust=1)) +
           labs(color  = "Method", fill = "Method", linetype = "Method") +
           ylab("MAE") +
           xlab("Country") +
           ggtitle(paste0("MAE - Mean SD Proportion - Cluster ", i, " - Window Length ", unique(window_summary_mae$window_length)[j]))
    )
  }
}
dev.off()


##### SD UPPER Q - COUNT ####

pdf("output/mae_upper_q_SD_count-line.pdf", width = 12, height = 6)
for (i in 1:N_clusters) {
  
  plot(ggplot(data = window_summary_mae %>%
                filter(measure == "SD_upper_q_count",
                       cluster == i) %>%
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
         
         labs(color  = "Method", linetype = "Method") +
         ylab("MAE Mean SD Count") +
         xlab("Country") +
         ggtitle(paste0("MAE - Upper Q SD Count - Cluster ", i, " by Window Length"))
  )
}
dev.off()


pdf("output/mae_upper_q_SD_count.pdf",width=12,height=6)
for (i in 1:N_clusters) {
  for (j in 1:length(unique(window_summary_mae$window_length))) {
    
    plot(ggplot(data = window_summary_mae %>%
                  filter(cluster == i,
                         window_length == unique(window_summary_mae$window_length)[j],
                         measure == "SD_upper_q_count"),
                
                aes(x = m,
                    y = mae,
                    fill = method)) +
           
           geom_col(position=position_dodge()) +
           
           geom_hline(data = window_summary_mae %>% 
                        filter(cluster == i,
                               window_length == unique(window_summary_mae$window_length)[j],
                               measure == "SD_upper_q_count") %>%
                        group_by(v, method, window_length) %>% 
                        summarise(mae = mean(mae)),
                      
                      aes(yintercept = mae, 
                          col = method,
                          linetype = method),
                      size = 1) +
           
           facet_wrap(~v, scale = "free") +
           scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                         "rw" = "#8DA0CB",
                                                         "dlm" = "#E78AC3")) +
           scale_color_manual(name = "Method", values = c("hmm" = "#FC8D62",
                                                          "rw" = "#FFD92F",
                                                          "dlm" = "#E5C494")) +
           
           theme(axis.text.x=element_text(angle=60, hjust=1)) +
           labs(color  = "Method", fill = "Method", linetype = "Method") +
           ylab("MAE") +
           xlab("Country") +
           ggtitle(paste0("MAE - Upper Q SD Count - Cluster ", i, " - Window Length ", unique(window_summary_mae$window_length)[j]))
    )
  }
}
dev.off()



##### SD UPPER Q - PROPORTION ####

pdf("output/mae_upper_q_SD_proprotion-line.pdf", width = 12, height = 6)
for (i in 1:N_clusters) {
  
  plot(ggplot(data = window_summary_mae %>%
                filter(measure == "SD_upper_q_prop",
                       cluster == i) %>%
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
         
         labs(color  = "Method", linetype = "Method") +
         ylab("MAE Mean SD Proportion") +
         xlab("Country") +
         ggtitle(paste0("MAE - Upper Q SD Proportion - Cluster ", i, " by Window Length"))
  )
}
dev.off()


pdf("output/mae_upper_q_SD_proportion.pdf",width=12,height=6)
for (i in 1:N_clusters) {
  for (j in 1:length(unique(window_summary_mae$window_length))) {
    
    plot(ggplot(data = window_summary_mae %>%
                  filter(cluster == i,
                         window_length == unique(window_summary_mae$window_length)[j],
                         measure == "SD_upper_q_prop"),
                
                aes(x = m,
                    y = mae,
                    fill = method)) +
           
           geom_col(position=position_dodge()) +
           
           geom_hline(data = window_summary_mae %>% 
                        filter(cluster == i,
                               window_length == unique(window_summary_mae$window_length)[j],
                               measure == "SD_upper_q_prop") %>%
                        group_by(v, method, window_length) %>% 
                        summarise(mae = mean(mae)),
                      
                      aes(yintercept = mae, 
                          col = method,
                          linetype = method),
                      size = 1) +
           
           facet_wrap(~v, scale = "free") +
           scale_fill_manual(name = "Method", values = c("hmm" = "#66C2A5",
                                                         "rw" = "#8DA0CB",
                                                         "dlm" = "#E78AC3")) +
           scale_color_manual(name = "Method", values = c("hmm" = "#FC8D62",
                                                          "rw" = "#FFD92F",
                                                          "dlm" = "#E5C494")) +
           
           theme(axis.text.x=element_text(angle=60, hjust=1)) +
           labs(color  = "Method", fill = "Method", linetype = "Method") +
           ylab("MAE") +
           xlab("Country") +
           ggtitle(paste0("MAE - Upper Q SD Proportion - Cluster ", i, " - Window Length ", unique(window_summary_mae$window_length)[j]))
    )
  }
}
dev.off()


#######################.
