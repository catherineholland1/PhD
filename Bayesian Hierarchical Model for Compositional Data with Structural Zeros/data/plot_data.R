#####################.
## GLASS DATA PLOT ##
#####################.

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)

##################################.

setwd("data/")

##################################.

## DATA ####

glass <- read.table("", header = T)

glass$Type <- as.factor(sapply(glass$Name, function(x) substr(x, 1, 1)))


ratio_data <- glass %>%
  mutate(across(Na:Fe, ~ (.x/O))) %>% 
  select(Na:Fe, Item, Type, piece) %>%
  rename(Piece = piece)


sqrt_data <- glass %>%
  mutate(across(Na:Fe, ~ sqrt(.x/O))) %>% 
  select(Na:Fe, Item, Type, piece) %>%
  rename(Piece = piece)


##################################.

## PLOT ####

### PAIRS ####

ggpairs(data = sqrt_data,
        columns = 1:7,
        legend = 1, 
        aes(color = factor(Type), 
            labels = c("Glass Type")), 
        upper = NULL,
        lower = list(continuous = wrap("points", alpha = 0.5))) + 
  theme(legend.position = "bottom") +  
  scale_fill_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854"), 
                    name = 'Glass Type', 
                    labels = c('bulb', 'car window', 'headlamp', 'container', 'building window')) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854")) +
  theme(axis.text.x = element_text(angle = 45)) 


### BOXPLOT ####

pdf('output/full_data_boxplots.pdf', height = 5)
ratio_data %>%
  pivot_longer(cols = c(Na:Fe), names_to = 'element', values_to = 'value') %>%
  mutate(element = factor(element, levels = c('Na', 'Mg', 'Alu', 'Si', 'K', 'Ca', 'Fe')),
         Type = case_when(Type == 'b' ~ 'bulb',
                          Type == 'c' ~ 'car window',
                          Type == 'h' ~ 'headlamp',
                          Type == 'p' ~ 'container',
                          Type == 'w' ~ 'building window'),
         Type = factor(Type, levels = c('bulb', 'car window', 'headlamp', 'container', 'building window'))) %>%
  ggplot(., aes(x = factor(Type), y = value, fill = factor(Type))) +
  geom_boxplot() +
  facet_grid(~element) +
  labs(x = '', y = '') +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854"), 
                    name = 'Glass Type', 
                    labels = c('bulb', 'car window', 'headlamp', 'container', 'building window')) +
  theme(axis.text.x=element_text(angle=60, hjust=1),
        legend.position = 'bottom')
dev.off()

pdf('output/sqrt_full_data_boxplots.pdf', height = 5)
sqrt_data %>%
  pivot_longer(cols = c(Na:Fe), names_to = 'element', values_to = 'value') %>%
  mutate(element = factor(element, levels = c('Na', 'Mg', 'Alu', 'Si', 'K', 'Ca', 'Fe')),
         Type = case_when(Type == 'b' ~ 'bulb',
                          Type == 'c' ~ 'car window',
                          Type == 'h' ~ 'headlamp',
                          Type == 'p' ~ 'container',
                          Type == 'w' ~ 'building window'),
         Type = factor(Type, levels = c('bulb', 'car window', 'headlamp', 'container', 'building window'))) %>%
  ggplot(., aes(x = factor(Type), y = value, fill = factor(Type))) +
  geom_boxplot() +
  facet_grid(~element) +
  labs(x = '', y = '') +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854"), 
                    name = 'Glass Type', 
                    labels = c('bulb', 'car window', 'headlamp', 'container', 'building window')) +
  theme(axis.text.x=element_text(angle=60, hjust=1),
        legend.position = 'bottom')
dev.off()


pdf('output/grid_sqrt_full_data_boxplots.pdf', height = 4)
ratio_data %>%
  mutate(transform = 'untransformed') %>%
  rbind(sqrt_data %>% mutate(transform = 'square root')) %>%
  pivot_longer(cols = c(Na:Fe), names_to = 'element', values_to = 'value') %>%
  mutate(element = factor(element, levels = c('Na', 'Mg', 'Alu', 'Si', 'K', 'Ca', 'Fe')),
         Type = case_when(Type == 'b' ~ 'bulb',
                          Type == 'c' ~ 'car window',
                          Type == 'h' ~ 'headlamp',
                          Type == 'p' ~ 'container',
                          Type == 'w' ~ 'building window'),
         Type = factor(Type, levels = c('bulb', 'car window', 'headlamp', 'container', 'building window')),
         transform = factor(transform, levels = c('untransformed', 'square root'))) %>%
  ggplot(., aes(x = factor(Type), y = value, fill = factor(Type))) +
  geom_boxplot(outliers = F) +
  facet_grid2(transform~element, scales='free', independent='y') +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  labs(x = '', y = '') +
  scale_x_discrete(labels = c('','','','','')) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854"), 
                    name = 'Glass Type', 
                    labels = c('bulb', 'car window', 'headlamp', 'container', 'building window')) +
  theme(#axis.text.x=element_text(angle=60, hjust=1),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 8),
        legend.position = 'bottom',
        legend.margin = margin(t = -15))
dev.off()

pdf('output/scatterplot_si_ca.pdf', height = 4, width = 6)
sqrt_data %>%
  mutate(Type = case_when(Type == 'b' ~ 'bulb',
                          Type == 'c' ~ 'car window',
                          Type == 'h' ~ 'headlamp',
                          Type == 'p' ~ 'container',
                          Type == 'w' ~ 'building window'),
         Type = factor(Type, levels = c('bulb', 'car window', 'headlamp', 'container', 'building window'))) %>%
  ggplot(., aes(x = Si, y = Ca, col = Type)) +
  geom_point() +
  scale_color_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854"), 
                    name = 'Glass Type', 
                    labels = c('bulb', 'car window', 'headlamp', 'container', 'building window')) +
  theme(axis.text.x=element_text(angle=60, hjust=1),
        legend.position = 'bottom')
dev.off()

##################################.
