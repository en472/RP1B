
# Analysing Plasmid Diversity in Clinical Samples

## Libraries ===================================================

library(tidyverse)
library(ggplot2)
library(patchwork)

## Read in data ================================================

meta <- read.csv('meta_data.csv')

distances <- read.csv('distances.csv')

# drop first columns for formatting

meta <- meta[, 2:ncol(meta)]

distances <- distances[, 2:ncol(distances)]

## Clinical vs. Non-Clinical origin plasmid diversity by country ===

# Data formatting ---------------------------------------------------

# Extract italy, thai, and argentina data separately

argentina_only <- distances %>%
  filter(grepl('KLAR', id1) == TRUE & grepl('KLAR', id2) == TRUE)

thailand_only <- distances %>%
  filter(grepl('OH-DART', id1) == TRUE & grepl('OH-DART', id2) == TRUE)

italy_only <- distances %>%
  filter(grepl('SPARK', id1) == TRUE & grepl('SPARK', id2) == TRUE)

# join clinical/non-clinical labels to thai and italian data 
# (all argentinian data is clinical)

# select out relevant meta data columns
meta_sample_origin <- meta %>%
  select(id, Sample.type.1)

thailand_only <- thailand_only %>%
  # join labels for plasmid 1
  left_join(meta_sample_origin, by = c('id1' = 'id')) %>%
  # join labels for plasmid 2
  left_join(meta_sample_origin, by = c('id2' = 'id')) %>%
  # filter both only for human samples
  filter(grepl('human', Sample.type.1.x) == TRUE) %>%
  filter(grepl('human', Sample.type.1.y) == TRUE) %>%
  # now group by clinical or non-clinical origin
  # group hospital disease and hospital into the same category in this instance
  mutate(location_s1 = case_when(
    grepl('non-hospital', Sample.type.1.x) == TRUE ~ 'nh',
    .default = 'h'
  )) %>%
  mutate(location_s2 = case_when(
    grepl('non-hospital', Sample.type.1.y) == TRUE ~ 'nh',
    .default = 'h'
  )) %>%
  # now add final label
  mutate(Origin = case_when(
    location_s1 == 'h' & location_s2 == 'h' ~ 'Clinical',
    location_s1 == 'nh' & location_s2 == 'nh' ~ 'Non-Clinical',
    .default = 'Mix'
  )) %>%
  # select out relevant cols
  select(id1, id2, mash_distance, pling_distance, Origin)

# repeat for italian data
italy_only <- italy_only %>%
  # join labels for plasmid 1
  left_join(meta_sample_origin, by = c('id1' = 'id')) %>%
  # join labels for plasmid 2
  left_join(meta_sample_origin, by = c('id2' = 'id')) %>%
  # filter both only for human samples
  filter(grepl('human', Sample.type.1.x) == TRUE) %>%
  filter(grepl('human', Sample.type.1.y) == TRUE) %>%
  # now group by clinical or non-clinical origin
  # group hospital disease and hospital into the same category in this instance
  mutate(location_s1 = case_when(
    grepl('non-hospital', Sample.type.1.x) == TRUE ~ 'nh',
    .default = 'h'
  )) %>%
  mutate(location_s2 = case_when(
    grepl('non-hospital', Sample.type.1.y) == TRUE ~ 'nh',
    .default = 'h'
  )) %>%
  # now add final label
  mutate(Origin = case_when(
    location_s1 == 'h' & location_s2 == 'h' ~ 'Clinical',
    location_s1 == 'nh' & location_s2 == 'nh' ~ 'Non-Clinical',
    .default = 'Mix'
  )) %>%
  # select out relevant cols
  select(id1, id2, mash_distance, pling_distance, Origin)

argentina_only <- argentina_only %>%
  mutate(Origin = 'Clinical')


# plotting ----------------------------------------------------

thailand_only %>%
  ggplot(aes(y = mash_distance, x = Origin, fill = Origin)) +
  geom_violin() +
  labs(
    title = 'MASH Distance of Klebsiella Plamids in Humans - Thailand',
    y = 'Distance',
    x = ''
  ) +
  scale_fill_manual(values = c('#5AB1BB', '#A5C882', '#F7DD72'))

italy_only %>%
  ggplot(aes(y = mash_distance, x = Origin, fill = Origin)) +
  geom_violin() +
  labs(
    title = 'MASH Distance of Klebsiella Plamids in Humans - Italy',
    y = 'Distance',
    x = ''
  ) +
  scale_fill_manual(values = c('#5AB1BB', '#A5C882', '#F7DD72'))

# exclude, nothing to compare between - or do histogram
argentina_only %>%
  ggplot(aes(y = mash_distance, x = Origin, fill = Origin)) +
  geom_violin() +
  labs(
    title = 'MASH Distance of Klebsiella Plamids in Humans - Argentina',
    y = 'Distance',
    x = ''
  ) +
  scale_fill_manual(values = c('#5AB1BB'))


# percentages for write up ---------------
argentina_only %>%
  filter(Origin == 'Clinical') %>%
  filter(mash_distance == 0) %>%
  summarise(
    no_close = n(),
    prop = no_close / nrow(argentina_only %>% filter(Origin == 'Clinical')) * 100
  ) 



## Now Clinical Samples between countries =========================

# remove comparisons of plasmids from the same countries
country_comparisons <- distances %>%
  filter(!(grepl('SPARK', id1) == TRUE & grepl('SPARK', id2) == TRUE)) %>%
  filter(!(grepl('KLAR', id1) == TRUE & grepl('KLAR', id2) == TRUE)) %>%
  filter(!(grepl('OH-DART', id1) == TRUE & grepl('OH-DART', id2) == TRUE)) %>%
  # now label by origin
  # join labels for plasmid 1
  left_join(meta_sample_origin, by = c('id1' = 'id')) %>%
  # join labels for plasmid 2
  left_join(meta_sample_origin, by = c('id2' = 'id')) %>%
  # filter both only for human samples
  filter(grepl('human', Sample.type.1.x) == TRUE) %>%
  filter(grepl('human', Sample.type.1.y) == TRUE) %>%
  # now group by clinical or non-clinical origin
  # group hospital disease and hospital into the same category in this instance
  mutate(location_s1 = case_when(
    grepl('non-hospital', Sample.type.1.x) == TRUE ~ 'nh',
    .default = 'h'
  )) %>%
  mutate(location_s2 = case_when(
    grepl('non-hospital', Sample.type.1.y) == TRUE ~ 'nh',
    .default = 'h'
  )) %>%
  # now add final label
  mutate(Origin = case_when(
    location_s1 == 'h' & location_s2 == 'h' ~ 'Clinical',
    location_s1 == 'nh' & location_s2 == 'nh' ~ 'Non-Clinical',
    .default = 'Mix'
  )) %>%
  # now add country labels
  mutate(id1_country = case_when(
    grepl('SPARK', id1) == TRUE ~ 'Italy',
    grepl('KLAR', id1) == TRUE ~ 'Argentina',
    grepl('OH-DART', id1) == TRUE ~ 'Thailand'
  )) %>%
  mutate(id2_country = case_when(
    grepl('SPARK', id2) == TRUE ~ 'Italy',
    grepl('KLAR', id2) == TRUE ~ 'Argentina',
    grepl('OH-DART', id2) == TRUE ~ 'Thailand'
  )) %>%
  # label according to country
  mutate(Country = case_when(
    (id1_country == 'Argentina' & id2_country == 'Italy') | (id2_country == 'Argentina' & id1_country == 'Italy') ~ 'Argentina/Italy',
    (id1_country == 'Argentina' & id2_country == 'Thailand') | (id2_country == 'Argentina' & id1_country == 'Thailand') ~ 'Argentina/Thailand',
    (id1_country == 'Thailand' & id2_country == 'Italy') | (id2_country == 'Thailand' & id1_country == 'Italy') ~ 'Italy/Thailand'
  )) %>%
  # select out relevant cols
  select(id1, id2, mash_distance, pling_distance, Origin, Country)

# plot clinical samples
country_comparisons %>%
  filter(Origin == 'Clinical') %>%
  ggplot(aes(y = mash_distance, x = Country, fill = Country)) +
           geom_violin() +
           labs(
             title = 'MASH Distance of Klebsiella Plamids in Humans - Clinical',
             y = 'Distance',
             x = ''
           ) +
           scale_fill_manual(values = c('#6290C3', '#C2E7DA', '#F1FFE7'))
  

## Comparing with Stats and Bar Charts ================================

# calculate mean and standard deviation for each category

thai_table <- thailand_only %>%
  group_by(Origin) %>%
  mutate(meanMASH = mean(mash_distance)) %>%
  mutate(sdMASH = sd(mash_distance)) %>%
  mutate(meanPling = mean(pling_distance)) %>%
  mutate(sdPling = sd(pling_distance)) %>%
  select(Origin, meanMASH, sdMASH, meanPling, sdPling) %>%
  distinct() %>%
  mutate(Country = 'Thailand')

italy_table <- italy_only %>%
  group_by(Origin) %>%
  mutate(meanMASH = mean(mash_distance)) %>%
  mutate(sdMASH = sd(mash_distance)) %>%
  mutate(meanPling = mean(pling_distance)) %>%
  mutate(sdPling = sd(pling_distance)) %>%
  select(Origin, meanMASH, sdMASH, meanPling, sdPling) %>%
  distinct() %>%
  mutate(Country = 'Italy')

arg_table <- argentina_only %>%
  group_by(Origin) %>%
  mutate(meanMASH = mean(mash_distance)) %>%
  mutate(sdMASH = sd(mash_distance)) %>%
  mutate(meanPling = mean(pling_distance)) %>%
  mutate(sdPling = sd(pling_distance)) %>%
  select(Origin, meanMASH, sdMASH, meanPling, sdPling) %>%
  distinct() %>%
  mutate(Country = 'Argentina')

# stack together
bar_plot_table <- rbind(thai_table, italy_table, arg_table)

# plot
bar_plot_table %>%
  ggplot(aes(x = Origin, y = meanMASH, fill = Country)) +
  geom_bar(position = 'dodge', stat = 'identity', width = 0.5) +
  scale_fill_manual(values = c('#e76f51', '#e9c46a', '#2a9d8f')) +
  geom_errorbar(aes(ymin = meanMASH - sdMASH, ymax = meanMASH + sdMASH), width = 0.1, position = position_dodge(.5)) +
  theme_minimal() +
  labs(
    title = 'Mean MASH Distance Values by Clinical Origin and Country',
    y = 'mean MASH'
  )

# percentages for write up
country_comparisons %>%
  filter(Country == 'Argentina/Thailand') %>%
  filter(mash_distance < 0.01) %>%
  summarise(
    n = n(),
    prop = n / nrow(country_comparisons %>% filter(Country == 'Argentina/Thailand')) * 100
  )


## repeat plots for pling ======================================

# variation in pling distance across clinical and non-clinical samples by country

thailand_only %>%
  ggplot(aes(y = pling_distance, x = Origin, fill = Origin)) +
  geom_violin() +
  labs(
    title = 'Pling Distance of Klebsiella Plamids in Humans - Thailand',
    y = 'Distance',
    x = ''
  ) +
  scale_fill_manual(values = c('#5AB1BB', '#A5C882', '#F7DD72'))

italy_only %>%
  ggplot(aes(y = pling_distance, x = Origin, fill = Origin)) +
  geom_violin() +
  labs(
    title = 'Pling Distance of Klebsiella Plamids in Humans - Italy',
    y = 'Distance',
    x = ''
  ) +
  scale_fill_manual(values = c('#5AB1BB', '#A5C882', '#F7DD72'))

argentina_only %>%
  ggplot(aes(y = pling_distance, x = Origin, fill = Origin)) +
  geom_violin() +
  labs(
    title = 'Pling Distance of Klebsiella Plamids in Humans - Argentina',
    y = 'Distance',
    x = ''
  ) +
  scale_fill_manual(values = c('#5AB1BB'))

# clinical samples between countries
country_comparisons %>%
  filter(Origin == 'Clinical') %>%
  ggplot(aes(y = pling_distance, x = Country, fill = Country)) +
  geom_violin() +
  labs(
    title = 'Pling Distance of Klebsiella Plamids in Humans - Clinical',
    y = 'Distance',
    x = ''
  ) +
  scale_fill_manual(values = c('#6290C3', '#C2E7DA', '#F1FFE7')) 

# clinical samples between countries barplot
country_table <- country_comparisons %>%
  filter(Origin == 'Clinical') %>%
  group_by(Country) %>%
  mutate(meanMASH = mean(mash_distance)) %>%
  mutate(sdMASH = sd(mash_distance)) %>%
  mutate(meanPling = mean(pling_distance)) %>%
  mutate(sdPling = sd(pling_distance)) %>%
  select(Origin, meanMASH, sdMASH, meanPling, sdPling) %>%
  distinct() 

country_table %>%
  filter(Origin == 'Clinical') %>%
  ggplot(aes(x = Origin, y = meanMASH, fill = Country)) +
  geom_bar(position = 'dodge', stat = 'identity', width = 0.6) +
  scale_fill_manual(values = c('#6290C3', '#C2E7DA', '#F1FFE7')) +
  geom_errorbar(aes(ymin = meanMASH - sdMASH, ymax = meanMASH + sdMASH), width = 0.1, position = position_dodge(.6)) +
  theme_minimal() +
  labs(
    title = 'Comparison of MASH Distance of Clinical Isolates Between Countries',
    y = 'mean MASH'
  )

country_table %>%
  filter(Origin == 'Clinical') %>%
  ggplot(aes(x = Origin, y = meanPling, fill = Country)) +
  geom_bar(position = 'dodge', stat = 'identity', width = 0.6) +
  scale_fill_manual(values = c('#6290C3', '#C2E7DA', '#F1FFE7')) +
  geom_errorbar(aes(ymin = meanPling - sdPling, ymax = meanPling + sdPling), width = 0.1, position = position_dodge(.6)) +
  theme_minimal() +
  labs(
    title = 'Comparison of Pling Distance of Clinical Isolates Between Countries',
    y = 'mean Pling'
  )

# one way anova for stats

country_for_stats <- country_comparisons %>%
  filter(Origin == 'Clinical') %>%
  select(Country, Origin, mash_distance, pling_distance)

# anova - mash
model <- aov(mash_distance ~ Country, data = country_for_stats)

summary(model)

# tukey post-hoc
TukeyHSD(model, conf.level = .95)

# anova - pling
model <- aov(pling_distance ~ Country, data = country_for_stats)

summary(model)

# tukey post-hoc
TukeyHSD(model, conf.level = .95)

# bar plot by pling
bar_plot_table %>%
  ggplot(aes(x = Origin, y = meanPling, fill = Country)) +
  geom_bar(position = 'dodge', stat = 'identity', width = 0.5) +
  scale_fill_manual(values = c('#e76f51', '#e9c46a', '#2a9d8f')) +
  geom_errorbar(aes(ymin = meanPling - sdPling, ymax = meanPling + sdPling), width = 0.1, position = position_dodge(.5)) +
  theme_minimal() +
  labs(
    title = 'Mean Pling Distance Values by Clinical Origin and Country',
    y = 'mean Pling'
  )

# statistical tests between groups ===========================

# mash 

# t tests between clincal and non-clinical only samples in italy and then thailand

italy_clinical <- italy_only %>%
  filter(Origin == 'Clinical')

italy_non_clinical <- italy_only %>%
  filter(Origin == 'Non-Clinical')

# not significant
t.test(italy_clinical$mash_distance, italy_non_clinical$mash_distance)

t.test(italy_clinical$pling_distance, italy_non_clinical$pling_distance)


# thailand 

thai_clinical <- thailand_only %>%
  filter(Origin == 'Clinical')

thai_non_clinical <- thailand_only %>%
  filter(Origin == 'Non-Clinical')

# is actually significant p < 2.2e-16
t.test(thai_clinical$mash_distance, thai_non_clinical$mash_distance)

t.test(thai_clinical$pling_distance, thai_non_clinical$pling_distance)

# pling

# italy - not significant
t.test(italy_clinical$pling_distance, italy_non_clinical$pling_distance)

# thailand - siginficant to same level
t.test(thai_clinical$mash_distance, thai_non_clinical$mash_distance)

## ST307 ========================================================

# look within ST groups for clinically relevant plasmids

# Extract ST307 data from metadata
meta_ST307 <- meta %>%
  filter(hybass_mlst_ST == 'klebsiella_ST307' | hybass_mlst_ST == '307')

# ... and origin labels
origin_labels_ST307 <- meta_ST307 %>%
  select(id, Sample.type.1)

# add on labels via meta data
ST307_hospital_labels <- distances %>%
  filter(id1 %in% meta_ST307$id) %>%
  filter(id2 %in% meta_ST307$id) %>%
  left_join(origin_labels_ST307, by = c('id1' = 'id')) %>%
  left_join(origin_labels_ST307, by = c('id2' = 'id')) %>%
  # filter only for human samples
  filter(grepl('human', Sample.type.1.x) == TRUE) %>%
  filter(grepl('human', Sample.type.1.y) == TRUE) %>%
  # group hospital disease and hospital into the same category in this instance
  mutate(location_s1 = case_when(
    grepl('non-hospital', Sample.type.1.x) == TRUE ~ 'nh',
    .default = 'h'
  )) %>%
  mutate(location_s2 = case_when(
    grepl('non-hospital', Sample.type.1.y) == TRUE ~ 'nh',
    .default = 'h'
  )) %>%
  # now add final label
  mutate(Origin = case_when(
    location_s1 == 'h' & location_s2 == 'h' ~ 'Clinical',
    location_s1 == 'nh' & location_s2 == 'nh' ~ 'Non-Clinical',
    .default = 'Mix'
  )) %>%
  # select out relevant cols
  select(id1, id2, mash_distance, pling_distance, Origin)

# plot hospital vs non-hospital cases in a violin plot
# mash
ST307_hospital_mash_graph <- ST307_hospital_labels %>%
  ggplot(aes(y = mash_distance, x = Origin, fill = Origin)) +
  geom_violin() +
  labs(
    title = 'ST307 MASH Distance of Clinical Origin',
    y = 'MASH Distance',
    x = 'Origin'
  ) +
  scale_fill_manual(values = c('#5AB1BB', '#A5C882', '#F7DD72'))

# pling
ST307_hospital_pling_graph <- ST307_hospital_labels %>%
  ggplot(aes(y = pling_distance, x = Origin, fill = Origin)) +
  geom_violin() +
  labs(
    title = 'ST307 Pling Distance of Clinical Origin',
    y = 'Pling Distance',
    x = 'Origin'
  ) +
  scale_fill_manual(values = c('#5AB1BB', '#A5C882', '#F7DD72'))

# plot together for final figure
ST307_hospital_mash_graph + ST307_hospital_pling_graph + 
  plot_annotation(title = 'ST307 Relatedness Distances within Hospital and Non-Hopsital derived Human Samples',  
                  theme = theme(plot.title = element_text(size = 20, hjust = 0.5))) +
  plot_layout(guides = 'collect')

# calculate plotting table for bar chart
ST307_table <- ST307_hospital_labels %>%
  group_by(Origin) %>%
  mutate(meanMASH = mean(mash_distance)) %>%
  mutate(sdMASH = sd(mash_distance)) %>%
  mutate(meanPling = mean(pling_distance)) %>%
  mutate(sdPling = sd(pling_distance)) %>%
  select(Origin, meanMASH, sdMASH, meanPling, sdPling) %>%
  distinct() %>%
  pivot_longer(!Origin, names_to = 'measure', values_to = 'value') %>%
  rowwise() %>%
  mutate(tool = case_when(
    grepl('MASH', measure) == TRUE ~ 'MASH',
    grepl('Pling', measure) == TRUE ~ 'Pling'
  )) %>%
  mutate(maths = case_when(
    grepl('mean', measure) == TRUE ~ 'mean',
    grepl('sd', measure) == TRUE ~ 'sd'
  )) %>%
  select(!measure) %>%
  pivot_wider(id_cols = c(Origin, tool), names_from = maths, values_from = value)

# plot - skip this mash and pling values dont match
ST307_table %>%
  ggplot(aes(x = Origin, y = mean, fill = tool)) +
  geom_bar(position = 'dodge', stat = 'identity', width = 0.5) +
  scale_fill_manual(values = c('#e76f51', '#e9c46a', '#2a9d8f')) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1, position = position_dodge(.5)) +
  theme_minimal() +
  labs(
    title = 'Mean Distance Values by Clinical Origin in ST307',
    y = 'Mean Distance'
  )

# calculate plotting table for bar chart
ST307_table <- ST307_hospital_labels %>%
  group_by(Origin) %>%
  mutate(meanMASH = mean(mash_distance)) %>%
  mutate(sdMASH = sd(mash_distance)) %>%
  mutate(meanPling = mean(pling_distance)) %>%
  mutate(sdPling = sd(pling_distance)) %>%
  select(Origin, meanMASH, sdMASH, meanPling, sdPling) %>%
  distinct()

# plot
ST307_table %>%
  ggplot(aes(x = Origin, y = meanMASH, fill = Origin)) +
  geom_bar(position = 'dodge', stat = 'identity', width = 0.5) +
  scale_fill_manual(values = c('#5AB1BB', '#A5C882', '#F7DD72')) +
  geom_errorbar(aes(ymin = meanMASH - sdMASH, ymax = meanMASH + sdMASH), width = 0.1, position = position_dodge(.5)) +
  theme_minimal() +
  labs(
    title = 'Mean MASH Distance Values by Clinical Origin in ST307',
    y = 'MASH Distance'
  )

# now pling
ST307_table %>%
  ggplot(aes(x = Origin, y = meanPling, fill = Origin)) +
  geom_bar(position = 'dodge', stat = 'identity', width = 0.5) +
  scale_fill_manual(values = c('#5AB1BB', '#A5C882', '#F7DD72')) +
  geom_errorbar(aes(ymin = meanPling - sdPling, ymax = meanPling + sdPling), width = 0.1, position = position_dodge(.5)) +
  theme_minimal() +
  labs(
    title = 'Mean Pling Distance Values by Clinical Origin in ST307',
    y = 'Pling Distance'
  )

# percentages for write up
ST307_hospital_labels %>%
  filter(Origin == 'Clinical') %>%
  filter(mash_distance < 0.01) %>%
  summarise(
    n = n(),
    prop = n / nrow(ST307_hospital_labels %>% filter(Origin == 'Clinical')) * 100
  )


# ST307 stats --------------------------------

# two sided t tests prep values

ST307_clinical_only <- ST307_hospital_labels %>%
  filter(Origin == 'Clinical')

ST307_non_clinical_only <- ST307_hospital_labels %>%
  filter(Origin == 'Non-Clinical')

# significant 0.05 level
t.test(ST307_clinical_only$mash_distance, ST307_non_clinical_only$mash_distance)

# non signficant
t.test(ST307_clinical_only$pling_distance, ST307_non_clinical_only$pling_distance)


## ST661 ========================================================

# Compare transmission between cows and humans

# Extract ST307 data from metadata
meta_ST661 <- meta %>%
  filter(hybass_mlst_ST == 'klebsiella_ST661')

# ... and origin labels
origin_labels_ST661 <- meta_ST661 %>%
  select(id, Sample.type.1, associated_species)

# most of these samples are from cow, with 2 pig

# label
# add on labels via meta data
ST661_species_labels <- distances %>%
  # filter for ST661 samples
  filter(id1 %in% meta_ST661$id) %>%
  filter(id2 %in% meta_ST661$id) %>%
  # join labels
  left_join(origin_labels_ST661, by = c('id1' = 'id')) %>%
  left_join(origin_labels_ST661, by = c('id2' = 'id')) %>%
  # filter out pig samples for clarity
  filter(!associated_species.x == 'pig') %>%
  filter(!associated_species.y == 'pig') %>%
  # group hospital disease and hospital into the same category in this instance
  mutate(species_s1 = case_when(
    grepl('human', associated_species.x) == TRUE ~ 'human',
    .default = 'cow'
  )) %>%
  mutate(species_s2 = case_when(
    grepl('human', associated_species.y) == TRUE ~ 'human',
    .default = 'cow'
  )) %>%
  # now add final label
  mutate(Origin = case_when(
    species_s1 == 'human' & species_s2 == 'human' ~ 'Human',
    species_s1 == 'cow' & species_s2 == 'cow' ~ 'Cow',
    .default = 'Mix'
  )) %>%
  # select out relevant cols
  select(id1, id2, mash_distance, pling_distance, Origin)

# plot
# mash
ST661_species_mash_graph <- ST661_species_labels %>%
  ggplot(aes(y = mash_distance, x = Origin, fill = Origin)) +
  geom_violin() +
  labs(
    title = 'ST661 MASH Distance of Species Origin',
    y = 'MASH Distance',
    x = ''
  ) +
  scale_fill_manual(values = c('#084C61', '#DB504A', '#E3B505'))

# pling
ST661_species_pling_graph <- ST661_species_labels %>%
  ggplot(aes(y = pling_distance, x = Origin, fill = Origin)) +
  geom_violin() +
  labs(
    title = 'ST661 Pling Distance of Species Origin',
    y = 'Pling Distance',
    x = ''
  ) +
  scale_fill_manual(values = c('#084C61', '#DB504A', '#E3B505'))

# plot together for final figure
ST661_species_mash_graph + ST661_species_pling_graph + 
  plot_annotation(title = 'ST661 Relatedness Distances Between Plasmids Derived from Human and Cow Samples',  
                  theme = theme(plot.title = element_text(size = 20, hjust = 0.5))) +
  plot_layout(guides = 'collect')

# make table for bar chart
ST661_table <- ST661_species_labels %>%
  group_by(Origin) %>%
  mutate(meanMASH = mean(mash_distance)) %>%
  mutate(sdMASH = sd(mash_distance)) %>%
  mutate(meanPling = mean(pling_distance)) %>%
  mutate(sdPling = sd(pling_distance)) %>%
  select(Origin, meanMASH, sdMASH, meanPling, sdPling) %>%
  distinct()

# plot - MASH
ST661_table %>%
  ggplot(aes(x = Origin, y = meanMASH, fill = Origin)) +
  geom_bar(position = 'dodge', stat = 'identity', width = 0.5) +
  scale_fill_manual(values = c('#084C61', '#DB504A', '#E3B505')) +
  geom_errorbar(aes(ymin = meanMASH - sdMASH, ymax = meanMASH + sdMASH), width = 0.1, position = position_dodge(.5)) +
  theme_minimal() +
  labs(
    title = 'Mean MASH Distance Values by Clinical Origin in ST661',
    y = 'MASH Distance'
  )

# plot Pling
ST661_table %>%
  ggplot(aes(x = Origin, y = meanPling, fill = Origin)) +
  geom_bar(position = 'dodge', stat = 'identity', width = 0.5) +
  scale_fill_manual(values = c('#084C61', '#DB504A', '#E3B505')) +
  geom_errorbar(aes(ymin = meanPling - sdPling, ymax = meanPling + sdPling), width = 0.1, position = position_dodge(.5)) +
  theme_minimal() +
  labs(
    title = 'Mean Pling Distance Values by Clinical Origin in ST661',
    y = 'Pling Distance'
  )

# stats ST661 --------------------------------

# sort data for testing
ST661_for_stats <- ST661_species_labels %>%
  select(Origin, mash_distance, pling_distance)

# anova - mash
model <- aov(mash_distance ~ Origin, data = ST661_for_stats)

summary(model)

# tukey post-hoc
TukeyHSD(model, conf.level = .95)

# anova - pling
model <- aov(pling_distance ~ Origin, data = ST661_for_stats)

summary(model)

# tukey post-hoc
TukeyHSD(model, conf.level = .95)

## table associated group on all data, ST307, and ST661 data

meta %>%
  group_by(Sample.type.1) %>%
  summarise(
    n = n()
  )

meta %>%
  group_by(associated_species) %>%
  summarise(
    n = n()
  )

meta_ST307 %>%
  group_by(associated_species) %>%
  summarise(
    n = n()
  )

meta_ST661 %>%
  group_by(associated_species) %>%
  summarise(
    n = n()
  )

# find smallest values for write up

thailand_only %>%
  filter(Origin == 'Non-Clinical') %>%
  arrange(mash_distance) %>%
  head(5)

# how many italy non clinical only
italy_only %>%
  filter(Origin == 'Non-Clinical')

# balance of thailand data
thailand_only %>%
  group_by(Origin) %>%
  summarise(
    n = n()
  )
