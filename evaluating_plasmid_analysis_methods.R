
# Evaluation of MASH vs. Pling for classifying plasmids =========

# look within ST307 (abundant in sample and relevant clinically)

# libraries
library(tidyverse)
library(ggplot2)

## Read in data ========================================================

meta <- read.csv('meta_data.csv')

distances <- read.csv('distances.csv')

# drop first columns for formatting

meta <- meta[, 2:ncol(meta)]

distances <- distances[, 2:ncol(distances)]

# Filter for ST307 samples ===========================================

meta_ST307 <- meta %>%
  filter(hybass_mlst_ST == 'klebsiella_ST307' | hybass_mlst_ST == '307')

distances_ST307 <- distances %>%
  filter(id1 %in% meta_ST307$id) %>%
  filter(id2 %in% meta_ST307$id)

# plot mash vs pling - group where tools disagree =================

distances_ST307 %>%
  ggplot(aes(x = pling_distance, y = mash_distance)) +
  geom_point() +
  labs(
    title = 'MASH vs. Pling Distance in ST307 data',
    y = 'MASH Distance',
    x = 'Pling Distance'
  )

# overlay with hub status (according to pling) ========================================

# Label pairwise samples as hub or not

# extract out labels from metadata
hub_labels_ST307 <- meta_ST307 %>%
  select(id, pling_type)

ST307_hub_labelled <- distances_ST307 %>%
  left_join(hub_labels_ST307, by = c('id1' = 'id')) %>%
  left_join(hub_labels_ST307, by = c('id2' = 'id')) %>%
  # label according to hub/ non-hub
  mutate(hub_s1 = case_when(
    pling_type.x == 'hub' ~ 'hub',
    .default = 'non-hub'
  )) %>%
  mutate(hub_s2 = case_when(
    pling_type.y == 'hub' ~ 'hub',
    .default = 'non-hub'
  )) %>%
  # now add final label
  mutate(hub = case_when(
    hub_s1 == 'hub' & hub_s2 == 'hub' ~ 'Both Hubs',
    hub_s1 == 'non-hub' & hub_s2 == 'non-hub' ~ 'Neither Hubs',
    .default = 'Mix'
  )) %>%
  # select out relevant cols
  select(id1, id2, mash_distance, pling_distance, hub)

# replot but highlight hub plasmids - association with hub plasmids

ST307_hub_labelled %>%
  ggplot(aes(x = pling_distance, y = mash_distance, colour = hub)) +
  geom_point() +
  labs(
    title = 'MASH vs. Pling Distance in ST307 data',
    y = 'MASH Distance',
    x = 'Pling Distance',
    colour = 'Hub Presence'
  )

## Calculate how many are in each hub category ==========================

# across all samples
hub_proportion_all_samples <- ST307_hub_labelled %>%
  group_by(hub) %>%
  summarise(
    all = n(),
    all_prop = all / nrow(ST307_hub_labelled) * 100
  ) %>%
  select(!all)

# find number of samples in the cluster
cluster_sample_no <- ST307_hub_labelled %>%
  filter(mash_distance > 0.05 & pling_distance < 10) %>%
  nrow()

# across clustered samples
hub_proportion_cluster_samples <- ST307_hub_labelled %>%
  filter(mash_distance > 0.05 & pling_distance < 10) %>%
  group_by(hub) %>%
  summarise(
    cluster = n(),
    cluster_prop = cluster / cluster_sample_no * 100
  ) %>%
  select(!cluster)

# join together - how much of each group is within the cluster
hub_proportion_ST307 <- hub_proportion_all_samples %>%
  left_join(hub_proportion_cluster_samples, by = 'hub') %>%
  # calculate percentage decrease
  rowwise() %>%
  mutate(change = cluster_prop / all_prop)

# calculate by name, overall
ST307_id_list <- stack(ST307_hub_labelled[c('id1', 'id2')]) %>%
  unique() %>%
  left_join(meta_ST307, by = c('values'= 'id')) %>%
  select(values, pling_type)

ST307_id_list %>%
  filter(pling_type == 'hub') %>%
  summarise(
    n = n(),
    prop = n / nrow(ST307_id_list) * 100
  )

# calculate by name, in cluster
ST307_id_list_cluster <- ST307_hub_labelled %>%
  filter(mash_distance >= 0.05 & pling_distance < 10) %>%
  select(id1, id2)

ST307_id_list_cluster_names <- stack(ST307_id_list_cluster[c('id1', 'id2')]) %>%
  unique() %>%
  left_join(meta_ST307, by = c('values' = 'id')) %>%
  select(values, pling_type)

ST307_id_list_cluster_names %>%
  filter(pling_type == 'hub') %>%
  summarise(
    n = n(),
    prop = n / nrow(ST307_id_list_cluster_names) * 100
  )

## Identify which hub genes are most represented in the cluster ======

# select out cluster samples and store in separate data frame

hub_cluster_samples_ST307 <- ST307_hub_labelled %>%
  filter(mash_distance > 0.05 & pling_distance < 10)

# extract all of the names present in both pairwise columns

highly_represented_hub_samples <- stack(hub_cluster_samples_ST307[c('id1', 'id2')]) %>%
  # join with meta data to identify hub plasmids
  left_join(hub_labels_ST307, by = c('values' = 'id')) %>%
  # filter only for hub plasmids
  filter(pling_type == 'hub') %>%
  select(!ind)

# table as a proportion of the total number of plasmids in the cluster

hub_proportion_in_cluster_table <- highly_represented_hub_samples %>%
  group_by(values) %>%
  summarise(
    n = n(),
    prop. = n / nrow(highly_represented_hub_samples)
  ) %>%
  arrange(desc(prop.)) %>%
  head(5)

# plot hub proportion in cluster as a graph

ggplot(hub_proportion_in_cluster_table, aes(x = reorder(values, -n), y = prop.)) +
  geom_bar(stat = 'identity', fill = '#5AB1BB') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = round(prop., 3)), vjust = -0.5) +
  ylab('Proportion of appearance in pairwise comparisons') +
  xlab('') +
  labs(
    title = 'Overrepresentation of ST307 Hub Plasmids in MASH/Pling Outlier Cluster'
  ) +
  ylim(0, 0.75) 


# Select out largest disagreements for most common hub plasmids =====
# for comparison in ziplign

# select out spark1218 relations
SPARK1218_samples <- ST307_mash_pling %>%
  filter(id1 == 'SPARK_1218_C1_5'| id2 == 'SPARK_1218_C1_5') %>%
  # order by 'difference' - pling = 1 and mash ordered by highest distance
  #filter(pling_distance == '1') %>%
  arrange(desc(mash_distance)) %>%
  head(10)

# export results table
write.csv(SPARK1218_samples, 'SPARK1218_samples.csv')

# OH-DART 109360
OHDART109360_samples <- ST307_mash_pling %>%
  filter(id1 == 'OH-DART_109360-KN1_3'| id2 == 'OH-DART_109360-KN1_3') %>%
  # order by 'difference' - pling = 1 and mash ordered by highest distance
  #filter(pling_distance == '1') %>%
  arrange(desc(mash_distance)) %>%
  head(10)

# export
write.csv(OHDART109360_samples, 'OHDART109360_samples.csv')

# SPARK 1205
SPARK1205_samples <- ST307_mash_pling %>%
  filter(id1 == 'SPARK_1205_C1_4'| id2 == 'SPARK_1205_C1_4') %>%
  # order by 'difference' - pling = 1 and mash ordered by highest distance
  filter(pling_distance == '2') %>%
  arrange(desc(mash_distance)) %>%
  head(10)

# export
write.csv(SPARK1205_samples, 'SPARK1205_samples.csv')

# SPARK 1205 but no limit on pling
SPARK1205_samples_noplinglimit <- ST307_mash_pling %>%
  filter(id1 == 'SPARK_1205_C1_4'| id2 == 'SPARK_1205_C1_4') %>%
  # order by 'difference' - pling = 1 and mash ordered by highest distance
  #filter(pling_distance == '2') %>%
  arrange(desc(mash_distance)) %>%
  head(10)

# export
write.csv(SPARK1205_samples_noplinglimit, 'SPARK1205_samples_noplinglimit.csv')
  

