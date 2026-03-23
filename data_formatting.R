
# Data Processing ================================================

## Libraries:

library(tidyverse)

## Meta data ==================================================

meta <- read.csv('4315-metadata.csv')

# filter for italian, thai, and argentinian data only by study ID
meta <- meta %>%
  filter(grepl('SPARK', meta$id) == TRUE | 
           grepl('OH-DART', meta$id) == TRUE |
           grepl('KLAR', meta$id) == TRUE) %>%
  # filter only for klebsiella species
  filter(hybass_mlst_scheme == 'klebsiella')

# get counts of different klebsiella species
meta %>%
  group_by(hybridass_kleborate3_species) %>%
  summarise(
    n = n()
  )

# select specific columns relevant for this analysis
meta <- meta %>%
  select(id, Sample.type.1, associated_species, hybass_mlst_ST, pling_type)

# check for duplicate ID entries - no, 2547 unique ids and 2547 rows in df
length(unique(meta$id))

# check for NA values - yes, but these are in associated species column - not critical
anyNA(meta)

# fill in missing Argentina sample origin data - all were clinical human samples

meta$Sample.type.1[grepl('KLAR', meta$id) == TRUE] <- 'human - hospital'

meta$associated_species[grepl('KLAR', meta$id) == TRUE] <- 'human'


# output --------------------------------------------------------

write.csv(meta, 'meta_data.csv')

## MASH distance matrix =====================================================

mash <- read.table('4315_mash_matrix.tab', header = TRUE, row.names = 1)

# change . to - in column headers ----------------------------------

# extract column names
mash_colnames <- colnames(mash)

# loop through and replace text
for(i in 1:length(mash_colnames)){
  str <- mash_colnames[i]
  new_str <- str_replace_all(str, '\\.', '-')
  mash_colnames[i] <- new_str
}

# now replace old headers
colnames(mash) <- mash_colnames

# transpose the matrix ---------------------------------------------

# make empty results data frame
mash_transposed <- NA

# begin iterating through columns 
for(col in 1:ncol(mash)){
  
  # make a temporary data frame based on the column
  tmp <- data.frame(
    sample1 = rownames(mash)[col:ncol(mash)],
    sample2 = colnames(mash)[col],
    distance = mash[[col]][col:ncol(mash)] 
  )
  
  # add to results
  mash_transposed <- rbind(mash_transposed, tmp)
  
} 

# clean up ----------------------------------------------------------

# remove top row (empty)
mash_transposed <- mash_transposed[-1, ]

# now remove self-self comparisons
mash_transposed <- mash_transposed %>%
  rowwise() %>%
  filter(sample1 != sample2)

# filter for italy, thailand, and argentina ------------------------------------

mash_filtered <- mash_transposed %>%
  filter(grepl('SPARK', sample1) == TRUE | grepl('OH-DART', sample1) == TRUE | grepl('KLAR', sample1) == TRUE) %>%
  filter(grepl('SPARK', sample2) == TRUE | grepl('OH-DART', sample2) == TRUE | grepl('KLAR', sample2) == TRUE)

## Pling distance matrix =========================================

pling <- read_tsv('4315_pling_distance_matrix.tsv')

# remove self-self comparisons from pling distance matrix
pling <- pling %>%
  rowwise() %>%
  filter(plasmid_1 != plasmid_2)

# filter for italy, thailand, and argentina
pling_filtered <- pling %>% 
  filter(grepl('SPARK', plasmid_1) == TRUE | grepl('OH-DART', plasmid_1) == TRUE | grepl('KLAR', plasmid_1) == TRUE) %>%
  filter(grepl('SPARK', plasmid_2) == TRUE | grepl('OH-DART', plasmid_2) == TRUE | grepl('KLAR', plasmid_2) == TRUE)

## Combine mash and pling distances =================================

# change colnames to be clearer
colnames(mash) <- c('sample1', 'sample2', 'mash_distance')

colnames(pling) <- c('sample1', 'sample2', 'pling_distance')

# re-order names across the rows so that they match
mash <- mash %>%
  rowwise() %>%
  mutate(
    id1 = min(sample1, sample2),
    id2 = max(sample1, sample2)
  ) %>%
  ungroup() %>%
  select(id1, id2, mash_distance)

pling <- pling %>%
  rowwise() %>%
  mutate(
    id1 = min(sample1, sample2),
    id2 = max(sample1, sample2)
  ) %>%
  ungroup() %>%
  select(id1, id2, pling_distance)

# join mash and pling distances into one data frame
distances <- left_join(pling, mash, by = c('id1', 'id2'))


# output ------------------------------------------------------

write.csv(distances, 'distances.csv')
