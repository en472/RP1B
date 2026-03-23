# RP1B

This repository contains all of the required scripts for the data processing, analysis, and plotting/statistics of the content in Repearch Project 1B. 

## data_formatting.R
A basic R script which filters for Klebsiella specific samples in the meta data, and combines the Pling and MASH distance values for each available pairwise comparison. It also contains a small section of code for transposing the MASH distance matrix from 4315 x 4315 symmetrical matrix into a three column data frame.

## plasmid_diversity_clinical_samples.R
R script for the formatting, labelling, and plotting of the first section of analysis (see report). Also contains the statistical tests performed. Looks at clinical vs. non-clinical klebsiella isolate plasmids within countries, between countries, within ST307, and then human-only vs. cow-only klebsiella isolate plasmids in ST661.

## evaluating_plasmid_analysis_methods.R
R script for comparing the pairwise distances output by MASH and Pling. 

## TetATetR-Sectioning.py
A small python script to read in .fasta sequence files, and cut them at specific places. The cut sites and filenames were input manually each time it was used. 
