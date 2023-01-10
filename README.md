# UrbanARGs

This repository contains all analysis scripts for the manuscript "Metagenome-assembled genomes indicate that antimicrobial resistance genes are highly prevalent among urban bacteria and multidrug and glycopeptide resistances are ubiquitous in most taxa" (submitted to Frontiers in Microbiology on 06.09.2022, accepted on 02.01.2023).

Authors: Stefanía Magnúsdóttir, Joao Pedro Saraiva, Alexander Bartholomäus, Majid Soheili, Rodolfo Brizola Toscan, Junya Zhang, Ulisses Nunes da Rocha


# Descriptions
## 01_mudoger_gotupick.sh
The submission script used to cluster the MAGs to genomic operational taxonomic units (gOTUs) using the gOTU picking script from Module 5 in MuDoGeR.

## 02_completeness_contamination.R
Plots the completeness against contamination (from CheckM) of all MAGs that clustered to gOTU species. Summarizes medium and high quality information.

## 03_count_args_per_mag.R
Counts the number of unique ARGs that were identifed per gOTU species according to DeepARG.

## 04_count_gotus_per_arg_class.R
Counts the number of gOTU species that have at least one ARG in each ARG class.

## 05_arg_class_prevalence_phylum.R
Calculates the prevalence of each ARG class at the phylum level based on ARGs found in the gOTU species.

## 06_gotus_args_per_arg_class.R
Calculates the number of gOTU species that have each individual ARG and plots the results as boxplots per ARG class.

## 07_nmi_calculations.R
Calculates the normalized mutual information (NMI) and chi-squared test p-value for every ARG class and taxon combination. Reports back those combinations that were identified as non-random (NMI >=0.5, p-value <0.05)

## 08_count_args_per_mag_identity50.R
Same as 03_count_args_per_mag.R, but uses only DeepARG results with percent identity >=50.

## 09_arg_class_prevalence_phylum_identity50.R
Same as 05_arg_class_prevalence_phylum.R, but uses only DeepARG results with percent identity >=50.

## 10_gotus_args_per_arg_class_identity50.R
Same as 06_gotus_args_per_arg_class.R, but uses only DeepARG results with percent identity >=50.

## 11_count_gotus_per_arg_class_identity50.R
Same as 04_count_gotus_per_arg_class.R, but uses only DeepARG results with percent identity >=50.

## 12_arg_class_prevalence_table.R
Calculates the prevalence of each ARG class per gOTU species. Used for the phylogenetic tree heatmap.

## 13_gotu_tree.R
Prepares the GTDB identifier intput for phyloT to generate the phylogenetic tree. Reads in the tree file and plots the tree with the prevalence heatmap from 
12_arg_class_prevalence_table.R and the reported non-random taxons from 07_nmi_calculations.R.

## 15_plasmids_risk_virulence.R
Takes in the PlasFlow and virulence factor alignment data in addition to the DeepARG output to analyze the occurrence of ARGs on plasmids vs. chromosomes and calculate the number of virulence factors per gOTU (species).
