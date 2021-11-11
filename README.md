# README

This repository contains the dataflow and scripts used to process the COI metabarcode reads in the paper Robinson et al., 2021 (submitted to BioRxiv).

## Infiles

diatom_invert_combined.csv contains the ESV ids, sample names, and taxonomic assignments.  

Sites.csv

metadata_diatom.csv

metadata_invert.csv

original_globi.csv

original_globi2.csv

## R Scripts

Fig1_Richness.R calcultes ESV richness and effective number of ESVs.  Uses diatom_invert_combined.csv as an infile.  Creates Fig1_Richness.pdf

Fig2_NMDS_FigS5_Dissimilarities.R creates NMDS plots with fitted environmental parameters.  Uses diatom_invert_combined.csv, Sites.csv, metadata_diatom.csv, and metadata_invert.csv as infiles.  Creates Scree_diat.pdf, Scree_macroinvert.pdf, stressplot_diat.pdf, stressplot_macroinvert.pdf, BetaDispersion_diatom.pdf, BetaDispersion_macroinvertebrate.pdf, Fig2_NMDS.jpeg, values for Table S4, and FigS5_dissimilarities.jpeg.

Fig3_part1_searchGlobi.R to query GloBI and reformat for cheddar.  Uses diatom_invert_combined.csv, original_globi.csv and original_globi2.csv as infiles.  Creates globi.csv from the first search of genus & species queries, globi2.csv from the second search where species were collapsed into genus queries, a directory called 'waterloo' that contains site-specific results that are needed for running cheddar in the next script.  

Fig3_part2_makeNetworks_Fig4_circle_FigS5_centrality.R creates trophic networks.  Uses infiles produced from Fig3_part1_searchGlobi.R script.  Calculates values for Table 2, Fig3_foodwebs.pdf, Fig4_circlegraphs.pdf, data for Table S5-S9, Fig5_centrality.pdf

FigS1_Map.R creates a site map.  Creates Sites.csv and FigS1_Map.jpg .

FigS2_Rarefaction.R shows rarefied sampling curves.  Uses diatom_invert_combined.csv as an infile.  Creates FigS2_Rarefaction.jpg .

FigS3_FigS4_Heatmaps.R plots rarefaction curves.  Uses cat.csv as an infile.  Creates FigS3_Diatom_Heatmap.jpg and FigS4_Macroinvertebrate_Heatmap.jpg .

FigS6_IndicSpecies.R conducts indicator value and correlation analyses.  Uses diatom_invert_combined.csv and Sites.csv as infiles.  Creates FigS6_Indicators.jpg .

## References

Robinson, C.V., Porter, T.M., Maitland, V.C., Wright, M.T.G., Hajibabaei, M. (2021) Multi-marker metabarcoding resolves subtle variations in freshwater condition: Bioindicators, ecological traits, and trophic interactions.  Submitted to BioRxiv.

## Acknowledgements

I would like to acknowledge funding from the Canadian government from the Genomics Research and Development Initiative (GRDI), Metagenomics-Based Ecosystem Biomonitoring (Ecobiomics) project.

Last updated: November 10, 2021
