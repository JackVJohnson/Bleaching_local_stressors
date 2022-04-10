# Bleaching_local_stressors

data and codes to run all the analysis requried for the paper "Local anthropogenic stress does not exacerbate coral bleaching under global climate change"

Codes run in the following order. Skip to step 6 for the analysis. The previous steps are data collection. 

1.	RC_with_Ecoregions.R
- Input Reef check data (supplied from reefcheck.org) and combine with ERG shape files from Veron et al (2015) cited in the main paper. 
- Output is the RC_ERG.csv where reef check data are combined with coral ecoregions.

2.	RC_ERG_CoRTAD.R
- Input RC_ERG.csv with CoRTAD data from https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:NCEI-CoRTADv6. The code extracts sea surfcae tmeperature (SST) values from CoRTAD which are used in the analyses.
- Output is the RC_ERG_SST_DHW.csv.

3.	RC_ERG_CoRTAD_SEDAC.R
- Input RC_ERG_SST_DHW.csv with sedac human pop https://sedac.ciesin.columbia.edu/data/collection/gpw-v4 to extract the human population density within 25, 50, and 100km of each reef check survey.
- Output is the RC_ERG_SST_DHW_Pop.csv

4.	RC_ERG_CoRTAD_SEDAC_MPA.R
- Input RC_ERG_SST_DHW_Pop.csv with world mpa shape file https://www.iucn.org/theme/protected-areas/our-work/world-database-protected-areas which determines whether a coral reef surveyed by reef check were inside the jurisdiction of an MPA.
- Output is the RC_ERG_SST_DHW_Pop_MPA.csv

5.	RC_ERG_CoRTAD_SEDAC_MPA_WP.R 
- Input RC_ERG_SST_DHW_Pop_MPA.csv with world pop tiff files https://www.worldpop.org/ for another measure of human population (from a different source) at 25, 50, and 100km distance from each reef check survey.
- Output is Bleaching_data_new_pop.csv

6.	Global_analysis.R
- Input the Bleaching_data_new_pop.csv to run the global analysis for Bayesian models and produce coefficient plot

7.	Ecoregion_analysis.R
- Input Bleaching_data_new_pop.csv to run ecoregion analyses and obrain ecoregional p-values 

8.	Ecoregion_maps.R
- Input the ecoregion_p_values.csv and Ecoregion_shape_files to produce the ecoregion maps in the paper. 

9.	Survey_maps.R
- Input Bleaching_data_new_pop.csv and World_shape_files to produce the distirbution and density of reef check surveys throughout the world used in this study. 
