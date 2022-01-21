# Bleaching_local_stressors

data and codes to run all the analysis requried for the paper "Local anthropogenic stress does not exacerbate coral bleaching under global climate change"

Codes run in the following order. Skip to step 6 for the analysis. The previous steps are data collection. 

1.	RC_with_Ecoregions.R
-	Input Reef check data and combine with ERG shape files from Veron et al (2015)

2.	RC_ERG_CoRTAD.R
-	RC_ERG.csv with CoRTAD data from https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:NCEI-CoRTADv6

3.	RC_ERG_CoRTAD_SEDAC.R
-	RC_ERG_SST_DHW.csv with sedac human pop https://sedac.ciesin.columbia.edu/data/collection/gpw-v4 

4.	RC_ERG_CoRTAD_SEDAC_MPA.R
-	RC_ERG_SST_DHW_Pop.csv with world mpa shape file https://www.iucn.org/theme/protected-areas/our-work/world-database-protected-areas 

5.	RC_ERG_CoRTAD_SEDAC_MPA_WP.R 
-	RC_ERG_SST_DHW_Pop_MPA.csv with world pop tiff files https://www.worldpop.org/ 

6.	Global_analysis.R
-	Bleaching_data_new_pop.csv

7.	Ecoregion_analysis.R
-	Bleaching_data_new_pop.csv

8.	Ecoregion_maps.R
-	ecoregion_p_values.csv

9.	Survey_maps.R
-	Bleaching_data_new_pop.csv
