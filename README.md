# Ichthyophonus_timeseries
Data and code related to manuscript titled: Elevated prevalence of Ichthyopohnus sp. in a collapsed population of Pacific herring and submitted to ICES J of Marine Science 2021

Data: 
1) AK Herring disease samples for R: From yearly disease surveillance conducted in Sitka Sound and Prince William Sound, led by Paul Hershberger at the USGS Marrowstone Marine Field Station

2) Historical measurements 2009-2019.csv: Measurements of ventricle area and schizont areas using stitched images of histological sections of herring hearts. Hearts were collected during yearly disease surveillance in Sitka Sound and Prince William Sound

3) Mature_prop_at_age_PWS_Sitka.csv: Estimates from ADFG of the proportion of mature fish in each age class for each year between 2007 and 2019. 

4) PWS_Recruitment.csv: Estimates of the recruitment biomass (thousand metric tons) and mature biomass (thousand metric tons) from PWS 2019 stock assessment. 

6) Sitka Recruitment.csv: Estimates of the recruitment biomass (thousand metric tons) and mature biomass (thousand metric tons) from Sitka 2020 forecast model. 

7) Summarized PDO NPGO.csv: Estimates of NPGO and PDO for Oct-March and April-September. April through September estimates are lagged by a year (i.e., for the previous year), because the response variables from the disease surveys were taken in late March and early April.

8) SummarizedUpwelling.csv: Mean upwelling indices for PWS and Sitka Sound for March-October and April-September

9) alk_Cordova_fishmethods.txt: Age length keys for Prince William Sound for the years 2007-2019, calculated using the fishmethods package (see attached R code). 

10) alk_Sitka_fishmethods.txt: Age length keys for Sitka Sound for the years 2007-2019, calculated using the fishmethods package (see attached R code). 

Code:
1) Age-Length analysis Cordova_fishmethods.R: Code to create and validate an age-length key for PWS Pacific herring using age-length data from PWS. Note that code contains a link to access online data

2) Age-Length analysis Sitka_fishmethods.R: Code to create and validate an age-length key for Sitka Sound Pacific herring using age-length data from PWS.

3) Age comparison disease surveillance data.R: Compares scale annuli ages of Pacific herring collected for disease surveillance in PWS with those predicted by the PWS Age-length key using Bland-Altman methods

4) Population Prevalence.R: Code to generate estimates of population-level Ichthyophonus prevalence in Sitka Sound and PWS from 2007-2019. Also calculates age-length relationships for each year class over time. 

5) GLMMs.R: Code to conduct generalized linear mixed models to identify environmental and demographic correlates of Ichthyophonus sp. infection prevalence and severity


