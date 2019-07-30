## Tropical Pacific Chlorophyll Algorithm (TPCA): Source code and matchups (Pittman et al., 2019) 

### Contact:

Nicholas Pittman, PhD Candidate: Biogeochemistry & Satellite remote sensing. 
Institute of Marine and Antarctic Studies (IMAS), University of Tasmania, Australia. 
Australian Research Council Center of Excellence for Climate Extremes (CLEX)

email: nic.pittman@utas.edu.au

### Contents

- chl_tpca_algorithms.py
  - Set Python functions to process Rrs data from SeaWiFS, MODIS-Aqua into the TPCA [6] algorithm, built upon CI [4] and OCx [5]. Functions include:
    - blended_chl (Linear blending function [4,6])
    - calculate_chl_ocx (Calculate Chl OCx [4, 5, 6])
    - calculate_chl_ci (Calculate Chl CI [4,6])
    - calculate_seawifs_chl (Calculate TPCA chl for SeaWiFS with Rrs443, Rrs490, Rrs510, Rrs555, Rrs670)
    - calcuate_modis_chl (Calculate TPCA chl for MODIS-Aqua with Rrs443, Rrs488, Rrs547, Rrs667)
    - calculate_meris_chl (Calculate TPCA chl (Default NASA implementation with Rrs443, Rrs490, Rrs510, Rrs560, Rrs665)
- example_seawifs_download.py
  - Example script which uses the requests library to download L3M Daily 2000-01-01 Seawifs wavelengths for Rrs443,490,510,555,670 and the chlor_a file into a new directory: seawifs_data. Cuts the tropical Pacific out of these files, processes the TPCA algorithm and makes 3 plots; TPCA, chlor_a and the difference between the two.
- example_seawifs_matchups.py
  - Example script which produces the TPCA algorithm for SeaWiFS and uses tropical_pacific_matchups/seawifs_matchups.csv to produce chlorophyll estimates for the tropical Pacific and uses chl_statistics to assess model performance.
    - **Note** * The diagnostics produced by this script do not identically reproduce Table 3, SeaWiFS rank 2. The Rrs values provided in the matchup databases are an average of each wavelength in the 45 pixel matchup window. Table 3 was produced instead by calculating chl for each of the 45 pixels, and then averaging the 45 chlorophyll concentrations. This produces slightly different values than seen in the paper. A python Pickle file of the un-averaged Rrs values can be provided on request for accurate reproduction. 
- chl_statistics.py
  - Contains two functions:
    - plot_linear_trend - Function for plotting differences between two chlorophyll variables.
    - check_bias - Function for calculating % wins, bias, and also returns slope, r2 and intercept.  Uses the plot_linear_trend function and prints linear plots to assess model performance.
- requirements.txt 
  - For a conda environment, built using Python 3.7.3. 
- tropical_pacific_matchups/
  - seawifs_matchups.csv
  - modis_matchups.csv
  - meris_matchups.csv

#### Matchup details 

- Satellite to in situ matchup files as used in Pittman et al., 2019 [6] have been openly provided for as per the Journal of Geophysical Research: Oceans guidelines.

- Three files are provided. Each file is a *.csv matchup database for one of the three sensors analysed; SeaWiFS, MODIS-Aqua and MERIS.

- Chlor_a and relevant Rrs data for the three sensors was downloaded from the NASA ocean color portal: https://oceandata.sci.gsfc.nasa.gov/

- Monthly MEI (Multivariate ENSO Index) has matched from: https://www.esrl.noaa.gov/psd/enso/mei.old/

- More details about the matchup process and sensor details / sources. R,eprocessing versions are 2018.0 for SeaWiFS and MODIS-Aqua, and 2012.1 for MERIS.

- In situ fluorometric data has been compiled from three unique sources [1,2,3]
- Each sensor has a different number of matchups due to time period in orbit, different orbits, cloud cover and quality control.
- These matchup files are those used in [6]; 2 day radius and 1 pixel. This results in a total of 5 days (day, day, observation, day, day) and a grid of 9 pixels, with the observation located in the centre pixel. The matchups provided are the mean of these 45 pixels. For example in space:

| Sat  |  Sat   | Sat  |
| ---- | :----: | ---- |
| Sat  | **Ob** | Sat  |
| Sat  |  Sat   | Sat  |

### Matchup database

The *.csv files contain 28 fields including:

| Field              | Contains / Format                                            | Decimal Places |
| ------------------ | ------------------------------------------------------------ | -------------- |
| obs_date           | Observation date (YYYY-MM-DD)                                | NaN            |
| obs_lat            | Observation Latitude (decimal °)                             | 3              |
| obs_lon            | Observation Longitude (decimal °; processed into 0-360 rather than -180 to 180 as the dateline becomes problematic in the Pacific) | 3              |
| obs_source         | Source of in situ observation (See references [1,2,3])       | NaN            |
| chl_type           | Chlorophyll type (Fluoroescence or HPLC)                     | NaN            |
| in_situ_chl        | Observed Chlorophyll in mg m<sup>-3</sup>. An average of the surface 20m (first optical depth) where possible. If more than 1 observation occurred in a day and 0.1°, those observations were averaged. However, if an HPLC observation was in this 1 day 0.1° pixel, it was used in preference over fluorescence data as per O'Reilly et al., 1998. | 4              |
| NASA_chlor_a       | Chlor_a product (blended OCI estimates from Chl<sub>OCx </sub> and Chl<sub>CI</sub>) downloaded from the NASA ocean colour portal, an average of the matchup window. | 4              |
| TPCA_chl           | Tropical Pacific Chlorophyll Algorithm chlorophyll estimates. A modified version of the NASA_Chlor_a product, with updated OCx coefficients and CI to OCx blending window using the method described in  Pittman et al., (2019) [6]. Derived from the rrs fields. | 4              |
| chl_ci             | NASA Level 3 derived chl<sub>CI</sub> from the Hu et al., 2012 Color Index method, using CI. | 5              |
| chl_ocx            | NASA Level 3 derived chl<sub>OCx</sub> from the Hu et al., 2012 Color Index method. | 5              |
| CI                 | CI derived from  the Hu et al., 2012 method [4].             | 5              |
| MBR                | Max Band Ratio (MBR) calculated from max(blue/green) [5].    | 5              |
| max_blue_rrs       | Max of the (443nm, 490nm and 510nm) wavelengths. MODIS-Aqua excludes 510nm due to no wavelength here. | 5              |
| rrs443 / 443 / 443 | Rrs nearest to 443nm (Blue) for all three sensors (SeaWiFS, MODIS-Aqua, MERIS) | 5              |
| rrs490 / 488 / 490 | Rrs nearest to 490nm (Blue) for all three sensors (SeaWiFS, MODIS-Aqua, MERIS) | 5              |
| rrs510 / NaN / 510 | Rrs nearest to 510nm (Blue) for SeaWiFS and MERIS. MODIS-Aqua does not use this wavelength. | 5              |
| rs555 / 547 / 560  | Rrs nearest to 555nm (Green) for SeaWiFS and MERIS. MODIS-Aqua uses 547nm rather than 555nm [4]. | 5              |
| rrs670 / 667 / 665 | Rrs nearest to 670 (Red) for all three sensors (SeaWiFS, MODIS-Aqua, MERIS) | 5              |
| MEI                | MEI (Multivariate ENSO Index). Obtained from: https://www.esrl.noaa.gov/psd/enso/mei/. MEI >= 1 is defined as El Nino and <= -1 La Nina, and neutral between. | 3              |
| day_radius         | Day radius. Files provided are all ±2, a total matchup window of 5 days with the observation day in the middle) | 0              |
| pixel_radius       | Pixel radius. Files provided are all 1, a total matchup area of 9 pixels, with he observation in the centre pixel. A radius of approximately ± 15km) | 0              |
| sat_start_date     | Satellite start date; first day of the matchup window        | NaN            |
| sat_end_date       | Satellite end date; last day of the matchup window           | NaN            |
| sat_start_lat      | Satellite lowest latitude in matchup box.                    | 3              |
| sat_end_lat        | Satellite highest longitude in matchup box.                  | 3              |
| sat_start_lon      | Satellite lowest longitude in matchup area.                  | 3              |
| sat_end_lon        | Satellite highest longitude in matchup area.                 | 3              |
| validation_set     | The database was split in half randomly to prevent over-fitting. 0 indicates the training data, 1 is for testing / validation. | 0              |



#### References

1. Boyer, T.P., Baranova, O.K., Coleman, C., Garcia, H.E., Grodsky, A., Locarnini, R.A., Mishonov, A.V., Paver, C.R., Reagan, J.R., Seidov, D., et al. World Ocean Database 2018. A. V. Mishonov, Technical Editor, NOAA Atlas NESDIS 87.
2. Strutton, P.G., Evans, W., and Chavez, F.P. (2008). Equatorial Pacific chemical and biological variability, 1997–2003. Global Biogeochemical Cycles 22.
3. Valente, A., Sathyendranath, S., Brotas, V., Groom, S., Grant, M., Taberner, M., Antoine, D., Arnone, R., Balch, W.M., Barker, K., et al. (2016). A compilation of global bio-optical in situ data for ocean-colour satellite applications. Earth System Science Data 18.
4. Hu, C., Lee, Z., and Franz, B. (2012). Chlorophyll a algorithms for oligotrophic oceans: A novel approach based on three-band reflectance difference. Journal of Geophysical Research: Oceans *117*.
5. O’Reilly, J.E., Maritorena, S., Mitchell, B.G., Siegel, D.A., Carder, K.L., Garver, S.A., Kahru, M., and McClain, C. (1998). Ocean color chlorophyll algorithms for SeaWiFS. Journal of Geophysical Research: Oceans *103*, 24937–24953.
4. Pittman, N. A., Strutton, P.G., Johnson, R., Matear R., Chavez, F.P. (2019). An assessment and improvement of tropical Pacific Ocean Color algorithms. Submitted to Journal of Geophysical Research: Oceans 
