# Data and code for the paper: 'Combining climatic and genomic data improves range-wide tree height growth prediction in a forest tree'

**Authors:** Juliette Archambeau<sup>1</sup>, Marta Benito Garzón<sup>1</sup>, Frédéric Barraquand<sup>2</sup>, Marina de Miguel<sup>1,3</sup>, Christophe Plomion<sup>1</sup> and Santiago C. González-Martínez<sup>1</sup>

**1** INRAE, Univ. Bordeaux, BIOGECO, F-33610 Cestas, France

**2** CNRS, Institute of Mathematics of Bordeaux, F-33400 Talence, France

**3** EGFV, Univ. Bordeaux, Bordeaux Sciences Agro, INRAE, ISVV, F-33882, Villenave d’Ornon, France

**Corresponding author:** Juliette Archambeau, juli.archambeau@orange.fr

Published in *The American Naturalist* (manuscript # 60477R2), [paper available here](https://www.journals.uchicago.edu/doi/abs/10.1086/720619?journalCode=an).

## Data in the DRYAD repository

### Raw data

#### Height, climate and population structure

**In the dataset `HeightClimateSoilData_33121obs_25variables.csv`.**


The 33,121 height observations comes from 5 common gardens of the CLONAPIN network, consisting of 34 populations (i.e. provenances) and 523 clones (i.e. genotypes). 

Annual and monthly climatic data were extracted from the [EuMedClim database at 1-km resolution](https://gentree.data.inra.fr/climate/). The scripts for the extraction can be found here: [monthly data for climate in the test sites](https://github.com/JulietteArchambeau/HeightPinpinClonapin/blob/master/scripts/DataMunging/01_ExtractingMonthlyClimateEuMedClimData.R) and [annual data for climate in the provenances](https://github.com/JulietteArchambeau/HeightPinpinClonapin/blob/master/scripts/DataMunging/01_ExtractingAnnualClimateEuMedClimData.R). We calculated six variables that describe extreme and average temperature and precipitation conditions in the test sites during the year preceding the measurements ([script available here](https://github.com/JulietteArchambeau/HeightPinpinClonapin/blob/master/scripts/DataMunging/02_CalculatingSiteMonthlyClimaticVariables.R)) and four variables that describe the mean temperature and precipitation in the provenance locations over the period from 1901 to 2009, representing the climate under which provenances have evolved ([script available here](https://github.com/JulietteArchambeau/HeightPinpinClonapin/blob/master/scripts/DataMunging/02_CalculatingProvAnnualClimaticVariables.R)).

Missing data are indicated with *NA*. 

Meaning of the columns:

  1. `obs`: unique code differentiating each observation (i.e. each height measurement)
  2. `tree`: tree identity
  3. `site`: test site (common garden)
  4. `clon`: clone (i.e. genotype)
  5. `prov`: provenance (i.e. population)
  6. `latitude_site`: latitude of the test site (in degrees)
  7. `longitude_site`: longitude of the test site (in degrees)
  8. `latitude_prov`: latitude of the provenance (in degrees)
  9. `longitude_prov`: longitude of the provenance (in degrees)
  10. `age`: tree age when the measurement was taken (in months)
  11. `height`: tree height (in mm)
  12. `pre_summer_min_site` (*min.presummer* in the manuscript): minimum of the monthly precipitation during summer -June to September- in the test sites during the year preceding the measurement (°C)
  13. `pre_mean_1y_site` (*mean.pre* in the manuscript): mean of the monthly precipitation in the test sites during the year preceding the measurement (mm)
  14. `tmn_min_1y_site` (*min.tmn* in the manuscript):  minimum of the monthly minimum temperatures in the test sites during the year preceding the measurement (°C)
  15. `tmx_max_1y_site` (*max.tmx* in the manuscript): maximum of the monthly maximum temperatures in the test sites during the year preceding the measurement (°C)
  16. `pre_max_1y_site` (*max.pre* in the manuscript): maximum of the monthly precipitation in the test sites during the year preceding the measurement (mm)
  17. `tmx_mean_1y_site` (*mean.tmax* in the manuscript): mean of monthly maximum temperatures in the test sites during the year preceding the measurement (°C)
  18. `bio1_prov` (*mean.temp* in the manuscript): average of the annual daily mean temperature in the provenances over the period from 1901 to 2009 (°C).
  19. `bio5_prov` (*max.temp* in the manuscript): average of the maximum temperature of the warmest month in the provenances over the period from 1901 to 2009 (°C).
  20. `bio12_prov` (*min.pre* in the manuscript): average of the precipitation of the driest month in the provenances over the period from 1901 to 2009 (mm).
  21. `bio14_prov` (*mean.pre* in the manuscript): average of the annual precipitation in the provenances over the period from 1901 to 2009 (mm).
  22. `Q1`: proportion of assignment to the northern African (NA) gene pool for each clone.
  23. `Q2`: proportion of assignment to the Corsican (C) gene pool for each clone.
  24. `Q3`: proportion of assignment to the central Spain (CS) gene pool for each clone.
  25. `Q4`: proportion of assignment to the French Atlantic (FA) gene pool for each clone.
  26. `Q5`: proportion of assignment to the Iberian Atlantic (IA) gene pool for each clone.
  27. `Q6`: proportion of assignment to the south-eastern Spain (SES) gene pool for each clone.
  28. `max.Qvalue`: proportion of assignment to the main gene pool for each clone.
  29. `max.Q`: main gene pool for each clone.
  30. `P1`: qualitative variable indicating the assignment of each observation to the test or train dataset of the P1 partition.
  31. `P2`: qualitative variable indicating the assignment of each observation to the test or train dataset of the P2 partition.
  32. `P3`: qualitative variable indicating the assignment of each observation to the test or train dataset of the P3 partition.
  


#### Genomic data 

**In the dataset `GenomicData_5165SNPs_523clones.csv`.**

This file contains the genotype (noted as 0, 1 or 2) of each clone. 

5,165 SNPs in rows and 523 clones in columns. 

Missing data are indicated with *NA*. 


#### piMASS outputs from [de Miguel et al. (2022)](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.16367?casa_token=1nNTc88Iy40AAAAA%3ALd4EOK5ehk_cEHIkw5A9l8nk0NPzUzlYPX8eAjVCikIjHP0WJ1kxoHJSZjMLFsZcP-8wdbNuNrlOfp1jzw)

**Datasets `height_all_sites_res.mcmc.txt`, `height_iberian_atlantic_res.mcmc.txt`, `height_french_atlantic_res.mcmc.txt` and `height_mediterranean_res.mcmc.txt`.**

These files correspond to the piMASS outputs of the Bayesian variable selection regression (BVSR) implemented in the piMASS software ([Guan & Stephens 2011](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-5/issue-3/Bayesian-variable-selection-regression-for-genome-wide-association-studies-and/10.1214/11-AOAS455.full)) and which allows the identification of SNPs associated with the phenotypes. Here, the phenotypes are BLUPs for height estimated in:

  - the five CLONAPIN common gardens for the file `height_all_sites_res.mcmc.txt`.
  - the Iberian Atlantic common gardens (Asturias and Portugal) for the file  `height_iberian_atlantic_res.mcmc.txt`.
  - the French Atlantic common garden (Bordeaux) for the file `height_french_atlantic_res.mcmc.txt`.
  - the Mediterranean common gardens (Madrid and Cáceres) for the file `height_mediterranean_res.mcmc.txt`.

According to the [piMASS manual](https://www.haplotype.org/download/pimass-manual.pdf), the output files contain:

  - `rs`: SNP ID.
  - `chr`: chromosome (no information on it in Miguel et al. 2022, so there are only '0').
  - `pos`: position.
  - `postc`: estimates of the posterior inclusion probabilities based on simple counting.
  - `postrb`: estimates of the posterior inclusion probabilities based on Rao-Blackwellization.
  - `beta`: the naive estimates of the posterior effect size.
  - `betarb`: **Rao-Blackwellized estimates of the posterior effect size, that are used in the present study to differentiate height-associated SNPs from SNPs not associated with height (i.e. 'neutral' SNPs).**

### Intermediate files

#### Climatic similarity matrices

**File `VarCovMatSites.csv`**: covariance matrix describing the climatic similarity among test sites during the year preceding the measurements (Figure S6). 

**File `VarCovMatProvenancesP1.csv`**: covariance matrix describing the climatic similarity among provenances over the period 1901 to 2009 (Figure S9).

The climatic similarity matrices are calculated in the script `3_CreateClimaticSimilarityMatrices.R` and are then used in models *M3* to *M6*.

#### Gene pool-specific GRMs

**Files `GRM_AX.csv`**

with *X* being the gene pool ID:

  - 1 for Northern Africa
  - 2 for Corsica
  - 3 for Central Spain
  - 4 for the Atlantic part of France
  - 5 for the Atlantic part of the Iberian peninsula
  - 6 for South-eastern Spain


The gene pool-specific genomic relationship matrices (GRM) are calculated in the script `4_CalculateGenePoolSpecificGRM.R` and are then used when fitting *model M5*.

#### Counts of height-associated positive-effect alleles (PEAs)

**File `CountPEAs.csv`**

Meaning of the columns:

  1. `clon`: Clone (i.e. genotype) ID.
  2. `count_all_350`: PEA counts based on the BVSR run over the five common gardens (i.e. see file `height_all_sites_res.mcmc.txt`).
  3. `count_fratl_350`: PEA counts based on the BVSR run over the Iberian Atlantic common gardens (Asturias and Portugal) (i.e. see `height_iberian_atlantic_res.mcmc.txt`).
  4. `count_ibatl_350`: PEA counts based on the BVSR run over in the French Atlantic common garden (Bordeaux) (i.e. see `height_french_atlantic_res.mcmc.txt`).
  5. `count_med_350`:  PEA counts based on the BVSR run over the the Mediterranean common gardens (Madrid and Cáceres)(i.e. see `height_mediterranean_res.mcmc.txt`).
  

The count of global and regional PEAs (gPEAs and rPEAs, respectively) are calculated in the script `5_CalculateCountsPEAs.R` and are then used in models *M7*, *M8*, *M11* and *M12*.

#### Test datasets of the three partitions

**Files `TestP1prepared.csv`, `TestP2prepared.csv` and `TestP3prepared.csv`:** test datasets of the P1, P2 and P3 partitions, respectively.

Meaning of the columns:

  - *Columns 1 to 21.* same columns as in the file `HeightClimateSoilData_33121obs_25variables.csv`.
  - *Columns 22 to 27.* same columns as the columns `Q1`, `Q2`, `Q3`, `Q4`, `Q5` and  `Q6` of the file `HeightClimateSoilData_33121obs_25variables.csv`.
  - *Columns 28 to 32.* same columns as in the file `HeightClimateSoilData_33121obs_25variables.csv`.
  - *Column 33.* `age.sc`: tree age scaled based on the variance and mean of the train dataset of the corresponding partition.
  - *Columns 34 to 45.* `Q1` to `Q6` and `clon1` to `clon6`: columns required to estimate one varying intercept for each gene pool based on the gene pool-specific GRMs in the *brms* R package. 
  - *Column 46.* `site_age`: column for the varying intercept associated with the covariance matrix describing the climatic similarity among test sites during the year preceding the measurements.
  - *Column 47.* `prov_clim`: column for the varying intercept associated with the covariance matrix describing the climatic similarity among provenances over the period 1901 to 2009.
  - *Columns 48 to 51.* same columns as in the file `CountPEAs.csv`.
  - *Column 52.* `rPEA`: counts of regional PEAs.
  - *Column 53.* `rPEA.sc`: counts of regional PEAs scaled based on the variance and mean of the train dataset of the corresponding partition.
  - *Column 54.* `gPEA.sc`: counts of global PEAs scaled based on the variance and mean of the train dataset of the corresponding partition.
  - *Columns 55 to 61.* climatic variables (describing the climate at the location of the test sites and provenances) scaled based on the variance and mean of the train dataset of the corresponding partition.

The test datasets of the three partitions are formatted in the script `9_PrepareTestDatasetsForModelEvalutation.R`.

## Scripts in Zenodo

The code included in the present *zenodo* repository was run on *R version 3.6.3* and *RStudio version 1.1.463* and constitutes the code necessary to replicate the analyses of the present study.

Here is what the different scripts are for:

  - `1_DetailsExperimentalDesign.R`: to create the tables S1, S2 and S3 describing the experimental design of the study.
  - `2_ExploreCorrelationsClimaticVariables.R`: to build the figures S4, S5, S7 and S8 describing the correlations among the variables related to the climatic conditions in the common gardens and at the location of the provenances and the variables related to the population genetic structure (the proportion of gene pool assignement for each clone).
  - `3_CreateClimaticSimilarityMatrices.R`: to create the covariance matrices describing the climatic similarity among test sites during the year preceding the measurements (shown in Figure S6 and used in models *M3* to *M6*) and the climatic similarity among provenances (shown in Figure S9 and used in model *M6*).
  - `4_CalculateGenePoolSpecificGRM.R`: to create the gene pool-specific genomic relationship matrices (GRMs) that were then used to calculate gene pool-specific total genetic variance and broad-sense heritabilities in model *M6*.
  - `5_CalculateCountsPEAs.R`: to calculate the counts of global and regional height-associated positive-effect alleles (gPEAs and rPEAs, respectively) based on the GWAS outputs from [Miguel et al. (2022)](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.16367?casa_token=1nNTc88Iy40AAAAA%3ALd4EOK5ehk_cEHIkw5A9l8nk0NPzUzlYPX8eAjVCikIjHP0WJ1kxoHJSZjMLFsZcP-8wdbNuNrlOfp1jzw).
  - `6_FitBayesianModelsP1Partition.R`: to fit the hierarchical Bayesian height-growth models on the P1 partition with the *brms* R package. 
  - `7_FitBayesianModelsP2Partition.R`: to fit the hierarchical Bayesian height-growth models on the P2 partition with the *brms* R package. 
  - `8_FitBayesianModelsP3Partition.R`: to fit the hierarchical Bayesian height-growth models on the P3 partition with the *brms* R package. 
  - `9_PrepareTestDatasetsForModelEvalutation.R`: to format the test datasets which are then used to evaluate the model predictive ability.
  - `10_EvaluateModelPerformance.R`: to evaluate the goodness-of-fit and predictive ability of the models.
  - `11_CreateFigure2.R`: to build Figure 2 related to the genetic and plastic bases of height-growth variation and their potential underlying drivers.
  - `12_CreateFigures4andS10.R`: to build Figures 4 and S10 aimed at visualizing the out-of-sample and in-sample model predictive ability, respectively.
  - `13_CreateTablesFiguresPosteriorDistributions.R`: to create the tables and figures of the model posterior distributions found in Section 6 of the Supplementary Information.
  - `14_ExtracteGenePoolH2_TableS27S28andFigureS13.R`: to extract from model *M5* the gene pool-specific total genetic variances and heritabilities.
  - `15_CreateFigureS14.R`: to build Figure S14 showing the posterior distributions of the provenance, gene pool and provenance climate intercepts from models *M4*, *M5* and *M6*.
  - `16_CreateFiguresS15toS20.R`: to build Figures S15 to S20 showing the posterior distribution of gene pool intercepts and the site-specific intercepts of the gPEAs, rPEAs, the minimum precipitation during the driest month and the maximum temperature of the warmest month in models *M7* to *M12*.
  - `17_QstFstAnalysis.R`: to perform the $Q_{ST}-F_{ST}$ analysis mentioned in Section 5.2 of the Discussion.
  - `18_CheckDistributionHeightCaceresMadrid_FigS21.R`: to build Figure S21 aimed at checking the height distributions in Cáceres and Madrid (see Section 5.3 of the Discussion). 


An extended (but much more messy!) version of the code can be found in the github repository [https://github.com/JulietteArchambeau/HeightPinpinClonapin](https://github.com/JulietteArchambeau/HeightPinpinClonapin), where data sorting from a larger phenotypic database (inlcuding all traits measured in the common gardens of the [CLONAPIN network](https://www6.bordeaux-aquitaine.inrae.fr/biogeco/Ressources/In-situ-dispositifs-de-terrain-observation-experimentation/Tests-de-provenances/CLONAPIN)), climate and soil data extraction and some additional analyses and visualizations are included.