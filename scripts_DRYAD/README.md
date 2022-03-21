# Data and code for the paper: 'Combining climatic and genomic data improves range-wide tree height growth prediction in a forest tree'

**Authors:** Juliette Archambeau$^1$, Marta Benito Garzón$^1$, Frédéric Barraquand$^2$, Marina de Miguel Vega$^{1,3}$, Christophe Plomion$^1$ and Santiago C. González-Martínez$^1$

**1** INRAE, Univ. Bordeaux, BIOGECO, F-33610 Cestas, France

**2** CNRS, Institute of Mathematics of Bordeaux, F-33400 Talence, France

**3** EGFV, Univ. Bordeaux, Bordeaux Sciences Agro, INRAE, ISVV, F-33882, Villenave d’Ornon, France


**Published in *The American Naturalist* ** (manuscript # 60477R2)

## Data in the DRYAD repository

### Raw data

#### Height and climatic data 

**In the dataset `HeightClimateSoilData_33121obs_25variables.csv`.**


The 33,121 height observations comes from 5 common gardens of the CLONAPIN network, consisting of 34 populations (i.e. provenances) and 523 clones (i.e. genotypes). 

Annual and monthly climatic data were extracted from the [EuMedClim database at 1-km resolution](https://gentree.data.inra.fr/climate/). The scripts for the extraction can be found here: [monthly data for climate in the test sites](https://github.com/JulietteArchambeau/HeightPinpinClonapin/blob/master/scripts/DataMunging/01_ExtractingMonthlyClimateEuMedClimData.R) and [annual data for climate in the provenances](https://github.com/JulietteArchambeau/HeightPinpinClonapin/blob/master/scripts/DataMunging/01_ExtractingAnnualClimateEuMedClimData.R). We calculated six variables that describe extreme and average temperature and precipitation conditions in the test sites during the year preceding the measurements ([script available here](https://github.com/JulietteArchambeau/HeightPinpinClonapin/blob/master/scripts/DataMunging/02_CalculatingSiteMonthlyClimaticVariables.R)) and four variables that describe the mean temperature and precipitation in the provenance locations over the period from 1901 to 2009, representing the climate under which provenances have evolved ([script available here](https://github.com/JulietteArchambeau/HeightPinpinClonapin/blob/master/scripts/DataMunging/02_CalculatingProvAnnualClimaticVariables.R)).

Missing data are indicated with *NA*. 

Meaning of the columns:

  - `obs`: unique code differentiating each observation (i.e. each height measurement)
  - `tree`: tree identity
  - `site`: test site (common garden)
  - `clon`: clone (i.e. genotype)
  - `prov`: provenance (i.e. population)
  - `latitude_site`: latitude of the test site (in degrees)
  - `longitude_site`: longitude of the test site (in degrees)
  - `latitude_prov`: latitude of the provenance (in degrees)
  - `longitude_prov`: longitude of the provenance (in degrees)
  - `age`: tree age when the measurement was taken (in months)
  - `height`: tree height (in mm)
  - `pre_summer_min_site` (*min.presummer* in the manuscript): minimum of the monthly precipitation during summer -June to September- in the test sites during the year preceding the measurement (°C)
  - `pre_mean_1y_site` (*mean.pre* in the manuscript): mean of the monthly precipitation in the test sites during the year preceding the measurement (mm)
  - `tmn_min_1y_site` (*min.tmn* in the manuscript):  minimum of the monthly minimum temperatures in the test sites during the year preceding the measurement (°C)
  - `tmx_max_1y_site` (*max.tmx* in the manuscript): maximum of the monthly maximum temperatures in the test sites during the year preceding the measurement (°C)
  - `pre_max_1y_site` (*max.pre* in the manuscript): maximum of the monthly precipitation in the test sites during the year preceding the measurement (mm)
  - `tmx_mean_1y_site` (*mean.tmax* in the manuscript): mean of monthly maximum temperatures in the test sites during the year preceding the measurement (°C)
  - `bio1_prov` (*mean.temp* in the manuscript): average of the annual daily mean temperature in the provenances over the period from 1901 to 2009 (°C).
  - `bio5_prov` (*max.temp* in the manuscript): average of the maximum temperature of the warmest month in the provenances over the period from 1901 to 2009 (°C).
  - `bio12_prov` (*min.pre* in the manuscript): average of the precipitation of the driest month in the provenances over the period from 1901 to 2009 (mm).
  - `bio14_prov` (*mean.pre* in the manuscript): average of the annual precipitation in the provenances over the period from 1901 to 2009 (mm).
  - `Q1`: proportion of assignment to the northern African (NA) gene pool for each clone.
  - `Q2`: proportion of assignment to the Corsican (C) gene pool for each clone.
  - `Q3`: proportion of assignment to the central Spain (CS) gene pool for each clone.
  - `Q4`: proportion of assignment to the French Atlantic (FA) gene pool for each clone.
  - `Q5`: proportion of assignment to the Iberian Atlantic (IA) gene pool for each clone.
  - `Q6`: proportion of assignment to the south-eastern Spain (SES) gene pool for each clone.
  - `max.Qvalue`: proportion of assignment to the main gene pool for each clone.
  - `max.Q`: main gene pool for each clone.
  - `P1`: qualitative variable indicating the assignment of each observation to the test or train dataset of the P1 partition.
  - `P2`: qualitative variable indicating the assignment of each observation to the test or train dataset of the P2 partition.
  - `P3`: qualitative variable indicating the assignment of each observation to the test or train dataset of the P3 partition.
  


#### Genomic data 

**In the dataset `GenomicData_5165SNPs_523clones.csv`.**

This file contains the genotype (noted as 0, 1 or 2) of each clone. There are 5,165 SNPs in rows and 523 clones in columns. Missing data are indicated with *NA*. 


#### piMASS outputs from de [Miguel et al. (2022)](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.16367?casa_token=1nNTc88Iy40AAAAA%3ALd4EOK5ehk_cEHIkw5A9l8nk0NPzUzlYPX8eAjVCikIjHP0WJ1kxoHJSZjMLFsZcP-8wdbNuNrlOfp1jzw)

**Datasets `height_all_sites_res.mcmc.txt`, `height_iberian_atlantic_res.mcmc.txt`, `height_french_atlantic_res.mcmc.txt` and `height_mediterranean_res.mcmc.txt`.**

These files correspond to the piMASS outputs of the Bayesian variable selection regression [VSR] implemented in the piMASS software ([Guan & Stephens 2011](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-5/issue-3/Bayesian-variable-selection-regression-for-genome-wide-association-studies-and/10.1214/11-AOAS455.full)) and which allows the identification of SNPs associated with the phenotypes. Here, the phenotypes are BLUPs for height estimated in:

  - the five CLONAPIN common gardens for the file `height_all_sites_res.mcmc.txt`.
  - the Iberian Atlantic common gardens (Asturias and Portugal) for the file  `height_iberian_atlantic_res.mcmc.txt`.
  - the French Atlantic common garden (Bordeaux) for the file `height_french_atlantic_res.mcmc.txt`.
  - the Mediterranean common gardens (Madrid and Cáceres) for the file `height_mediterranean_res.mcmc.txt`.

Accroding the [piMASS manual](https://www.haplotype.org/download/pimass-manual.pdf), the output files contain:

  - `rs`: SNP ID.
  - `chr`: chromosome (no information on it in Miguel et al. 2022, so there are only '0').
  - `pos`: position.
  - `postc`: estimates of the posterior inclusion probabilities based on simple counting.
  - `postrb`: estimates of the posterior inclusion probabilities based on Rao-Blackwellization.
  - `beta`: the naive estimates of the posterior effect size.
  - `betarb`: **Rao-Blackwellized estimates of the posterior effect size, that are used in the present study to differentiate height-associated SNPs from SNPs not associated with height (i.e. 'neutral' SNPs).**

### Intermediate files

#### Gene pool-specific GRMs

**Files `GRM_AX.csv`, with *X* being the gene pool number.**

The gene pool-specific genomic relationship matrices (GRM) are calculated in the script `2_CalculateGenePoolSpecificGRM.R` and are then used when fitting *model M5*.

## Scripts in Zenodo

Available in the following *zenodo* repository:

The code included in the present *zenodo* repository was run on *R version 3.6.3* and *RStudio version 1.1.463* and constitutes the code necessary to replicate the analyses of the present Am Nat paper.

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


An extended version of the code can be found in the github repository [https://github.com/JulietteArchambeau/HeightPinpinClonapin](https://github.com/JulietteArchambeau/HeightPinpinClonapin), where data sorting from a larger phenotypic database (inlcuding all traits measured in the common gardens of the [CLONAPIN network](https://www6.bordeaux-aquitaine.inrae.fr/biogeco/Ressources/In-situ-dispositifs-de-terrain-observation-experimentation/Tests-de-provenances/CLONAPIN)), climate and soil data extraction,additional exploratory analyses, model output analyses and visualizations are included.