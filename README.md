
# Comparing Bayesian hierarchical models of maritime pine sapling height growth in a common garden experiment that combines genomic, phenotypic and climatic data.

This repository contains the scripts needed to reproduce the analyses in Archambeau et al. 2021 *Combining climatic and genomic data improves
range-wide tree height growth prediction in a forest tree*

Preprint avaibale here: https://www.biorxiv.org/content/10.1101/2020.11.13.382515v2


**Paper abstract**

Population response functions based on climatic and phenotypic data from common gardens have long been the gold standard for predicting quantitative trait variation in new environments. However, prediction accuracy might be enhanced by incorporating genomic information that captures the neutral and adaptive processes behind intra-population genetic variation. We used five clonal common gardens containing 34 provenances (523 genotypes) of maritime pine (*Pinus pinaster* Aiton) to determine whether models combining climatic and genomic data capture the underlying drivers of height-growth variation, and thus improve predictions at large geographical scales. The plastic component explained most of the height-growth variation, probably resulting from population responses to multiple environmental factors. The genetic component stemmed mainly from climate adaptation, and the distinct demographic and selective histories of the different maritime pine gene pools. Models combining climate-of-origin and gene pool of the provenances, and positive-effect height-associated alleles (PEAs) captured most of the genetic component of height-growth and better predicted new provenances compared to the climate-based population response functions. Regionally-selected PEAs were better predictors than globally-selected PEAs, showing high predictive ability in some environments, even when included alone in the models. These results are therefore promising for the future use of genome-based prediction of quantitative traits.


**Analyses**

Below, you can find the `html` corresponding to the `.md` files in the folder `reports/`:

1. In the folder `ExplanatoryVariables`: 

    - [Climatic similarity matrices](https://juliettearchambeau.github.io/HeightPinpinClonapin/ClimSimMatrices.html)
    - [Calculating counts of gPEAs and rPEAs](https://juliettearchambeau.github.io/HeightPinpinClonapin/CalculatingGPEAandRPEA.html)
    - [Calculating the genomic relationship matrice with all gene pools](https://juliettearchambeau.github.io/HeightPinpinClonapin/GenomicRelationshipMatrices.html)
    - [Gene pool-specific genomic relationship matrices](https://juliettearchambeau.github.io/HeightPinpinClonapin/GenePoolSpecificKinshipMatrix.html)
    
2. In the folder `ExploringData`:

    - [Details on the experimental design](https://juliettearchambeau.github.io/HeightPinpinClonapin/DetailsExperimentalDesign.html)
    - [Exploratory analyses](https://juliettearchambeau.github.io/HeightPinpinClonapin/ExploratoryAnalyses.html)

3. In the folder `ModelOutputsAnalysis`:

    - [Model Posterior Distributions - P1 partition](https://juliettearchambeau.github.io/HeightPinpinClonapin/ModelPosteriorDistributionsP1.html)
    - [Model posterior distributions - P2 partition](https://juliettearchambeau.github.io/HeightPinpinClonapin/ModelPosteriorDistributionsP2.html)
    - [Model posterior distributions - P3 partition](https://juliettearchambeau.github.io/HeightPinpinClonapin/ModelPosteriorDistributionsP3.html)
    - [Posterior Predictive Checks and Residuals - P2 partition](https://juliettearchambeau.github.io/HeightPinpinClonapin/PosteriorPredictiveChecks.html)
    - [Predictive Models - Figures](https://juliettearchambeau.github.io/HeightPinpinClonapin/PosteriorPredictiveFigures.html)
    - [Broad-sense heritability estimates](https://juliettearchambeau.github.io/HeightPinpinClonapin/EstimatingHeritability.html)
    