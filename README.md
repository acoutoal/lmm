# lmm_eigen
Implements several fast linear mixed model approaches for correlated samples / related individuals.
The methods are based on the Generalised Least Squares approach.

# gen_sim
Simulates genotype and phenotype data subject to population structure with k subpopulations. The genotype differences in allele frequencies can be specificied for each sub-population. Phenotypes consists of ancestry-related random effects, genetic effects, and random error 

# gwas_lmm_demo
A demo of the above functions using a simulated GWAS with sub-population structure. We compare different methods for correcting population structure and assess the impact on the effect size estimates of these methods. 

# LMM Lite
We also implement the method described in the following paper:
FaST linear mixed models for genome-wide association studies
C Lippert, J Listgarten, Y Liu, CM Kadie, RI Davidson… - Nature …, 2011 - nature.com
