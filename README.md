# DRDNet (Disease risk-associated pseudo-dynamic networks)
A statistical framework for recovering dynamic networks from static data with applications to tissue-specific gene network reconstruction from GTEx data

# simulations
In the simulations folder, we provide the codes to generate the results for three cases in the main text of the paper: 1. when the index(agent) is fixed and evenly spaced; 2. when the index(agent) is known random variable; 3. when the index(agent) is unobserved and imputed(estimated). For each case, we consider four combinations of sample sizes (n=100, n=200) and variances of measurement errors (0.1 and 0.5). They are named: n100v1_one, n100v5_one, n200v1_one, n200v5_one, espectively. All of the cases are under the generation of covariate followed by Bernoulli distribution. We also provide the codes to take into account the case where the covariate is from Normal distribution, named n100v1one_nor, n100v5one_nor, n200v1one_nor, n200v5one_nor.

# real data applications
In the real data application folder, we provide two code files: one, named data_process_seperate, is to process the data and calculate the imputed index(agent); the other one, named cvm_code_seperate, is to run the model and make network plots. The datasets used for the analyses described in this manuscript were obtained from dbGaP at http://www.ncbi.nlm.nih.gov/gap through dbGaP accession number phs000424.v7.p2.

