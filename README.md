# Gen Score Pipeline

Note: run `org.genvisis.cnv.gwas.utils.GeneScorePipeline -h` to see what is required to run the pipeline.

## Prior to Analysis

Create these files:


## Results

| Result Column | Explanation |
| ------------- | ------------- |
| STUDY | study name  |
| DATAFILE | Data file used for results line |
| INDEX-THRESHOLD |  |
| FACTOR | .pheno filename used for results line |
| BASE-R-SQUARED | r2 for model ran using genetics only no covariates |
| R-SQR | r2 for model ran using genetics with covariates |
| R-DIFF | difference between the two r2; larger is better |
| P-VALUE | from software |
| EXCEL-SIG | p-value calculated in excel |
| BETA | effect estimate for the risk score |
| SE | standard error of the effect estimate |
| NUM | total number used in analysis |
| CASES | total number of cases used in the analysis (if continuous pheno will be NA) |
| CONTROLS | total number of controls used in the analysis (if continuous pheno will be NA) |
| PAIRED-T-P-VALUE | used with TRIOS; average of parent compared to child (N/A if not a trio) |
| WILCOXON-SIGNED-RANK-P-VALUE | used with TRIOS; average of parent compared to child (N/A if not a trio) |
| PAIRED-STAT-NUM-TRIOS | used with TRIOS; average of parent compared to child (0 if not a trio) |
| #sigInMeta | <5x10-08 (number of variants in the metafile that were below the 5x10^-8 threshold [or whatever threshold you specified in the command]) |
| #indexVariantsInMeta | Number we have in our dataset |
| #indexVariantsInDataset | Number we have in our dataset |
| B-F-SCORE | How much is captured if T < U : missingness; ideal =1 |
| INVCHI-SCORE | How much is captured if T < U : missingness; ideal =1 |
