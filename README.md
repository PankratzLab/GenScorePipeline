# Gen Score Pipeline

Note: run `org.genvisis.cnv.gwas.utils.GeneScorePipeline -h` to see what is required to run the pipeline.


## Prior to analysis

1. Create files:
    - Genetic dataset
    - data.txt file that points to genetic dataset (tab-delimited file – D <tab> label <tab> location of data/map files <tab> extension of data files <tab> extension of map files <tab> full path to IDs file [any file; the first two columns are read as IDs])

    - Effect files saved as .meta files 
      - Meta files must contain ‘MarkerName’, ‘Effect’, ‘Freq’, ‘Pvalue’, ‘Chr’ and ‘Position’ (MarkerName should be reflective of name in the genetic dataset) (column names are not fixed - see appendix for notes on usable column names).
      - If ‘Chr’ and ‘Position’ are not present, there is a pre-process step available to look those up using DBSnp databases. To obtain chr and position for a list of rsIDs go to http://genreport.umn.edu/#topOfPage
      - For a single SNP the effect (beta value) should be 1.

    - .pheno files 
      - The first 3 columns should be ‘fid’, ‘iid’, and ‘case_control’ 
      - Then add columns for any covariates (e.g., age, sex, etc.)

2. Create a location directory for GSP (e.g., GSP)
    - Make a Data directory inside the GSP directory
      - Place .pheno files and a data.txt file here
    - Place .meta files inside the GSP location directory

  
## Run GenScorePipeline 
  
Run GSP from the directory above the GSP directory.
- To run interactively: <br> 
  module load R/3.5.0 <br>
  jcp org.genvisis.cnv.gwas.utils.GeneScorePipeline workDir=/ rLibsDir=/

- To submit as job (.qsub) saved above GSP directory: <br>
  echo "start GSP at: " `date` <br>
  cd DIRECTORYNAME <br>
  module load R/3.5.0 <br>
  java -Djava.awt.headless=true -Xmx16G -jar plab-internal.jar org.genvisis.cnv.gwas.utils.GeneScorePipeline workDir=./GSP rLibsDir=<DIRECT PATH TO  PERSONAL R LIBRARY> <br>
  echo "end GSP at: " `date`

  
## Review your results

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
                                              
                                              
## Appendix 1: Column Names
	
Column names can be any from a list of allowed aliases, e.g., ‘SNP’ instead of ‘MarkerName’.
Column aliases are defined in the Genvisis source code in the org.pankratzlab.common.Aliases class, using:
  MARKER_NAMES
  CHRS
  POSITIONS
  EFFECTS
  ALLELE_FREQS
  PVALUES
  STD_ERRS
