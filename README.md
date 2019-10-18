# Codes and Data for DGRP lifespan analyses

### Data

#### Individual level data

* life.all.csv: Individual level lifespan data for all individuals across environments, the columns (comma separated) are: line,sex,vial,lifespan,block,temperature,vial-ID
* life.lnS.all.csv: Individual (pseudo) level lifespan microenvironmental variability data across environments, the columns (comma separated) are: line,sex,vial(pseudo),lifespan ln(sigma),block,temperature,replicate(pseudo)-ID

#### Data used for GWAS

* *.pheno: the data used for GWAS (http://dgrp2.gnets.ncsu.edu): the file name contains sex.temperature.phenotype (mean or lns) with .pheno as the extension. When there are two temperatures separated by a period, the phenotype provided is the difference between the two environments used for mapping plasticity QTLs. Within each comma separated text file, the first column is line ID and the second column is phenotype.

#### Data used for xQTL mapping in the AIP

* fastq.gz: Please download the data from SRA (https://www.ncbi.nlm.nih.gov/sra) with the accession number for the BioProject PRJNA577841.

#### Data for the RNAi tests

* 18C: Formatted_18Degree_RNAi.csv
* 25C: Formatted_25Degree_RNAi.csv
* 28C: Formatted_28Degree_RNAi.csv

### Codes

#### Code for quantitative genetic analyses (ANOVAs)

* fullMixed.sas: the variance components are estimated using Proc Mixed. This code is for full model. Within each environment, this SAS program is modified to the respective reduced model.
* fullGLM.sas: the signficance of the terms are tested using Proc GLM. This code is for full model. Within each environment, this SAS program is modified to the respective reduced model.

#### Code for GWAS in the DGRP

* We use the DGRP server (http://dgrp2.gnets.ncsu.edu) for GWAS, including annotations.

#### Code for xQTL mapping in the AIP

The xQTL mapping consists of several steps:

1) mapping: bwa.bash

2) post-alignment processing using GATK

3) allele counts and filtering: filterVariants.pl and pileupCov.pl

4) statistical tests: bsaReps.R, bsaRepsDiff.R, bsaRepsSexDiff.R for the tests within each environment, between environments but within the same sex, and between sexes within the same environment respectively.
