# subsetting_cancer_imputed_genotype_bialletic_dosage
This script is intended to subset TCGA cancer imputed genotype bialletic dosage files 

## Prerequisite
- [R 3.0+](https://www.r-project.org)
- [data.table](https://github.com/Rdatatable/data.table): install.packages('data.table')
- [dplyr](https://github.com/tidyverse/dplyr): install.packages('dplyr')
- [tidyr](http://tidyr.tidyverse.org): install.packages('tidyr)
- [bigrquery](https://github.com/r-dbi/bigrquery): devtools::install_github("rstats-db/bigrquery")
- [google cloud account](https://cloud.google.com/) 


## Run
```
		Rscript subsetting_tcga_cancer_bialletic_dosage.R

```

## outputs 
- {TCGA project name}_filtered_biallelic_chr1.csv.gz