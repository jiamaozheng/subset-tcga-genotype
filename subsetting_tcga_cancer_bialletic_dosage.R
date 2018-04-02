# Author: Jiamao Zheng (jiamaoz@yahoo.com) 04/02/2018 

# latest stable version
# install.packages("bigrquery")
# install.packages('devtools')
# devtools::install_github("rstats-db/bigrquery")

library(data.table)
library(tidyr)
library(dplyr)
library(bigrquery)

# google cloud project ID 
project <- "cancer-predictdb"
cat('google cloud project ID: ', project, '\n\n')

# TCGA projects 
cat('querying all TCGA project names from BigQuery RNAseq data \n\n')
# https://bigquery.cloud.google.com/table/isb-cgc:TCGA_hg19_data_v0.RNAseq_Gene_Expression_UNC_RSEM
all_projects_sql <- "SELECT UNIQUE(project_short_name) FROM [isb-cgc:TCGA_hg19_data_v0.RNAseq_Gene_Expression_UNC_RSEM]"
all_projects = query_exec(all_projects_sql, project = project, use_legacy_sql = T)
colnames(all_projects) = c('project_name')

# table used to hold cancer sample size 
summary_table <- data.frame()

# subset dosage files 
for (i in 1:22) {
  cat('Chromosome ', i, '\n')
  
  chr = fread(paste('zcat < /Volumes/im-lab/nas40t2/jiamao/Projects/cancer_PredictDB/TCGA/TCGA/data/hg19/genotype/genotype/dosages/chr', i, '.filtered.biallelic.txt.gz', sep = ''))
  varIDs = chr$varID
  row.names(chr) = varIDs
  chr = chr %>% select(-1)

  for (project_name in all_projects) {
    cat('working on ', project_name, '\n')
    
    # individual ids 
    id_sql <- paste("SELECT UNIQUE(case_barcode) FROM [isb-cgc:TCGA_hg19_data_v0.RNAseq_Gene_Expression_UNC_RSEM] WHERE project_short_name = '", project_name, "'", sep = '')
    ids =  query_exec(id_sql, project = project, useLegacySql = T)
    colnames(ids) = c('id')
  
    aa = ''
    k = 1 
    for (a in ids$id) {
      for (b in colnames(chr)) {
        if (substr(b, 1, 12)==a) {
            aa[k]=b 
            k = k + 1 
        }
      }
    }
    
    chr = chr %>% select(one_of(aa))
    chr = t(chr)
    colnames(chr) = varIDs
    chr = cbind(data.frame(row.names(chr)), chr)
    colnames(chr)[1] = 'Id'

    cat('Number of individuals: ', nrow(chr), '\n')
    
    output_dir <- paste('/Volumes/im-lab/nas40t2/jiamao/Projects/cancer_PredictDB/TCGA/TCGA/data/hg19/genotype/genotype/splited_dosages/', project_name, sep = '')
    if (!dir.exists(output_dir)){
      dir.create(output_dir)
    } else {
      cat(output_dir, " Dir already exists!\n")
    }
    
    # get tcga cancer sample size information 
    if (i==1 && project_name == 'TCGA-BRCA') {
      summary_table = rbind(summary_table, t(data.frame(c(project_name, nrow(chr)))))
    }
    
    cat('writing data into disk at the path - ', output_dir, '\n')
    write.csv(chr, file = gzfile(paste(output_dir, '/', project_name, '_filtered_biallelic', '_chr', i, '.csv.gz', sep = '')),  quote = FALSE, row.names = F)
    cat('done!\n\n')
  }
} 

cat('writing cancer sample size into disk at the path - /Volumes/im-lab/nas40t2/jiamao/Projects/cancer_PredictDB/TCGA/TCGA/data/hg19/genotype/genotype/splited_dosages/. \n\n')
colnames(summary_table) = c('project_name', 'n')
write.csv(summary_table, file = '/Volumes/im-lab/nas40t2/jiamao/Projects/cancer_PredictDB/TCGA/TCGA/data/hg19/genotype/genotype/splited_dosages/cancer_sample_size.csv', row.names = F, quote = F)
cat('done!\n\n')
