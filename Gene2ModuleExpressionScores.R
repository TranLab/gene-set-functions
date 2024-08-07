#Gene2ModuleExpressionScores


#Description
#This function collapses a gene expression matrix to module-level expression scores

#Usage
#Gene2ModuleExpressionScores(gene_expression_dat, module_list = c("lowBTMs", "highBTMs", "BloodGen3Module","MonacoModules"), summary_stat = c(mean, median)) 

#Arguments
#gene_expression_dat    gene expression matrix of normalized logCPM values (not counts) or ExpressionSet object containing such a matrix with rownames as HUGO symbols or Ensembl IDs
#module_list            name of module set to use. Can be "lowBTMs", "highBTMs", "BloodGen3Module" or "MonacoModules".
#summary_stat           mean or median

#Value
#Output is a dataframe with rownames as module names and column values representing the mean or median or member genes within a module

library(tidyr)
library(Biobase)
library(curl)
library(data.table)

Gene2ModuleExpressionScores <- function(gene_expression_dat, module_list = c("lowBTMs", "highBTMs", "BloodGen3Module","MonacoModules"), summary_stat = c(mean, median)) {
  temp <- tempfile(fileext = ".rds")
  url <- paste0("https://github.com/TranLab/ModuleLists/blob/main/",module_list, ".rds?raw=true")
  myModuleList <- readRDS(url(url, method="libcurl"))
  temp_df <- c()
  for(i in names(myModuleList)){
    temp_df[[i]] <- data.frame("module" = i, "GeneSymbol" = myModuleList[[i]])
  }
  temp_df <- bind_rows(temp_df)
  
  
  if(class(gene_expression_dat)[1] == "ExpressionSet"){
    exprs_dat_temp <- Biobase::exprs(gene_expression_dat) %>%
      as.data.frame() %>%
      rownames_to_column(var = "geneid")
  }
  if(class(gene_expression_dat)[1] == "matrix"){
    exprs_dat_temp <- gene_expression_dat %>%
      as.data.frame() %>%
      rownames_to_column(var = "geneid")
  }
  

  if(grepl("ENSG", exprs_dat_temp$geneid[1])){
    library('biomaRt')
    ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    exprs_dat_temp <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), mart = ensembl) %>%
      dplyr::rename(geneid = "ensembl_gene_id") %>%
      dplyr::rename(GeneSymbol = "hgnc_symbol") %>%
      distinct(geneid, .keep_all = TRUE) %>%
      right_join(., exprs_dat_temp,
                 by = "geneid") %>%
      dplyr::select(-c(geneid))
  }else
    {
    exprs_dat_temp <- exprs_dat_temp %>%
      dplyr::rename(GeneSymbol = "geneid")
    }
  exprs_dat_temp2 <- exprs_dat_temp %>%
    left_join(., temp_df,
              by="GeneSymbol") %>%
    dplyr::select(module, everything()) %>%
    dplyr::select(-c(GeneSymbol)) %>%
    filter(!is.na(module)) %>%
    as.data.table()
    exprs_dat_temp2  <- exprs_dat_temp2[, lapply(.SD, summary_stat, na.rm=TRUE), by=module ]
    exprs_dat_temp2 <- exprs_dat_temp2 %>%
      as.data.frame() %>%
      column_to_rownames(var = "module")
  return(exprs_dat_temp2)
  rm(temp_df, myModuleList,exprs_dat_temp, exprs_dat_temp2)
}

  
