# gene-set-functions

Convenience functions for gene set and module analysis for transcriptomic data

# Gene2ModuleExpressionScores function

## Description

This function collapses a gene expression matrix to module-level expression scores

## Usage

Gene2ModuleExpressionScores(gene_expression_dat, module_list = c("lowBTMs", "highBTMs", "BloodGen3Module","MonacoModules"), summary_stat = c(mean, median)) 

| Arguments   | Description |
| :--- | :- |
| gene_expression_dat       | gene expression matrix of normalized logCPM values (not counts) or ExpressionSet object containing such a matrix with rownames as HUGO symbols or Ensembl ID. |
| module_list | Name of module set to use. Can be "lowBTMs", "highBTMs", "BloodGen3Module" or "MonacoModules". |
| summary_stat       | Mean or median. |  
       
## Value

Output is a dataframe with rownames as module names and column values representing the mean or median or member genes within a module


# NamedGeneRankList2GseaTable

## Description

This function takes a named ranked list of genes and applies the fgseaMultilevel function from the fgsea package using multiple gene sets relevant to immunology and blood transcriptomics.

## Usage

NamedGeneRankList2GseaTable(rankedgenes, geneset = c("all", "bloodmodules", "MSigDB"), output_directory = getwd(), filename_prefix = "myFile", fixed_seed = TRUE, ...)

| Arguments   | Description |
| :--- | :- |
| rankedgenes       | A named vector of a ranking metric such as log2 fold-change with names being HUGO gene symbols. Can also be two column dataframe or tibble       where the first column has gene names and the second column is the ranking metric. Vector automatically gets sorted in descending order of ranking metric. |
| geneset | Can be "bloodmodules" which includes low- and high-annotation blood transcription modules (Li et al. Nat Immunol 2014 [PMID 24336226]; Kazmin et al PNAS 2017 [PMID 28193898]); Monaco modules (Monaco et al. Cell Reports 2019 [PMID 30726743]); and BloodGen3Modules (Rinchai et al. Binformatics 2021 [PMID 33624743] or "MSigDB" which includes most collections from the Molecular Signature Database v7.4 (https://www.gsea-msigdb.org/gsea/msigdb/) or "all" which includes both.|
| output_directory       | User-defined directory path.|  
| filename_prefix       | User-defined pre-fix for filename where results will be saved, defaults to "myFile". |
| fixed_seed       | If TRUE, a fixed seed will be set for reproducibility of exact NES and significance values. | 
| ...       | All other arguments passed onto fgseaMultilevel. |

## Value

Output is a single table of all gene sets. Columns are module_type (e.g. "highBTMs" or "MonacoModules") followed all the columns in the table generated by fgseaMultilevel (see "Value" in fgseaMultilevel help) 