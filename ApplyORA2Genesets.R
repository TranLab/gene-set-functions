#ApplyORA2Genesets

#Description
#This function takes a list of genes and applies the fora function from the fgsea package using multiple gene sets relevant to immunology and blood transcriptomics.

#Usage
#ApplyORA2Genesets <- function(genelist, geneset = c("all", "bloodmodules", "MSigDB"), output_directory = getwd(), filename_prefix = "myFile", fixed_seed = TRUE, ...)

#Arguments
#genelist             A vector of HGNC gene symbols representing genes of interest
#universe             A universe from which 'genelist' was selected
#geneset              Can be "bloodmodules" which includes low- and high-annotation blood transcription modules (Li et al. Nat Immunol 2014 [PMID 24336226]; Kazmin et al PNAS 2017 [PMID 28193898]); Monaco modules (Monaco et al. Cell Reports 2019 [PMID 30726743]); and BloodGen3Modules (Rinchai et al. Binformatics 2021 [PMID 33624743] or "MSigDB" which includes most collections from the Molecular Signature Database v7.4 (https://www.gsea-msigdb.org/gsea/msigdb/) or "all" which includes both.  
#output_directory     User-defined directory path
#filename_prefix      User-defined pre-fix for filename where results will be saved, defaults to "myFile"
#fixed_seed           If TRUE, a fixed seed will be set for reproducibility of exact NES and significance values.
# ...                 All other arguments passed onto fgseaMultilevel

#Value
#Output is a single table of all gene sets. Columns are module_type (e.g. "highBTMs" or "MonacoModules") followed all the columns in the table generated by fgseaMultilevel (see "Value" in fgseaMultilevel help) 


ApplyORA2Genesets <- function(genelist, geneset = c("all", "bloodmodules", "MSigDB"), universe, output_directory = getwd(), filename_prefix = "myFile", fixed_seed = TRUE, ...) {
  
  library(tidyverse)
  library(fgsea)
  library(data.table)
  library(curl)
  
  if(fixed_seed == TRUE){
    set.seed(12345)
  }
  if(class(genelist) == "data.frame" | class(genelist) == "tbl" | class(genelist) == "data.table"){
    genelist <- deframe(genelist)
  } 
  genelist <- sort(genelist, decreasing = TRUE)
  #for rank lists  that still have original "MARCH" or "SEPT" gene names, replace with "MARCHF" or "SEPTIN"
  if(length(names(genelist)[grepl("MARCHF", names(genelist))]) == 0){
    names(genelist) <- gsub("MARCH", "MARCHF", names(genelist))
    }
  if(length(names(genelist)[grepl("SEPTIN", names(genelist))]) == 0){
    names(genelist) <- gsub("SEPT", "SEPTIN", names(genelist))
    }
  #Load the pathways into a named list
  if(geneset == "all"){
    allGeneSets <- c("MonacoModules", "highBTMs", "lowBTMs", "BloodGen3Module",
                     "MSigDB_Hallmark_v7.4", "MSigDB_C2_biocarta_v7.4", "MSigDB_C2_kegg_v7.4", "MSigDB_C8_all_v7.4","MSigDB_C7_all_v7.4","MSigDB_C5_GO_bp_v7.4", "MSigDB_C5_GO_mf_v7.4")
    }
  if(geneset == "bloodmodules"){
    allGeneSets <- c("MonacoModules", "highBTMs", "lowBTMs", "BloodGen3Module")
    }
  if(geneset == "MSigDB"){
    allGeneSets <- c("MSigDB_Hallmark_v7.4", "MSigDB_C2_biocarta_v7.4", "MSigDB_C2_kegg_v7.4", "MSigDB_C8_all_v7.4","MSigDB_C7_all_v7.4","MSigDB_C5_GO_bp_v7.4", "MSigDB_C5_GO_mf_v7.4")
                     }
  if(!geneset %in% c("all", "bloodmodules", "MSigDB")){
    print("geneset argument can only be 'all', 'bloodmodules', or 'MSigDB'")
    stop()
    }
  myModuleList <- oraRes <- c()
  for(i in allGeneSets){
    temp <- tempfile(fileext = ".rds")
    url <- paste0("https://github.com/TranLab/ModuleLists/blob/main/", i, ".rds?raw=true")
    myModuleList[[i]] <- readRDS(url(url, method="libcurl"))
    #for module lists that still have original "MARCH" or "SEPT" gene names, replace with "MARCHF" or "SEPTIN"
    if(length(unlist(myModuleList[[i]])[grepl("MARCHF", unlist(myModuleList[[i]]))]) == 0){
      for(j in 1:length(myModuleList[[i]])){
        myModuleList[[i]][[j]] <- gsub("MARCH", "MARCHF", myModuleList[[i]][[j]])
      }
      }
    if(length(unlist(myModuleList[[i]])[grepl("SEPTIN", unlist(myModuleList[[i]]))]) == 0){
      for(j in 1:length(myModuleList[[i]])){
        myModuleList[[i]][[j]] <- gsub("SEPT", "SEPTIN", myModuleList[[i]][[j]])
      }
      }
    print(paste0("Running ORA for ", i, " signatures"))
    #for info on the SampleSize argument: https://www.biostars.org/p/479821/
    #run fgsea with 1000 permutations
    
    oraRes[[i]] <- fora(pathways=myModuleList[[i]], genes=genelist, universe=universe, ...) %>%
      as_tibble() %>%
      # dplyr::select(pathway, overlap, size, leadingEdge, pval, padj) %>% 
      # mutate(leadingEdge = gsub("^c\\(|\\)$", "", leadingEdge)) %>%
      # mutate(leadingEdge = gsub('"', "", leadingEdge)) %>%
      arrange(padj) %>%
      as.data.frame()
    }
  ORAtab <- oraRes %>%
    dplyr::bind_rows(.id = "module_type")
  writexl::write_xlsx(oraRes, paste0(output_directory, filename_prefix, " ORA results ", geneset, ".xlsx"))
  return(oraRes)
  close(url)
  }
