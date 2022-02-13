################################################################################
# Required libraries
################################################################################

library(clusterProfiler)
library(ggnewscale)
library(enrichplot)
library(ggplot2)
library(DOSE)
library(splitstackshape)
library(dplyr)
library(org.Hs.eg.db) # Genome wide annotation for Human
library(genefilter)
library(reshape2)
library(msigdbr)
library(hpar) # The Human Protein Atlas
library(ggridges)
library(ggupset)
library(plotly)
library(forcats)
library(stringr)
library(UniprotR) # Connect to Uniprot to retrieve information about proteins
library(ggpubr)
library(pdp)
library(gridExtra)


################################################################################
# Functions
################################################################################


################################################################################
# Function 1: upload.data.function 
#             upload the proteomics matrix under study
################################################################################

upload.data.function <- function(data, 
                                 header = TRUE, 
                                 sep = "\t",
                                 dec = ".", 
                                 remove.na = TRUE) {
  
  prot_list <- read.csv(data, header = header, sep = sep, dec = dec)
  if (remove.na == TRUE) {
    prot_list <- prot_list %>% na.omit()
  }
  
  # save uploaded data 
  dir.create(file.path(getwd(), "uploaded.data"), showWarnings = FALSE)
  write.table(prot_list, "uploaded.data/uploaded.data.csv", row.names=FALSE, quote=FALSE,
              sep = "\t")
  
  return(prot_list)
}


################################################################################
# Function 2: convert.ids.function
#             Biological Id TRanslator using `bitr` 
################################################################################

# First column needs to be the id.
convert.ids.function <- function(data, 
                                 fromType, 
                                 toType = "ENTREZID",
                                 OrgDb = org.Hs.eg.db,
                                 rm.duplicates = TRUE) {
  
  entrez_uniprot_dict <- bitr(data[, 1], 
                              fromType = fromType, 
                              toType = toType,
                              OrgDb = OrgDb, 
                              drop = TRUE) # drop NAs
  
  # Merge the data and the created dictionary by "fromType" common column:
  entrez_uniprot_expression <- merge(entrez_uniprot_dict, 
                                     data, 
                                     by.x = colnames(entrez_uniprot_dict[1]), 
                                     by.y = colnames(data[1]),
                                     all=TRUE) 
  
  # Retain all columns except "fromType" so we just have the entrez id relative
  # abundances and remove those ids not translated:
  entrez_expression <- entrez_uniprot_expression[, !names(entrez_uniprot_expression) 
                                                 %in% fromType] 
  entrez_expression <- entrez_expression[!is.na(entrez_expression[[toType]]),]
  
  # Remove duplicates:
  if (rm.duplicates == TRUE) {
    entrez_expression <- entrez_expression[
      which(duplicated(entrez_expression[[toType]]) == FALSE),] 
  }
  
  return(entrez_expression)
}


################################################################################
# Function 3: peptide.function
#             returns a list containing all expression dataframes made up of
#             detected proteins (ENTREZ) depending on the number of peptides
################################################################################

peptide.function <- function(entrez_expression, 
                             peptides = seq(1:5), 
                             Np.col) {
  
  entrez.expr.peptides <- list()
  prot_count <- vector()
  num.peptides <- vector()
  plot.num.peptides <- vector()
  
  for (i in 1:length(peptides)) {
    entrez.expr <- paste0("entrez.expression_", i)
    
    if (grepl(">", peptides[i]) == TRUE) {
      
      assign(entrez.expr, 
             entrez_expression[entrez_expression[[Np.col]] > i-1, ])
      
    } else if (grepl("-", peptides[i]) == TRUE) {
      
      pep <- eval(parse(text = (str_replace(peptides[i], "-", ":")) ))
      
      assign(entrez.expr, 
             entrez_expression[entrez_expression[[Np.col]] %in% pep, ])
      
    } else{
      
      assign(entrez.expr, 
             entrez_expression[entrez_expression[[Np.col]] == i, ])
      
    }
    
    # protein expression
    entrez.expr.peptides[[i]] <- get(entrez.expr)
    # count vectors
    prot_count[i] <- nrow(get(entrez.expr))
    
    if (grepl(">", peptides[i]) == TRUE) {
      
      num.peptides[i] <- paste0(eval(parse(text = substr(peptides[i], 2, str_length(peptides[i]))))+1, 
                                ".or.more.peptides")
      
    } else if (grepl("-", peptides[i]) == TRUE) {
      
      num.peptides[i] <- paste0(str_replace(peptides[i], "-", "_"), 
                                ".peptides")
      
    } else {
      
      if (peptides[i] == 1) {
        num.peptides[i] <- paste0(peptides[i], ".peptide")
      } else {
        num.peptides[i] <- paste0(peptides[i], ".peptides")
      }
      
    }
    
    plot.num.peptides[i] <- peptides[i]
    
    # create a folder for each protein abundance matrix
    dir.create(file.path(getwd(), num.peptides[i]), showWarnings = FALSE)
    prot.abundance.matrix <- convert.ids.function(entrez.expr.peptides[[i]],
                                                  fromType = "ENTREZID",
                                                  toType = "UNIPROT")
    write.table(prot.abundance.matrix, paste0(num.peptides[i], "/protein.abundance.matrix.", num.peptides[i], ".csv"), 
                row.names=FALSE, quote=FALSE, sep = "\t")
  }
  
  # dataframe containing values to plot
  count_dataframe <- data.frame(row.names = seq(1:length(peptides)))
  count_dataframe$peptides <- plot.num.peptides
  count_dataframe$proteins <- prot_count
  # save in summary folder
  dir.create(file.path(getwd(), "summary"), showWarnings = FALSE)
  # create detected proteins folder
  dir.create(file.path(paste0(getwd(), "/summary"), "detected.proteins"), showWarnings = FALSE)
  write.csv(count_dataframe, "summary/detected.proteins/detected.proteins.csv", 
            row.names=FALSE, quote=FALSE)
  
  # remove unnecessary variables
  rm(list = entrez.expr, prot_count, num.peptides)
  
  entrez_expression_list <- list(entrez.expr.peptides, 
                                 count_dataframe)
  
  return(entrez_expression_list)
}


################################################################################
# Function 4: uniprot.peptide.function
#             returns a list containing all expression dataframes made up of
#             detected proteins (UniProt) depending on the number of peptides
################################################################################

uniprot.peptide.function <- function(data, 
                             peptides = seq(1:5), 
                             Np.col) {
  
  entrez.expr.peptides <- list()
  prot_count <- vector()
  num.peptides <- vector()
  plot.num.peptides <- vector()
  
  for (i in 1:length(peptides)) {
    entrez.expr <- paste0("entrez.expression_", i)
    
    if (grepl(">", peptides[i]) == TRUE) {
      
      assign(entrez.expr, 
             data[data[[Np.col]] > i-1, ])
      
    } else if (grepl("-", peptides[i]) == TRUE) {
      
      pep <- eval(parse(text = (str_replace(peptides[i], "-", ":")) ))
      
      assign(entrez.expr, 
             data[data[[Np.col]] %in% pep, ])
      
    } else{
      
      assign(entrez.expr, 
             data[data[[Np.col]] == i, ])
      
    }
    
    # protein expression
    entrez.expr.peptides[[i]] <- get(entrez.expr)
    # count vectors
    prot_count[i] <- nrow(get(entrez.expr))
    
    if (grepl(">", peptides[i]) == TRUE) {
      
      num.peptides[i] <- paste0(eval(parse(text = substr(peptides[i], 2, str_length(peptides[i]))))+1, 
                                ".or.more.peptides")
      
    } else if (grepl("-", peptides[i]) == TRUE) {
      
      num.peptides[i] <- paste0(str_replace(peptides[i], "-", "_"), 
                                ".peptides")
      
    } else {
      
      if (peptides[i] == 1) {
        num.peptides[i] <- paste0(peptides[i], ".peptide")
      } else {
        num.peptides[i] <- paste0(peptides[i], ".peptides")
      }
      
    }
    
    plot.num.peptides[i] <- peptides[i]
    
    
  }
  
  # dataframe containing values to plot
  count_dataframe <- data.frame(row.names = seq(1:length(peptides)))
  count_dataframe$peptides <- plot.num.peptides
  count_dataframe$proteins <- prot_count
  
  
  entrez_expression_list <- list(entrez.expr.peptides, 
                                 count_dataframe)
  
  return(entrez_expression_list)
}


################################################################################
# Function 5: t.test.function
#             performs a differential abundance analysis using a 
#             Welch two sample t-test function followed by a FDR
#             multiple testing correction
################################################################################

t.test.function <- function(entrez_expression, 
                            control.cols, 
                            case.cols,
                            log.scale = TRUE) {
  
  groups <- factor(rep(c("control", "case"), c(length(control.cols),
                                               length(case.cols))))
  

  # perform t-test
  ################
  
  ttest_result <- rowttests(as.matrix(entrez_expression[,control.cols[1]:case.cols[length(case.cols)]]), 
                            factor(groups), na.rm = TRUE) 
  ttest_result <- ttest_result %>% mutate(
    adjusted.p.value = p.adjust(ttest_result$p.value, "fdr"))
  
  # dm column output from t-test function indicates the difference of the group
  # (controls vs cases) means, this is, the fold change of the two groups under
  # study, especifically it is the log2foldchange
  # Rename the column to be clearer
  
  colnames(ttest_result)[which(colnames(ttest_result)=="dm")] <- "fold_change"
  
  # Join the previously calculated statistics 
  entrez_expression <- cbind(entrez_expression, ttest_result)
  
  if (log.scale == TRUE) {
    # To deal with 0 and negative values:
    abs_log <- function(x) {
      x[x==0] <- 1
      si <- sign(x)
      si * log2(si*x)
    }
    
    # use entrez_expression[] to preserve data frame structure
    entrez_expression[, "fold_change"] <- 
      sapply(entrez_expression[, "fold_change"],
             abs_log)
    
  }
  
  # sort in descending order by fold change:
  entrez_expression <- entrez_expression %>% arrange(desc(fold_change))
  
  return(entrez_expression)
}


################################################################################
# Function 6: gsea.function
#             performs a functional profiling using a Gene Set Enrichment 
#             Analysis (GSEA)
################################################################################

gsea.function <- function(entrez_expression, 
                          metric = "fold_change", 
                          ont = "all", 
                          OrgDb = org.Hs.eg.db, 
                          gsea = "all", 
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "fdr",
                          minGSSize = 10,
                          maxGSSize = 500) {
  
  # vector containing desc arranged log2foldchange with respective entrez ids
  entrez_list <- entrez_expression[, metric]
  names(entrez_list) <- as.character(entrez_expression[, "ENTREZID"]) 
  
  
  if (gsea == "all") {
    gsea <- c("GO", "KEGG", "DISGENET", "PHENOTYPE", "CELL", "TISSUE")
  }
  
  
  #---------------------------
  # save gsea objects 
  gsea.results <- list()
  all.gsea.results <- list()
  # save category counts 
  count.df <- data.frame()
  count.gsea.results <- list()
  #---------------------------
  
  
  if ("GO" %in% gsea) {
    # for multiple ontologies
    if (ont == "all") {
      ont <- c("BP", "MF", "CC")
    }
    for (j in 1:length(ont)) {
      if (isEmpty(names(entrez_list)) == FALSE) {
        # Gene Ontology (GO)
        gsea_go <- gseGO(entrez_list,  
                         keyType = "ENTREZID",
                         OrgDb = OrgDb, 
                         ont = ont[j], 
                         pvalueCutoff = pvalueCutoff,
                         pAdjustMethod = pAdjustMethod,
                         minGSSize = minGSSize,
                         maxGSSize = maxGSSize)
        
        #In case no mapping has been performed (else):
        if (!is.null(gsea_go)) {
          gsea.results$GO[[ont[j]]] <- gsea_go
        } else {
          gsea.results$GO[[ont[j]]] <- setNames(data.frame(matrix(ncol = 11, nrow = 0)), 
                                                c("ID", "Description", "setSize", "enrichmentScore", 
                                                  "NES", "pvalue", "p.adjust", "qvalues", "rank", 
                                                  "leading_edge", "core_enrichment"))
        }
        
        
        go.name <- paste0("GO(", ont[j], ")")
        count.gsea.results[[go.name]] <- nrow(as.data.frame(gsea_go))
      }
      else {
        gsea.results$GO[[ont[j]]] <- setNames(data.frame(matrix(ncol = 11, nrow = 0)), 
                                              c("ID", "Description", "setSize", "enrichmentScore", 
                                                "NES", "pvalue", "p.adjust", "qvalues", "rank", 
                                                "leading_edge", "core_enrichment"))
        go.name <- paste0("GO(", ont[j], ")")
        count.gsea.results[[go.name]] <- 0
      }
    }
  }
  
  
  if ("DISGENET" %in% gsea) {
    if (isEmpty(names(entrez_list)) == FALSE) {
      # Gene Ontology (GO)
      gsea_dgn <- gseDGN(entrez_list, 
                         pvalueCutoff = pvalueCutoff,
                         pAdjustMethod = pAdjustMethod,
                         minGSSize = minGSSize,
                         maxGSSize = maxGSSize)
      
      #In case no mapping has been performed (else):
      if (!is.null(gsea_dgn)) {
        gsea.results$DISGENET <- gsea_dgn
      } else {
        gsea.results$DISGENET <- setNames(data.frame(matrix(ncol = 11, nrow = 0)), 
                                          c("ID", "Description", "setSize", "enrichmentScore", 
                                            "NES", "pvalue", "p.adjust", "qvalues", "rank", 
                                            "leading_edge", "core_enrichment"))
      }
      
      count.gsea.results$DISGENET <- nrow(as.data.frame(gsea_dgn))
    }
    else {
      gsea.results$DISGENET <- setNames(data.frame(matrix(ncol = 11, nrow = 0)), 
                                        c("ID", "Description", "setSize", "enrichmentScore", 
                                          "NES", "pvalue", "p.adjust", "qvalues", "rank", 
                                          "leading_edge", "core_enrichment"))
      count.gsea.results$DISGENET <- 0
    }
  }
  
  
  if ("KEGG" %in% gsea) {
    if (isEmpty(names(entrez_list)) == FALSE) {
      # Gene Ontology (GO)
      gsea_kegg <- gseKEGG(entrez_list,
                           organism = "hsa",
                           pAdjustMethod = pAdjustMethod,
                           minGSSize = minGSSize,
                           maxGSSize = maxGSSize)
      
      if (!is.null(gsea_kegg)) {
        gsea.results$KEGG <- gsea_kegg
      } else {
        gsea.results$KEGG <- setNames(data.frame(matrix(ncol = 11, nrow = 0)), 
                                      c("ID", "Description", "setSize", "enrichmentScore", 
                                        "NES", "pvalue", "p.adjust", "qvalues", "rank", 
                                        "leading_edge", "core_enrichment"))
      }
      
      gsea.results$KEGG <- gsea_kegg
      count.gsea.results$KEGG <- nrow(as.data.frame(gsea_kegg))
    }
    else {
      gsea.results$KEGG <- setNames(data.frame(matrix(ncol = 11, nrow = 0)), 
                                    c("ID", "Description", "setSize", "enrichmentScore", 
                                      "NES", "pvalue", "p.adjust", "qvalues", "rank", 
                                      "leading_edge", "core_enrichment"))
      count.gsea.results$KEGG <- 0
    }
  }
  
  
  ## Phenotype Human Ontology (PHO)
  if ("PHENOTYPE" %in% gsea) {
    if (isEmpty(names(entrez_list)) == FALSE) {
      # Download Human Phenotype Ontology
      human.phenotype.ontology <-  msigdbr(species = "Homo sapiens", 
                                           category = "C5", 
                                           subcategory = "HPO")
      
      # select entrezid and term columns to use in cluster profiler (TERM TO GENE)
      hpo_term.to.gene <- human.phenotype.ontology[, c("gs_id", "entrez_gene")]
      hpo_term.to.name <- human.phenotype.ontology[, c("gs_id", "gs_name")]
      
      gsea_PHENOTYPE <- GSEA(gene = entrez_list,
                             pAdjustMethod = pAdjustMethod,
                             minGSSize = minGSSize,
                             maxGSSize = maxGSSize,
                             pvalueCutoff = pvalueCutoff,
                             TERM2GENE = hpo_term.to.gene,
                             TERM2NAME = hpo_term.to.name)
      
      #In case no mapping has been performed (else):
      if (!is.null(gsea_PHENOTYPE)) {
        gsea.results$PHENOTYPE <- gsea_PHENOTYPE
      } else {
        gsea.results$PHENOTYPE <- setNames(data.frame(matrix(ncol = 8, nrow = 0)), 
                                           c("ID", "Description", "setSize", "enrichmentScore", 
                                             "NES", "pvalue", "p.adjust", "qvalues"))
      }
      
      
      count.gsea.results$PHENOTYPE <- nrow(as.data.frame(gsea_PHENOTYPE))
    }
    else {
      gsea.results$PHENOTYPE <- setNames(data.frame(matrix(ncol = 8, nrow = 0)), 
                                         c("ID", "Description", "setSize", "enrichmentScore", 
                                           "NES", "pvalue", "p.adjust", "qvalues"))
      count.gsea.results$PHENOTYPE <- 0 
    }
  }
  
  # CELL TYPE (GSEA)
  if ("CELL" %in% gsea) {
    if (isEmpty(names(entrez_list)) == FALSE) {
      
      
      cell.signature <-  msigdbr(species = "Homo sapiens", 
                                 category = "C8")
      
      cell.to.gene <- cell.signature[, c("gs_id", "entrez_gene")]
      cell.to.name <- cell.signature[, c("gs_id", "gs_name")]
      
      gsea_CELL <- GSEA(gene = entrez_list,
                        pAdjustMethod = pAdjustMethod,
                        minGSSize = minGSSize,
                        maxGSSize = maxGSSize,
                        pvalueCutoff = pvalueCutoff,
                        TERM2GENE = cell.to.gene,
                        TERM2NAME = cell.to.name)
      
      #In case no mapping has been performed (else):
      if (!is.null(gsea_CELL)) {
        gsea.results$CELL <- gsea_CELL
      } else {
        gsea.results$CELL <- setNames(data.frame(matrix(ncol = 8, nrow = 0)), 
                                      c("ID", "Description", "setSize", "enrichmentScore", 
                                        "NES", "pvalue", "p.adjust", "qvalues"))
      }
      
      
      count.gsea.results$CELL <- nrow(as.data.frame(gsea_CELL))
    }
    else {
      gsea.results$CELL <- setNames(data.frame(matrix(ncol = 8, nrow = 0)), 
                                    c("ID", "Description", "setSize", "enrichmentScore", 
                                      "NES", "pvalue", "p.adjust", "qvalues"))
      count.gsea.results$CELL <- 0 
    }
  }
  
  # Tissue Human Protein Atlas (hpar package)
  if ("TISSUE" %in% gsea) {
    if (isEmpty(names(entrez_list)) == FALSE) {
      
      data(hpaNormalTissue)
      # convert ids to entrezid
      hpaNormalTissue <- convert.ids.function(hpaNormalTissue, 
                                              fromType = "ENSEMBL",
                                              toType = "ENTREZID")
      # assign an id
      hpaNormalTissue <- transform(hpaNormalTissue,id=as.numeric(factor(Tissue)))
      
      tissue.to.gene <- hpaNormalTissue[, c("id", "ENTREZID")]
      tissue.to.name <- hpaNormalTissue[, c("id", "Tissue")]
      
      gsea_TISSUE <- GSEA(gene = entrez_list,
                          pAdjustMethod = pAdjustMethod,
                          minGSSize = minGSSize,
                          maxGSSize = maxGSSize,
                          pvalueCutoff = pvalueCutoff,
                          TERM2GENE = tissue.to.gene,
                          TERM2NAME = tissue.to.name)
      
      #In case no mapping has been performed (else):
      if (!is.null(gsea_TISSUE)) {
        gsea.results$TISSUE <- gsea_TISSUE
      } else {
        gsea.results$TISSUE <- setNames(data.frame(matrix(ncol = 8, nrow = 0)), 
                                        c("ID", "Description", "setSize", "enrichmentScore", 
                                          "NES", "pvalue", "p.adjust", "qvalues"))
      }
      
      count.gsea.results$TISSUE <- nrow(as.data.frame(gsea_TISSUE))
    }
    else {
      gsea.results$TISSUE <- setNames(data.frame(matrix(ncol = 8, nrow = 0)), 
                                      c("ID", "Description", "setSize", "enrichmentScore", 
                                        "NES", "pvalue", "p.adjust", "qvalues"))
      count.gsea.results$TISSUE <- 0 
    }
  }
  
  
  ## Transform dataframe into large formart using melt function
  count.df <- melt(as.data.frame(count.gsea.results), 
                   value.name = "Count", 
                   variable.name = "Category")
  
  
  # RESULTS
  gsea_results_list <- list(gsea.results,
                            count.df)
  
  
  return(gsea_results_list)
}


################################################################################
# Function 7: ora.function
#             performs a functional profiling using an Over Representation 
#             Analysis (ORA)  
################################################################################

ora.function <- function(entrez_expression, 
                         metric = "fold_change",
                         ont = "all",
                         FDR.name.col = "FDR",
                         regulation = "all",
                         OrgDb = org.Hs.eg.db,
                         ora = "all",
                         pvalueCutoff = 0.05, 
                         pAdjustMethod = "fdr",
                         minGSSize = 10,
                         maxGSSize = 500,
                         fdr = 0.05,
                         fc = c(1.5, 2)) {
  
  ## consider as DEGs those detected proteins with FDR<0.05:
  DEGs <- entrez_expression[entrez_expression[[FDR.name.col]]<fdr,]
  
  ## Define numeric vectors containing the ENTREZ id with their respective log2foldchange.
  DEGs_list <- DEGs[[metric]]
  names(DEGs_list) <- as.character(DEGs$ENTREZID) 
  
  # save ora objects and detected DEGs
  ora.results <- list()
  all.ora.results.reg <- list()
  all.ora.results <- list()
  DEGs_fc_reg <- list()
  DEGs_fc <- list()
  
  # save category counts 
  count.df <- data.frame()
  count.ora.results <- list()
  
  if (ora == "all") {
    ora <- c("GO", "KEGG", "DISGENET", "PHENOTYPE", "CELL", "TISSUE")
  }
  
  if (ont == "all") {
    ont <- c("BP", "MF", "CC")
  }
  
  if (regulation == "all") {
    regulation = c("UP", "DOWN", "UP.and.DOWN")
  }
  
  # fold change (i) loop
  for (i in 1:length(fc)) {
    
    # regulation (j) loop
    
    for (j in 1:length(regulation)){
      
      # define differentially expressed proteins:
      if (regulation[j] == "UP") {
        deg <- DEGs_list[DEGs_list > fc[i]]
      } else if (regulation[j] == "DOWN") {
        deg <- DEGs_list[DEGs_list < -fc[i]]
      } else if (regulation[j] == "UP.and.DOWN") {
        deg <- DEGs_list[abs(DEGs_list) > fc[i]]
      }
      
      # create fold change named list of differentially expressed proteins:
      fold <- noquote(paste0("fold_change=",fc[i]))
      reg <- paste0(regulation[j], ".regulated")
      
      
      # save the number of detected DEGs depending on the foldchange and regulation
      # count.ora.results is a list
      count.ora.results$regulation[j] <- regulation[j]
      count.ora.results$fold_change[j] <- fc[i]
      count.ora.results$DEGs[j] <- length(deg)
      
      
      ## define protein universe numeric vector
      universe <- entrez_expression[[metric]]
      names(universe) <- as.character(entrez_expression$ENTREZID)
      
      # perform ORA for each ora.test
      #-------------------------------------------------------------------------
      
      ## Gene Ontology (GO) ----------------------------------------------------
      
      if ("GO" %in% ora) {
        
        # specific GO annotations loop (k)
        
        for (k in 1:length(ont)) {
          
          ## ORA analysis only if any DEG has been detected
          if (isEmpty(names(deg)) == FALSE) {
            
            ora_go <- enrichGO(gene = names(deg),
                               universe = names(universe), #background genes
                               OrgDb = org.Hs.eg.db, 
                               keyType = "ENTREZID",
                               readable = T,
                               ont = ont[k],
                               minGSSize = minGSSize,
                               maxGSSize = maxGSSize,
                               pAdjustMethod = pAdjustMethod,
                               pvalueCutoff = pvalueCutoff)
            
            
            #In case no mapping has been performed (else):
            if (!is.null(ora_go)) {
              ora.results$GO[[ont[k]]] <- ora_go
            } else {
              ora.results$GO[[ont[k]]] <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), 
                                                   c("ID", "Description", "GeneRatio", "BgRatio", 
                                                     "pvalue", "p.adjust", "qvalue", "geneID", "Count"))
            }
            
            # save detected ora categories
            go.name <- paste0("GO(", ont[k], ")")
            count.ora.results[[go.name]][j] <- nrow(as.data.frame(ora_go))
            
          }
          
          else {
            
            # save ora enrichment 
            ora.results$GO[[ont[k]]] <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), 
                                                 c("ID", "Description", "GeneRatio", "BgRatio", 
                                                   "pvalue", "p.adjust", "qvalue", "geneID", "Count"))
            # save detected ora categories
            go.name <- paste0("GO(", ont[k], ")")
            count.ora.results[[go.name]][j] <- 0 
          }
          
        } # end specific GO annotations loop (k)
        
      } # END GO condition 
      
      ## DISGENET --------------------------------------------------------------
      
      if ("DISGENET" %in% ora) {
        
        if (isEmpty(names(deg)) == FALSE) {
          
          ora_DGN <- enrichDGN(gene = names(deg),
                               universe = names(universe), 
                               pAdjustMethod = pAdjustMethod,
                               minGSSize = minGSSize,
                               maxGSSize = maxGSSize,
                               pvalueCutoff = pvalueCutoff)
          
          #In case no mapping has been performed (else):
          if (!is.null(ora_DGN)) {
            ora.results$DISGENET <- ora_DGN
          } else {
            ora.results$DISGENET <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), 
                                             c("ID", "Description", "GeneRatio", "BgRatio", 
                                               "pvalue", "p.adjust", "qvalue", "geneID", "Count"))
          }
          
          # save detected ora categories
          count.ora.results$DISGENET[j] <- nrow(as.data.frame(ora_DGN))
          
        }
        
        else {
          
          # save ora enrichment
          ora.results$DISGENET <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), 
                                           c("ID", "Description", "GeneRatio", "BgRatio", 
                                             "pvalue", "p.adjust", "qvalue", "geneID", "Count"))
          # save detected ora categories
          count.ora.results$DISGENET[j] <- 0
          
        }
      } 
      
      ## KEGG ------------------------------------------------------------------
      
      if ("KEGG" %in% ora) {
        
        if (isEmpty(names(deg)) == FALSE) {
          
          ora_KEGG <- enrichKEGG(gene = names(deg),
                                 universe = names(universe),
                                 pAdjustMethod = pAdjustMethod,
                                 minGSSize = minGSSize,
                                 maxGSSize = maxGSSize,
                                 pvalueCutoff = pvalueCutoff)
          
          # save ora enrichment
          #In case no mapping has been performed (else):
          if (!is.null(ora_KEGG)) {
            ora.results$KEGG <- ora_KEGG
          } else {
            ora.results$KEGG <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), 
                                         c("ID", "Description", "GeneRatio", "BgRatio", 
                                           "pvalue", "p.adjust", "qvalue", "geneID", "Count"))
          }
          
          # save detected ora categories
          count.ora.results$KEGG[j] <- nrow(as.data.frame(ora_KEGG))
          
        }
        
        else {
          
          # save ora enrichment
          ora.results$KEGG <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), 
                                       c("ID", "Description", "GeneRatio", "BgRatio", 
                                         "pvalue", "p.adjust", "qvalue", "geneID", "Count"))
          # save detected ora categories
          count.ora.results$KEGG[j] <- 0 
          
        }
        
      } 
      
      ## PHENOTYPE (GSEA) ------------------------------------------------------
      
      if ("PHENOTYPE" %in% ora) {
        
        if (isEmpty(names(deg)) == FALSE) {
          
          # Download Human Phenotype Ontology
          human.phenotype.ontology <-  msigdbr(species = "Homo sapiens", 
                                               category = "C5", 
                                               subcategory = "HPO")
          
          # select entrezid and term columns to use in cluster profiler:
          # TERM2GENE and TERM2NAME
          
          hpo_term.to.gene <- human.phenotype.ontology[, c("gs_id", "entrez_gene")]
          hpo_term.to.name <- human.phenotype.ontology[, c("gs_id", "gs_name")]
          
          ora_PHENOTYPE <- enricher(gene = names(deg),
                                    universe = names(universe),
                                    pAdjustMethod = pAdjustMethod,
                                    minGSSize = minGSSize,
                                    maxGSSize = maxGSSize,
                                    pvalueCutoff = pvalueCutoff,
                                    TERM2GENE = hpo_term.to.gene,
                                    TERM2NAME = hpo_term.to.name)
          
          #In case no mapping has been performed (else):
          if (!is.null(ora_PHENOTYPE)) {
            ora.results$PHENOTYPE <- ora_PHENOTYPE
          } else {
            ora.results$PHENOTYPE <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), 
                                              c("ID", "Description", "GeneRatio", "BgRatio", 
                                                "pvalue", "p.adjust", "qvalue", "geneID", "Count"))
          }
          
          # save detected ora categories
          count.ora.results$PHENOTYPE[j] <- nrow(as.data.frame(ora_PHENOTYPE))
          
        }
        
        else {
          
          # save ora enrichment
          ora.results$PHENOTYPE <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), 
                                            c("ID", "Description", "GeneRatio", "BgRatio", 
                                              "pvalue", "p.adjust", "qvalue", "geneID", "Count"))
          # save detected ora categories
          count.ora.results$PHENOTYPE[j] <- 0 
          
        }
      }
      
      ## CELL type --------------------------------------------------------------
      
      if ("CELL" %in% ora) {
        
        if (isEmpty(names(deg)) == FALSE) {
          
          # download cell type ontology from msigdbr package
          cell.signature <-  msigdbr(species = "Homo sapiens", 
                                     category = "C8")
          
          
          cell.to.gene <- cell.signature[, c("gs_id", "entrez_gene")]
          cell.to.name <- cell.signature[, c("gs_id", "gs_name")]
          
          
          ora_CELL <- enricher(gene = names(deg),
                               universe = names(universe),
                               pAdjustMethod = pAdjustMethod,
                               minGSSize = minGSSize,
                               maxGSSize = maxGSSize,
                               pvalueCutoff = pvalueCutoff,
                               TERM2GENE = cell.to.gene,
                               TERM2NAME = cell.to.name)
          
          
          #In case no mapping has been performed (else):
          if (!is.null(ora_CELL)) {
            ora.results$CELL <- ora_CELL
          } else {
            ora.results$CELL <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), 
                                         c("ID", "Description", "GeneRatio", "BgRatio", 
                                           "pvalue", "p.adjust", "qvalue", "geneID", "Count"))
          }
          
          # save detected ora categories
          count.ora.results$CELL[j] <- nrow(as.data.frame(ora_CELL))
          
        }
        
        else {
          
          # save ora enrichment
          ora.results$CELL <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), 
                                       c("ID", "Description", "GeneRatio", "BgRatio", 
                                         "pvalue", "p.adjust", "qvalue", "geneID", "Count"))
          # save detected ora categories
          count.ora.results$CELL[j] <- 0 
          
        }
        
      }
      
      ## TISSUE type --------------------------------------------------------------
      
      
      # Tissue Human Protein Atlas (hpar package)
      if ("TISSUE" %in% ora) {
        
        if (isEmpty(names(deg)) == FALSE) {
          
          
          data(hpaNormalTissue)
          # convert ids to entrezid
          
          entrez.hpaNormalTissue <- convert.ids.function(hpaNormalTissue, 
                                                         fromType = "ENSEMBL",
                                                         toType = "ENTREZID")
          # assign an id
          entrez.hpaNormalTissue <- transform(entrez.hpaNormalTissue,id=as.numeric(factor(Tissue)))
          
          tissue.to.gene <- entrez.hpaNormalTissue[, c("id", "ENTREZID")]
          tissue.to.name <- entrez.hpaNormalTissue[, c("id", "Tissue")]
          
          
          
          
          ora_TISSUE <- enricher(gene = names(deg),
                                 universe = names(universe),
                                 pAdjustMethod = pAdjustMethod,
                                 minGSSize = minGSSize,
                                 maxGSSize = maxGSSize,
                                 pvalueCutoff = pvalueCutoff,
                                 TERM2GENE = tissue.to.gene,
                                 TERM2NAME = tissue.to.name)
          
          
          #In case no mapping has been performed (else):
          if (!is.null(ora_TISSUE)) {
            ora.results$TISSUE <- ora_TISSUE
          } else {
            ora.results$TISSUE <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), 
                                           c("ID", "Description", "GeneRatio", "BgRatio", 
                                             "pvalue", "p.adjust", "qvalue", "geneID", "Count"))
          }
          
          # save detected ora categories
          count.ora.results$TISSUE[j] <- nrow(as.data.frame(ora_TISSUE))
          
        }
        
        else {
          
          # save ora enrichment
          ora.results$TISSUE <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), 
                                         c("ID", "Description", "GeneRatio", "BgRatio", 
                                           "pvalue", "p.adjust", "qvalue", "geneID", "Count"))
          # save detected ora categories
          count.ora.results$TISSUE[j] <- 0 
          
        }
        
      }
      
      
      
      
      # ------------------------------------------------------------------------
      
      # save all ora results depending on regulation 
      all.ora.results.reg[[reg]] <- ora.results
      DEGs_fc_reg[[reg]] <- deg
      
    } # end regulation loop (j)
    
    # save all ora results depending on fold change
    all.ora.results[[fold]] <- all.ora.results.reg
    
    # save differentially expressed proteins
    DEGs_fc[[fold]] <- DEGs_fc_reg
    
    # save detected categories for each fold change and regulation
    count.df <- rbind(count.df, as.data.frame(count.ora.results))
    
  } # end fold change loop (i)
  
  # Transform dataframe into large formart using melt function
  count.df <- melt(count.df, id.vars = c("fold_change", "DEGs", "regulation"), 
                   value.name = "Count", variable.name = "Category") %>%
    arrange(fold_change)
  
  # RESULTS
  ora_results_list <- list(all.ora.results,
                           count.df,
                           DEGs_fc)
  
  return(ora_results_list)
  
} 



################################################################################
# Function 8: line.plot.function
#             performs several types of plots needed in the main function  
################################################################################

line.plot.function <- function(data, x, y = NULL, z = NULL,
                               x.title = NULL,
                               y.title = NULL,
                               main.title = NULL,
                               type = "multiple") {
  
  if (length(x) > 7) {
    angle = 45
    hjust = 1
  } else { 
    angle = 0
    hjust = 0.5
  }
  
  theme_set(theme_bw())
  
  if (type == "simple") {
    ggplot(data, 
           aes(x = factor(x, levels = x), 
               y = y, 
               group = 1)) + 
      geom_line(color = 1) + 
      geom_point(color = 1, size = 2) +
      geom_text(aes(label = y), hjust = 0, vjust = -0.5, size = 3, color = "gray43") + 
      ggtitle(main.title) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = angle, hjust = hjust)) +
      xlab(x.title) + 
      ylab(y.title) 
  }
  
  else if (type == "multiple") {
    data %>% 
      ggplot(aes_string(x = x, 
                        y = y, 
                        color = z,
                        group = z)) +
      geom_line(aes_string(linetype = z)) +
      geom_point() +
      geom_text(aes_string(label = y), 
                hjust = 0, vjust = -0.5, size = 3,
                show.legend = FALSE) +
      ggtitle(main.title) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = angle, hjust = hjust)) +
      xlab(x.title) + 
      ylab(y.title) 
  }
  else if (type == "continuous.distribution"){
    
    ggplot(data, aes_string(x=x)) + 
      geom_histogram(color="cornflowerblue", fill="lightblue") + 
      ggtitle(main.title) + ylab("Detected proteins") + xlab(x.title) +
      theme(plot.title = element_text(hjust = 0.5)) 
    
  }
  
  else if (type == "discrete.distribution"){
    
    ggplot(data, aes_string(x=x)) + 
      geom_bar(color="cornflowerblue", fill="lightblue") + 
      ggtitle(main.title) + ylab("Detected proteins") + xlab(x.title) +
      theme(plot.title = element_text(hjust = 0.5)) 
    
  }
}


################################################################################
# Function 9: PlotGoAnnotations
#             performs a barplot with the main GO annotations for each 
#             module (biological process, molecular function, cellular component)
################################################################################

PlotGoAnnotations <- function(GOObj, 
                              showCategory = 15) {
  
  # save each GO term domain categories:
  
  GO.biological.procceses <- GOObj[GOObj$GO.domain == "biological_process",
                                   c("GO.term.name", "GO.count")] %>%
    arrange(desc(GO.count)) %>% distinct(GO.term.name, .keep_all = TRUE)
  

  GO.molecular.functions <- GOObj[GOObj$GO.domain == "molecular_function",
                                  c("GO.term.name", "GO.count")] %>%
    arrange(desc(GO.count)) %>% distinct(GO.term.name, .keep_all = TRUE)
  
  
  GO.cellular.components <- GOObj[GOObj$GO.domain == "cellular_component",
                                  c("GO.term.name", "GO.count")] %>%
    arrange(desc(GO.count)) %>% distinct(GO.term.name, .keep_all = TRUE)
  

  
  
  GO.biological.procceses.plot <- ggplot(data = GO.biological.procceses[1:15,], 
                                         aes(x = reorder(GO.term.name, GO.count), 
                                             y = GO.count)) +
    geom_bar(stat = "identity", fill = 2, alpha = 0.75) + 
    scale_color_gradient(high = "blue", low = "red") +
    xlab("") + 
    ylab("Protein count") +
    ggtitle("Main Gene Ontology categories",
            subtitle = "A. Biological Processes") +
    theme_bw() +
    theme(text = element_text(size = 12, face = "bold", colour = "black"),
          axis.text.x = element_text(vjust = 2), plot.title = element_text(hjust = 0.5)) +
    coord_flip()
  
  
  GO.molecular.functions.plot <- ggplot(data = GO.molecular.functions[1:15,], 
                                        aes(x = reorder(GO.term.name, GO.count), 
                                            y = GO.count)) +
    geom_bar(stat = "identity", fill = 3, alpha = 0.75) + 
    scale_color_gradient(high = "blue", low = "red") +
    xlab("") + 
    ylab("Protein count") +
    ggtitle("",
            subtitle = "B. Molecular Functions") +
    theme_bw() +
    theme(text = element_text(size = 12, face = "bold", colour = "black"),
          axis.text.x = element_text(vjust = 2)) +
    coord_flip()
  
  
  GO.cellular.components.plot <- ggplot(data = GO.cellular.components[1:15,], 
                                        aes(x = reorder(GO.term.name, GO.count), 
                                            y = GO.count)) +
    geom_bar(stat = "identity", fill = 4, alpha = 0.75) + 
    scale_color_gradient(high = "blue", low = "red") +
    xlab("") + 
    ylab("Protein count") +
    ggtitle("", 
            subtitle = "C. Cellular Components") +
    theme_bw() +
    theme(text = element_text(size = 12, face = "bold", colour = "black"),
          axis.text.x = element_text(vjust = 2)) +
    coord_flip()
  
  
  #Plot Go Info
  ggpubr::ggarrange(GO.biological.procceses.plot, GO.molecular.functions.plot, 
            GO.cellular.components.plot, 
            nrow = 3 , ncol = 1 , align = "hv", label.x = "Protein count")

}


################################################################################
# Function 10: core.enrichment.translation
#              translate ENTREZ gene sets of the resulting enriched categories
#              to UnipProt IDs
################################################################################

core.enrichment.translation <- function(data, 
                                        col.name = "core_enrichment") {
  
  core.enriched.entrez <- lapply(data[[col.name]], function (x) strsplit(x,
                                                                         split='/', 
                                                                         fixed=TRUE))
  
  
  
  core.enriched.prot <- lapply(unlist(core.enriched.entrez,recursive = FALSE), 
                               function (x) bitr(x, 
                                                 fromType = "ENTREZID",
                                                 toType = "UNIPROT",
                                                 OrgDb = org.Hs.eg.db,
                                                 drop = TRUE))                        
  
  data[[col.name]] <- lapply(core.enriched.prot,
                             function (x) as.character(paste(x[,2], collapse = "/")))
  
  data[[col.name]] <- as.character(data[[col.name]])
  
  return(data)
  
}



################################################################################
# Function 11: main.function
#              EnrichProt main function. Performs the whole analysis:
#                 A. Differential abundance analysis
#                 B. Funcional profiling using GSEA
#                 C. Functional profiling using ORA
################################################################################

main.function <- function(data,
                          output.dir,
                          header = TRUE, 
                          sep = "\t",
                          dec = ".", 
                          remove.na = TRUE,
                          # id conversion
                          fromType, 
                          OrgDb = org.Hs.eg.db,
                          rm.duplicates = TRUE,
                          # peptides under study
                          peptides = seq(from=1, to=5, by = 1), 
                          Np.col,
                          pept.aggregate = TRUE,
                          # t.test function
                          t_test = FALSE,
                          control.cols = NULL, 
                          case.cols = NULL,
                          log.scale = TRUE,
                          fdr = 0.05,
                          # gsea function
                          gsea.test = TRUE,
                          metric = "fold_change", 
                          ont = "all", 
                          gsea = "all", 
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "fdr",
                          minGSSize = 10,
                          maxGSSize = 500,
                          # ora function
                          ora.test = TRUE,
                          FDR.name.col = "FDR",
                          regulation = "all",
                          ora = "all",
                          fc = c(1.5, 2),
                          showCategory.plots = 10) {
  
  
  ######################################
  # 1: Create EnrichProt results folder
  ######################################
  
  date <- str_replace_all(Sys.Date(), "-", ".")
  dir.create(file.path(output.dir, paste0(date, "_EnrichProt.results")), showWarnings = FALSE)
  setwd(file.path(output.dir, paste0(date, "_EnrichProt.results")))
  
  
  ####################
  # 2: Upload the data
  ####################
  
  
  prot_list <- upload.data.function(data = data,
                                    header = header, 
                                    sep = sep,
                                    dec = dec, 
                                    remove.na = remove.na)
  
  
  ##############################
  # 3: ID conversion to ENTREZID
  ##############################
  
  
  # In case protein ids are in ENTREZ format, there is no need for translation
  
  if (fromType != "ENTREZID") {
    entrez_expression <- convert.ids.function(data = prot_list, 
                                              fromType = fromType,
                                              toType = "ENTREZID",
                                              OrgDb = OrgDb,
                                              rm.duplicates = rm.duplicates)
  }   else {
    entrez_expression <- prot_list
  }
  
  
  ###############################
  # 4: Generate each protein list
  ###############################
  
  
  entrez_expression_res <- peptide.function(entrez_expression, 
                                            peptides = peptides,
                                            Np.col = Np.col)
  
  # save all protein expression lists
  entrez_expression_list <- entrez_expression_res[[1]]
  # save count data frame of the detected proteins 
  count_entrez_expression <- entrez_expression_res[[2]]
  
  uniprot_expression_res <- uniprot.peptide.function(prot_list,
                                                     peptides = peptides,
                                                     Np.col = Np.col)
  count_uniprot_expression <- uniprot_expression_res[[2]]
  
  detected.prots <- merge(count_uniprot_expression, count_entrez_expression,
                          by = "peptides")
  
  colnames(detected.prots) <- c("peptides", "UNIPROT", "ENTREZ")
  
  if (">0" %in% peptides) {
    
    detected.prots[peptides == ">0", 1] <- "All"
  }
  
  detected.prots <- melt(detected.prots, variable.name = "ID.source")
  detected.prots[, "peptides"] <- factor(detected.prots[, "peptides"],
                                         levels = unique(detected.prots[, "peptides"]))
  
  
  # Represent the number of detected proteins based on the minimum number of 
  # peptides required for the detection of these proteins:
  line.plot.function (data = detected.prots,
                      x = "peptides",
                      y = "value", 
                      z = "ID.source",
                      x.title = "Number of peptides to consider proteins as true hits",
                      y.title = "Protein count",
                      main.title = "Detected proteins according to peptide number",
                      type = "multiple")
  
  # detected proteins folder
  ggsave("summary/detected.proteins/detected.proteins.png", 
         width = 7.7, height = 5.5)
  
  
  #########################
  # 5: t-test for each list
  #########################
  
  
  num.peptides <- list()
  
  if (t_test == "TRUE") {
    
    # perform t.test for each protein expression list
    for (i in 1:length(entrez_expression_list)) {
      entrez_expression_list[[i]] <- t.test.function(entrez_expression_list[[i]],
                                                     control.cols = control.cols,
                                                     case.cols = case.cols,
                                                     log.scale = log.scale)
      
      # create statistics folder for each peptide folder
      
      if (grepl(">", peptides[i]) == TRUE) {
        
        num.peptides[i] <- paste0(eval(parse(text = substr(peptides[i], 2, str_length(peptides[i]))))+1, 
                                  ".or.more.peptides")
        
      } else if (grepl("-", peptides[i]) == TRUE) {
        
        num.peptides[i] <- paste0(str_replace(peptides[i], "-", "_"), 
                                  ".peptides")
        
      } else {
        
        if (peptides[i] == 1) {
          num.peptides[i] <- paste0(peptides[i], ".peptide")
        } else {
          num.peptides[i] <- paste0(peptides[i], ".peptides")
        }
        
      }
      
      ###
      
      dir.create(file.path(num.peptides[i], "differential.abundance"), showWarnings = FALSE)
      
      # create data  and plots folder 
      path1 <- paste0(num.peptides[i], "/differential.abundance/")
      dir.create(file.path(path1, "results"), showWarnings = FALSE)
      dir.create(file.path(path1, "plots"), showWarnings = FALSE)
      
      file.name.data <- paste0(path1, "results/")
      file.name.plots <- paste0(path1, "plots/")
      
      # save data with t-test results but previosly translating ENTREZIDs to UNIPROT:
      
      differential.abundance.file <- entrez_expression_list[[i]]
      differential.abundance.file <- convert.ids.function(differential.abundance.file,
                                                          fromType = "ENTREZID",
                                                          toType = fromType)
      
      # add GO annotation previously downloaded from BioMart-Ensembl:
      
      GO.annotation <- read.csv("../annotations/mart_export.txt", sep = "\t",
                                dec = ".")
      
      # change first column name to "UNIPROT":
      
      colnames(GO.annotation)[1] <- "UNIPROT"
      
      # annotate all identified proteins:
      
      annotated.diff.ab.file <- merge(differential.abundance.file,
                                      GO.annotation, by = "UNIPROT")
      
      # add a column "count" containing the ocurrences of each Go term name:
      
      occurrences <- table(annotated.diff.ab.file$GO.term.name)
      
      annotated.diff.ab.file <- annotated.diff.ab.file %>% add_count(GO.term.accession,
                                                                     name = "GO.count")
      

      # save annotated differentially abundant proteins 
      
      write.table(annotated.diff.ab.file, paste0(file.name.data, 
                                                      "differential.abundance.", 
                                                      num.peptides[i], 
                                                      ".csv" ), 
                  row.names = FALSE, quote = FALSE, sep = "\t")
      
      # plot most frequent GO categories for each module
      
      PlotGoAnnotations(annotated.diff.ab.file)
      ggsave(paste0(file.name.plots, "differential.abundance.", num.peptides[i], 
                    ".png"), height = 9, width = 9.3)
      
      # since annotation takes time, we perform it once with smallest fold change
      # value that will contain the larger number of proteins:
      
      # save differential abundance data with that minimum value of fold change:
      differential.abundance.file.fc <- annotated.diff.ab.file
      
      
      for (j in 1:length(fc)) {
        
        # create a specific fold change folder with its differential abundance
        fc.dir <- paste0("fold_change=", fc[j])
        dir.create(file.path(file.name.data, fc.dir), showWarnings = FALSE)
        diff.ab.dir <- paste0(file.name.data, fc.dir, "/")
        
        # save differential abundance data with that fold change:
        differential.abundance.file.fc.. <- differential.abundance.file.fc[abs(differential.abundance.file.fc[metric]) > fc[j], ]
        
        # save data:
        write.table(differential.abundance.file.fc.., 
                    paste0(diff.ab.dir, "differential.abundance.", 
                           num.peptides[i], ".fc=", fc[j], ".csv" ), 
                    row.names = FALSE, quote = FALSE, sep = "\t")
        
        # generate GO annotation plot
        PlotGoAnnotations(differential.abundance.file.fc..)
        
        ggsave(paste0(diff.ab.dir, "differential.abundance.GO.", num.peptides[i], 
                      ".fc=", fc[j], ".png" ), height = 9, width = 9.3)
        
        
      }
      
      # save distribution plots of the statistics:
      
      ## Fold change distribution
      ### define same color for each fold change
      my.color <- rep(seq(from = 2, to = length(fc)+1), each = 1, length.out = length(fc)*2) # 1 is black color
      fold <- data.frame(fc = c(fc, -fc))
      fold$`fold change` <- paste0("±", gsub("-","",fold$fc))
      line.plot.function(entrez_expression_list[[i]], x = "fold_change", type = "continuous.distribution", 
                         main.title = "Fold change distribution", x.title = "Fold change") +
        geom_vline(aes(xintercept = fc, color = `fold change`), fold, show.legend = TRUE, linetype="dashed")  
      ggsave(paste0(file.name.plots, paste0("fold_change.distribution.", num.peptides[i], ".png")), width = 8, height = 5)
      
      ## Adjusted p-value distribution
      line.plot.function(entrez_expression_list[[i]], x = "adjusted.p.value", type = "continuous.distribution", 
                         main.title = "Adjusted p-value distribution",
                         x.title = "Adjusted p-value") +
        geom_vline(aes(xintercept = pvalueCutoff), color = 2, linetype = "dashed", size = 0.5) +
        geom_text(aes(x = pvalueCutoff+0.005, y = Inf, hjust = -0.1, vjust = 1.8, 
                      label = as.data.frame(paste0("p-value = ", pvalueCutoff))), 
                  colour = 2, angle = 0, text = element_text(size=0.5))
      ggsave(paste0(file.name.plots, paste0("adjusted.p.value.distribution.", num.peptides[i], ".png")), width = 8, height = 5)
      
      ## t-statistic distribution
      line.plot.function(entrez_expression_list[[i]], x = "statistic", type = "continuous.distribution", 
                         main.title = "t-statistic distribution",
                         x.title = "t-statistic") 
      ggsave(paste0(file.name.plots, paste0("t_statistic.distribution.", num.peptides[i], ".png")), width = 8, height = 5)
      
      
      if (grepl(">", peptides[i]) == TRUE) {
        
        
        # peptide number distribution
        line.plot.function(entrez_expression_list[[i]], x = Np.col, type = "discrete.distribution", 
                           main.title = "Detected proteins distribution by peptide number",
                           x.title = "Number of peptides")
        ggsave(paste0(file.name.plots, paste0("peptides.distribution.", num.peptides[i], ".png")), width = 8, height = 5)
        
        maximum <- round(max(entrez_expression_list[[i]][[Np.col]])/6, 0)
        
        # zoom peptide number distribution
        line.plot.function(entrez_expression_list[[i]], x = Np.col, type = "discrete.distribution",
                           main.title = "Detected proteins distribution by peptide number",
                           x.title = "Number of peptides") + 
          scale_x_continuous(limits = c(i-1, maximum+1), 
                             breaks = seq(i, maximum, by = 2), 
                             minor_breaks = 1) 
        ggsave(paste0(file.name.plots, paste0("zoom.peptides.distribution.", num.peptides[i], ".png")), width = 8, height = 5)
        
      }
      
      
      # Volcano plot:
      
      
      if (log.scale == "TRUE") {
        
        abs_log <- function(x) {
          x[x==0] <- 1
          si <- sign(x)
          si * log2(si*x)
        }
        
        xlab.name <- "log2(fold change)"
        fold.plot <- fold
        fold.plot$fc <- abs_log(fold$fc)
          
        } else {
        xlab.name <- "fold change"
        fold.plot <- fold
      }
      ggplot(entrez_expression_list[[i]], 
             aes(x = fold_change, y = -log10(adjusted.p.value))) + 
        ylab("-log10(p-value)") + xlab(xlab.name) +
        geom_point() +
        geom_hline(yintercept = -log10(pvalueCutoff), col = "red") +
        geom_vline(aes(xintercept = fc, color = `fold change`), 
                   fold.plot, show.legend = TRUE, linetype = "dashed") 
      
      ggsave(paste0(file.name.plots, paste0("volcano.plot.", num.peptides[i], ".png")))
    }
    
  } # end t-test function
  
  
  ##################
  # 6: GSEA analysis
  ##################
  
  
  if (gsea.test == "TRUE") {
    
    # save GSEA results
    gsea_results_list <- list()
    count.gsea.result <- data.frame()
    
    # 6.1. Protein list loop (i)
    
    
    for (i in 1:length(entrez_expression_list)) {
      
      gsea.result <- gsea.function(entrez_expression_list[[i]],
                                   metric = metric, 
                                   ont = ont, 
                                   gsea = gsea, 
                                   pvalueCutoff = pvalueCutoff,
                                   pAdjustMethod = pAdjustMethod,
                                   minGSSize = minGSSize,
                                   maxGSSize = maxGSSize)
      
      gsea.result[[2]]$Peptides <- rep(peptides[i], length = nrow(gsea.result[[2]]))
      gsea_results_list[[i]] <- gsea.result
      
      if (grepl(">", peptides[i]) == TRUE) {
        
        num.peptides[i] <- paste0(eval(parse(text = substr(peptides[i], 2, str_length(peptides[i]))))+1, 
                                  ".or.more.peptides")
        
      } else if (grepl("-", peptides[i]) == TRUE) {
        
        num.peptides[i] <- paste0(str_replace(peptides[i], "-", "_"), 
                                  ".peptides")
        
      } else {
        
        if (peptides[i] == 1) {
          num.peptides[i] <- paste0(peptides[i], ".peptide")
        } else {
          num.peptides[i] <- paste0(peptides[i], ".peptides")
        }
        
      }
      
      
      
      dir.create(file.path(num.peptides[i], "GSEA.results"), showWarnings = FALSE)
      
      
      # 6.2. Ontology loop (j)
      
      
      # save gsea results for each peptide for each ontology
      for (j in 1:length(names(gsea_results_list[[i]][[1]]))) {
        path <- paste0(num.peptides[i], 
                       "/GSEA.results/")
        dest.file <- names(gsea_results_list[[i]][[1]][j])
        dir.create(file.path(path, dest.file), showWarnings = FALSE)
        if (names(gsea_results_list[[i]][[1]][j]) == "GO") {
          
          
          # 6.2.1. Specific Gene Ontology (k)
          
          
          # specific GO ontology loop:
          for (k in 1:length(names(gsea_results_list[[i]][[1]][[j]]))) {
            path <- paste0(num.peptides[i], 
                           "/GSEA.results/GO")
            dest.file <- names(gsea_results_list[[i]][[1]][[j]][k])
            dir.create(file.path(path, dest.file), showWarnings = FALSE)
            # create two folders:
            path2 <- paste0(path, "/", dest.file)
            dest.file2 <- "enriched.categories"
            dest.file3 <- "plots"
            dir.create(file.path(path2, dest.file2), showWarnings = FALSE)
            dir.create(file.path(path2, dest.file3), showWarnings = FALSE)
            
            file.name.categories <- paste0(num.peptides[i], 
                                           "/GSEA.results/GO/", 
                                           dest.file, "/", dest.file2, "/",
                                           names(gsea_results_list[[i]][[1]][j]),
                                           "(",
                                           names(gsea_results_list[[i]][[1]][[j]][k]),
                                           ")")
            
            file.name.plots <- paste0(num.peptides[i], 
                                      "/GSEA.results/GO/", 
                                      dest.file, "/", dest.file3, "/",
                                      names(gsea_results_list[[i]][[1]][j]),
                                      "(",
                                      names(gsea_results_list[[i]][[1]][[j]][k]),
                                      ")")
            title <- paste0(names(gsea_results_list[[i]][[1]][j]),
                            "(",
                            names(gsea_results_list[[i]][[1]][[j]][k]),
                            ") enriched categories")
            
            
            # 6.2.1.1. Specific Gene Ontology (k) categories:
            
            
            # Use NES (Normalized Enrichment Score) as a criteria to define
            # down and up-regulated enriched categories
            
            all.enriched.categories <- as.data.frame(gsea_results_list[[i]][[1]][[j]][[k]])
            all.enriched.categories <- core.enrichment.translation(all.enriched.categories)
            down.regulated.enriched.categories <- all.enriched.categories[all.enriched.categories$NES < 0,]
            up.regulated.enriched.categories <- all.enriched.categories[all.enriched.categories$NES > 0,]
            
            
            write.table(all.enriched.categories,
                        paste0(file.name.categories, ".gsea.UP.and.DOWN.regulated.enriched.categories.", num.peptides[i] , ".csv"), 
                        row.names=FALSE, quote=FALSE, sep = "\t")
            write.table(down.regulated.enriched.categories,
                        paste0(file.name.categories, ".gsea.DOWN.regulated.enriched.categories.", num.peptides[i] , ".csv"), 
                        row.names=FALSE, quote=FALSE, sep = "\t")
            write.table(up.regulated.enriched.categories,
                        paste0(file.name.categories, ".gsea.UP.regulated.enriched.categories.", num.peptides[i] , ".csv"), 
                        row.names=FALSE, quote=FALSE, sep = "\t")
            
            # Reduced GO terms using REVIGO (C=0.5)
            
            if (nrow(all.enriched.categories > 10)) {
              
              revigo.terms <- all.enriched.categories %>% arrange(p.adjust, desc(NES))
              
              simMatrix <- calculateSimMatrix(revigo.terms$ID,
                                              orgdb = "org.Hs.eg.db",
                                              ont = ont[k],
                                              method = "Rel")
              scores <- setNames(-log10(revigo.terms$qvalue), revigo.terms$ID)
              reducedTerms <- reduceSimMatrix(simMatrix,
                                              scores,
                                              threshold = 0.5,
                                              orgdb = "org.Hs.eg.db")
              
              # save reduced GO terms:
              write.table(reducedTerms,
                          paste0(file.name.categories, ".gsea.reduced.enriched.categories.", num.peptides[i] , ".csv"), 
                          row.names = FALSE, quote = FALSE, sep = "\t")
              
              
              if (!nrow(reducedTerms) < 50) {
                
                pdf(file=paste0(file.name.plots, ".reduced.terms.pdf"), width = 10, height = 6)
                
                treemap::treemap(reducedTerms[1:50,], title = "",
                                 index = c("parentTerm", "term"),
                                 fontsize.labels=c(12, 10),
                                 fontcolor.labels=c("white","black"),
                                 vSize = "score",
                                 fontface.labels=c(2,1),
                                 bg.labels=c("transparent"),
                                 type = "index",
                                 overlap.labels=1)
                dev.off()
                
              } else {
                
                pdf(file=paste0(file.name.plots, ".reduced.terms.pdf"), width = 10, height = 6)
                
                treemap::treemap(reducedTerms, title = "",
                                 index = c("parentTerm", "term"),
                                 fontsize.labels=c(12, 10),
                                 fontcolor.labels=c("white","black"),
                                 vSize = "score",
                                 fontface.labels=c(2,1),
                                 bg.labels=c("transparent"),
                                 type = "index",
                                 overlap.labels=1)
                dev.off()
                
              }
              
            }
          
            
            
            # 6.2.1.2. Specific Gene Ontology (k) plots:
            
            
            # GSEA results distribution plot
            a <- gsea_results_list[[i]][[1]][[j]][[k]]
            a.df <- data.frame(a)
            
            
            
            # in case a = 0, make it have 0 rows
            if (!nrow(a.df) == 0) {
              if (a.df[1,1] == 0) {
                a <- NULL
                a.df <- data.frame()
              }
            }
            
            if (!nrow(a.df) == 0) {
              
              edox <- setReadable(a, 'org.Hs.eg.db', 'ENTREZID')
              
              if (length(a.df$Description) >= 5) {
                l <- gseaplot2(a, geneSetID = 1:5, ES_geom = "line") 
                d <- gseaplot2(a, geneSetID = 1:5, ES_geom = "dot") 
              } else if (length(a.df$Description) <= 5) {
                l <- gseaplot2(a, geneSetID = 1:length(a@result$Description), ES_geom = "line") 
                d <- gseaplot2(a, geneSetID = 1:length(a@result$Description), ES_geom = "dot") 
              }
              ggsave(paste0(file.name.plots, ".line.distribution.png"), l, 
                     width = 10, height = 6)
              ggsave(paste0(file.name.plots, ".dot.distribution.png"), d, 
                     width = 10, height = 6)
              
              # ridgeplot
              ridgeplot(a, showCategory = showCategory.plots) +
                geom_density_ridges(jittered_points = TRUE,
                                    position = position_points_jitter(width = 0.05, height = 0),
                                    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
                theme_ridges() + xlab("Fold change") + theme(axis.title.x = element_text(hjust = 0.5))
              ggsave(paste0(file.name.plots, ".ridgeplot.png"), 
                     width = 10, height = 6)
              
              if ((length(a.df$Description) != 1) & (a@result$setSize[1] != 1)) {
                
                # cnetplot
                cnetplot(edox,
                         foldChange = edox@geneList,
                         circular = TRUE, colorEdge = TRUE, 
                         cex_label_gene = 0.6, cex_label_category = 0.7,
                         cex_category = 0.9,
                         showCategory = showCategory.plots) 
                ggsave(paste0(file.name.plots, ".cnetplot.png"), 
                       width = 10, height = 8)
                
                #cnetplot2
                cnetplot(edox, 
                         node_label="all",
                         foldChange = edox@geneList,
                         colorEdge = TRUE,
                         cex_label_category = 0.8,color_category = "yellow",
                         cex_label_gene = 0.5, 
                         showCategory = showCategory.plots,) 
                ggsave(paste0(file.name.plots, ".cnetplot2.png"), 
                       width = 11, height = 6.7)
                
              }
              
              # upsetplot: at least three categories to compare
              if (nrow(a.df) > 3) {
                upsetplot(a, showCategory = showCategory.plots) + 
                  ylab("Fold change distribution") + 
                  xlab("Gene sets overlap")
                ggsave(paste0(file.name.plots, ".upsetplot.png"), 
                       width = 12, height = 8)
              }
              # heatplot 
              heatplot(edox, foldChange = edox@geneList, showCategory = showCategory.plots) + 
                coord_equal()
              ggsave(paste0(file.name.plots, ".heatplot.png"), 
                     width = 10, height = 4)
              
              # barplot
              ggplot(a, showCategory = showCategory.plots, 
                     aes(NES, fct_reorder(Description, NES), fill = p.adjust)) + 
                geom_bar(stat = "identity", alpha = 1) + xlab("Normalized Enrichment Score") + 
                ylab(NULL) + scale_fill_gradient(low = "red", high = "blue") +
                ggtitle(title) +
                theme_ridges() + theme(axis.title.x = element_text(hjust = 0.5),
                                       plot.title = element_text(hjust = 0.5))
              ggsave(paste0(file.name.plots, ".barplot.png"), width = 10, height = 6)
              
              # UP vs DOWN dotplot:
              a.df$regulation <-  "Up-regulated"
              a.df$regulation[a.df$NES < 0] <- "Down-regulated" 
              
              if (nrow(a.df) > showCategory.plots) {
                
                ggplot(a.df[1:showCategory.plots,], aes(x = setSize, y = fct_reorder(Description, -p.adjust))) + 
                  geom_point(aes(size = setSize, color = p.adjust), alpha = 0.7) +
                  theme_bw(base_size = 14) +
                  scale_colour_gradient(low = "red", high = "blue") +
                  ylab(NULL) +
                  ggtitle(title)  + theme(plot.title = element_text(hjust = 0.5)) +
                  facet_grid(.~regulation) +
                  scale_size(range = c(2, 14))
                
              } else {
                
                ggplot(a.df, aes(x = setSize, y = fct_reorder(Description, -p.adjust))) + 
                  geom_point(aes(size = setSize, color = p.adjust), alpha = 0.7) +
                  theme_bw(base_size = 14) +
                  scale_colour_gradient(low = "red", high = "blue") +
                  ylab(NULL) +
                  ggtitle(title)  + theme(plot.title = element_text(hjust = 0.5)) + 
                  facet_grid(.~regulation) +
                  scale_size(range = c(3, 12))
                
              }
              
              ggsave(paste0(file.name.plots, ".dotplot.png"), 
                     width = 10, height = 8)
              
            }
            
          }
          
        } # end GO ontology
        
        
        # 6.2. Ontology loop (j)
        
        
        else {
          
          path2 <- paste0(path, "/", dest.file)
          dest.file2 <- "enriched.categories"
          dest.file3 <- "plots"
          dir.create(file.path(path2, dest.file2), showWarnings = FALSE)
          dir.create(file.path(path2, dest.file3), showWarnings = FALSE)
          
          file.name.categories <- paste0(num.peptides[i], 
                                         "/GSEA.results/", 
                                         names(gsea_results_list[[i]][[1]][j]),
                                         "/", dest.file2, "/",
                                         names(gsea_results_list[[i]][[1]][j]))
          
          file.name.plots <- paste0(num.peptides[i], 
                                    "/GSEA.results/", 
                                    names(gsea_results_list[[i]][[1]][j]),
                                    "/", dest.file3, "/",
                                    names(gsea_results_list[[i]][[1]][j]))
          
          title <- paste0(names(gsea_results_list[[i]][[1]][j]),
                          " enriched categories")
          
          
          # 6.2.1. Ontology (k) categories:
          
          
          # Use NES (Normalized Enrichment Score) as a criteria to define
          # down and up-regulated enriched categories
          
          all.enriched.categories <- as.data.frame(gsea_results_list[[i]][[1]][[j]])
          all.enriched.categories <- core.enrichment.translation(all.enriched.categories) 
          
          
          down.regulated.enriched.categories <- all.enriched.categories[all.enriched.categories$NES < 0,]
          up.regulated.enriched.categories <- all.enriched.categories[all.enriched.categories$NES > 0,]
          
          write.table(all.enriched.categories,
                      paste0(file.name.categories, ".gsea.UP.and.DOWN.regulated.enriched.categories.", num.peptides[i] , ".csv"), 
                      row.names=FALSE, quote=FALSE, sep = "\t")
          write.table(down.regulated.enriched.categories,
                      paste0(file.name.categories, ".gsea.DOWN.regulated.enriched.categories.", num.peptides[i] , ".csv"), 
                      row.names=FALSE, quote=FALSE, sep = "\t")
          write.table(up.regulated.enriched.categories,
                      paste0(file.name.categories, ".gsea.UP.regulated.enriched.categories.", num.peptides[i] , ".csv"), 
                      row.names=FALSE, quote=FALSE, sep = "\t")
          
          
          
          # 6.2.2. Ontology (k) plots:
          
          
          # GSEA result distribution plot
          a <- gsea_results_list[[i]][[1]][[j]]
          a.df <- data.frame(a)
          
          
          # in case a = 0, make it have 0 rows
          if (!nrow(a.df) == 0) {
            if (a.df[1,1] == 0) {
              a <- NULL
              a.df <- data.frame()
            }
          }
          
          if (!nrow(a.df) == 0) {
            
            edox <- setReadable(a, 'org.Hs.eg.db', "ENTREZID")
            
            if (length(a.df$Description) >= 5) {
              l <- gseaplot2(a, geneSetID = 1:5, ES_geom = "line")
              d <- gseaplot2(a, geneSetID = 1:5, ES_geom = "dot")
            } else if (length(a.df$Description) <= 5) {
              l <- gseaplot2(a, geneSetID = 1:length(a@result$Description), ES_geom = "line")
              d <- gseaplot2(a, geneSetID = 1:length(a@result$Description), ES_geom = "dot")
            }
            ggsave(paste0(file.name.plots, ".line.distribution.png"), l, 
                   width = 10, height = 6)
            ggsave(paste0(file.name.plots, ".dot.distribution.png"), d, 
                   width = 10, height = 6)
            
            # ridgeplot
            ridgeplot(a, showCategory = showCategory.plots) +
              geom_density_ridges(jittered_points = TRUE,
                                  position = position_points_jitter(width = 0.05, height = 0),
                                  point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
              theme_ridges() + xlab("Fold change") + theme(axis.title.x = element_text(hjust = 0.5))
            ggsave(paste0(file.name.plots, ".ridgeplot.png"), 
                   width = 10, height = 6)
            
            if ((length(a.df$Description) != 1) & (a@result$setSize[1] != 1)) {
              
              # cnetplot
              
              cnetplot(edox,
                       foldChange = edox@geneList,
                       circular = TRUE, colorEdge = TRUE, 
                       cex_label_gene = 0.6, cex_label_category = 0.7,
                       cex_category = 0.9, showCategory = showCategory.plots)
              ggsave(paste0(file.name.plots, ".cnetplot.png"), 
                     width = 10, height = 8)
              
              #cnetplot2
              cnetplot(edox, 
                       node_label="all",
                       foldChange = edox@geneList,
                       colorEdge = TRUE,
                       cex_label_category = 0.8,
                       cex_label_gene = 0.6,
                       showCategory = showCategory.plots) 
              ggsave(paste0(file.name.plots, ".cnetplot2.png"), 
                     width = 10, height = 8)
              
            }
            
            # upsetplot: at least three categories to compare
            if (nrow(a.df) > 3) {
              upset <- upsetplot(a, showCategory = showCategory.plots) + 
                ylab("Fold change distribution") + 
                xlab("Protein sets overlap")
              ggsave(paste0(file.name.plots, ".upsetplot.png"), upset, 
                     width = 12, height = 8)
            }
            
            # heatplot 
            heatplot(edox, foldChange = edox@geneList, showCategory = showCategory.plots) + 
              coord_equal()
            ggsave(paste0(file.name.plots, ".heatplot.png"), 
                   width = 10, height = 4)
            
            # barplot
            ggplot(a, showCategory = showCategory.plots, 
                   aes(NES, fct_reorder(Description, NES), fill = p.adjust)) + 
              geom_bar(stat = "identity", alpha = 1) + xlab("Normalized Enrichment Score") + 
              ylab(NULL) + scale_fill_gradient(low = "red", high = "blue") +
              ggtitle(title) +
              theme_ridges() + theme(axis.title.x = element_text(hjust = 0.5),
                                     plot.title = element_text(hjust = 0.5))
            ggsave(paste0(file.name.plots, ".barplot.png"), width = 10, height = 6)
            
            
            # UP vs DOWN dotplot:
            a.df$regulation <- "Up-regulated"
            a.df$regulation[a.df$NES < 0] <- "Down-regulated"
            
            if (nrow(a.df) > showCategory.plots) {
              
              ggplot(a.df[1:showCategory.plots,], aes(x = setSize, y = fct_reorder(Description, -p.adjust))) + 
                geom_point(aes(size = setSize, color = p.adjust), alpha = 0.7) +
                theme_bw(base_size = 14) +
                scale_colour_gradient(low = "red", high = "blue") +
                ylab(NULL) +
                ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)) + 
                facet_grid(.~regulation) +
                scale_size(range = c(2, 14))
              
            } else {
              
              ggplot(a.df, aes(x = setSize, y = fct_reorder(Description, -p.adjust))) + 
                geom_point(aes(size = setSize, color = p.adjust), alpha = 0.7) +
                theme_bw(base_size = 14) +
                scale_colour_gradient(low = "red", high = "blue") +
                ylab(NULL) +
                ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)) + 
                facet_grid(.~regulation) +
                scale_size(range = c(2, 14))
              
            }
            
            ggsave(paste0(file.name.plots, ".dotplot.png"), 
                   width = 10, height = 8)
            
          }
        }
      }
      
      # 6.3. Summary results
      # Join all count dataframes for each peptide entrez expression:
      count.gsea.result <- rbind(count.gsea.result,
                                 as.data.frame(gsea_results_list[[i]][[2]]))
      
    } # end peptide loop (i)
    
    
    
    
    # create gsea results in summary folder
    dir.create(file.path("summary", "gsea.enriched.categories"), showWarnings = FALSE)
    file.gsea.name <- paste0("summary/gsea.enriched.categories")
    
    write.table(count.gsea.result, paste0(file.gsea.name, "/gsea.enriched.categories.csv"), 
                row.names=FALSE, quote=FALSE, sep = "\t")
    
    count.gsea.result$Peptides <- factor(count.gsea.result$Peptides, levels = unique(count.gsea.result$Peptides))
    
    line.plot.function(data = count.gsea.result, 
                       x = "Peptides",
                       y = "Count",
                       z = "Category",
                       x.title = "Number of peptides",
                       y.title = "Detected categories",
                       main.title = "Detected categories according to the number of peptides",
                       type = "multiple")
    ggsave(paste0(file.gsea.name, "/gsea.enriched.categories.png"))
    
    
    # 6.4. Bubble chart
    
    # Save enriched categories for each ontology:
    
    # Ontology loop (i)
    if (ont == "all") {
      ont <- c("BP", "MF", "CC")
    }
    
    if ("GO" %in% gsea) {
      go.names <- paste0("GO\\(", ont, ")")
      gsea.names <- c(go.names, gsea[gsea != "GO"])
    } else {
      gsea.names <- gsea
    }
    
    bubble.cols <- c("ID", "Description", "setSize", "NES", "p.adjust")
    
    for (i in 1:length(gsea.names)) {
      
      all.files <- list.files(getwd(), recursive = TRUE)
      up.files <- grep(paste0(gsea.names[i], ".gsea.UP.regulated.enriched.categories.*peptide.*.csv"), 
                       all.files, value = TRUE)
      down.files <- grep(paste0(gsea.names[i], ".gsea.DOWN.regulated.enriched.categories.*peptide.*.csv"), 
                         all.files, value = TRUE)
      up.and.down.files <- grep(paste0(gsea.names[i], ".gsea.UP.and.DOWN.regulated.enriched.categories.*peptide.*.csv"), 
                                all.files, value = TRUE)
      
      # Peptide loop (j): we want to merge all categories for the different peptides
      
      UP.regulated.enriched.categories <- data.frame()
      DOWN.regulated.enriched.categories <- data.frame()
      UP.and.DOWN.regulated.enriched.categories <- data.frame()
      
      for (j in 1:length(peptides)) {
        
        # add a column indicating the peptide
        
        pept.up.enriched.cat <- read.csv(up.files[j], header = TRUE, sep = "\t", dec = ".")[, bubble.cols]
        pept.down.enriched.cat <- read.csv(down.files[j], header = TRUE, sep = "\t", dec = ".")[, bubble.cols]
        pept.up.and.down.enriched.cat <- read.csv(up.and.down.files[j], header = TRUE, sep = "\t", dec = ".")[, bubble.cols] 
        
        
        
        
        if (nrow(pept.up.enriched.cat) == 0) {
          pept.up.enriched.cat[nrow(pept.up.enriched.cat)+1,] <- NA
        }
        pept.up.enriched.cat$peptides <- factor(peptides[j])
        
        if (nrow(pept.down.enriched.cat) == 0) {
          pept.down.enriched.cat[nrow(pept.down.enriched.cat)+1,] <- NA
        }
        pept.down.enriched.cat$peptides <- factor(peptides[j])
        
        if (nrow(pept.up.and.down.enriched.cat) == 0) {
          pept.up.and.down.enriched.cat[nrow(pept.up.and.down.enriched.cat)+1,] <- NA
        }
        pept.up.and.down.enriched.cat$peptides <- factor(peptides[j])
        
        
        
        
        UP.regulated.enriched.categories <- rbind(UP.regulated.enriched.categories, pept.up.enriched.cat)
        DOWN.regulated.enriched.categories <- rbind(DOWN.regulated.enriched.categories, pept.down.enriched.cat)
        UP.and.DOWN.regulated.enriched.categories <- rbind(UP.and.DOWN.regulated.enriched.categories, pept.up.and.down.enriched.cat)
        
      } # end peptide loop (j)
      
      # create bubble chart for each regulation and save in summary folder 
      enriched.cat <- c("UP.regulated.enriched.categories", "DOWN.regulated.enriched.categories", 
                        "UP.and.DOWN.regulated.enriched.categories")
      
      enriched.cat.reg.names <- c("UP.regulated", "DOWN.regulated", "UP.and.DOWN regulated")
      
      # Regulation loop (r)
      
      for (r in 1:length(enriched.cat)) {
        
        if (nrow(get(enriched.cat[r])) != 0) {
          
          # arrange in ascending order by p.adjust:
          enriched.matrix <- get(enriched.cat[r]) %>% arrange(p.adjust) 
          
          # as some categories are enriched in different protein lists
          # we make sure we plot unique categories using showCategory.plots
          
          categories.subset.to.plot <- unique(enriched.matrix$Description)[1:showCategory.plots]
          
          if (paste0("GO\\(", ont[i], ")") %in% gsea.names[i]) {
            
            gsea.names[i] <- paste0(paste0("GO(", ont, ")"))
            
          }
          
          # do not plot in case of no enriched categories are found:
          
          if (nrow(enriched.matrix[is.na(enriched.matrix[, 1]),]) != nrow(enriched.matrix)) {
            
            x_label_values = levels(enriched.matrix$peptides)
            
            enriched.matrix <- enriched.matrix[complete.cases(enriched.matrix), ]
            
            enriched.matrix[enriched.matrix$Description %in% categories.subset.to.plot,] %>% 
              #mutate(Average_GeneRatio = factor(Average_GeneRatio, Average_GeneRatio)) %>%
              ggplot(aes(x = peptides, 
                         y = factor(Description, levels = rev(unique(Description))), 
                         size = setSize, 
                         color = p.adjust)) +
              geom_point(alpha = 0.7) +
              #Adjust the range of points size:
              scale_size(name = "setSize", breaks = function(x) unique(floor(pretty(seq(1, (max(x) + 1) * 1.1)))),
                         range = c(2, 14)) +
              scale_color_gradient(low = "red", high = "blue") + 
              theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=0)) +
              # drop makes sure all levels are plot despite there are no values for them
              scale_x_discrete(position = "top", drop = FALSE,
                               breaks = x_label_values,
                               labels = x_label_values)  +
              xlab("Detected proteins") +
              ylab(paste(gsea.names[i], enriched.cat.reg.names[r], "enriched categories")) 
              
              dir.create(file.path("summary/gsea.enriched.categories", "bubble.chart"), showWarnings = FALSE)
            
            ggsave(paste0("summary/gsea.enriched.categories/bubble.chart/", gsea.names[i], ".", enriched.cat[r], ".png"))
            
          } # in case no enriched categories have been found
          
          
        } # end condition
        
        
      } # end regulation loop (r)
      
      
    } # end ontology loop (i)
    
    
  } # end GSEA analysis
  
  
  ##################
  # 7: ORA analysis
  ##################
  
  
  if (ora.test == "TRUE") {
    # save ORA results
    ora_results_list <- list()
    count.ora.result <- data.frame()
    
    
    # 7.1. Protein list loop (i)
    
    
    for (i in 1:length(entrez_expression_list)) {
      
      ora.result <- ora.function(entrez_expression_list[[i]], 
                                 metric = metric,
                                 ont = ont,
                                 FDR.name.col = FDR.name.col,
                                 ora = ora,
                                 regulation = regulation,
                                 OrgDb = OrgDb,
                                 fdr = fdr,
                                 pvalueCutoff = pvalueCutoff, 
                                 pAdjustMethod = pAdjustMethod,
                                 minGSSize = minGSSize,
                                 maxGSSize = maxGSSize,
                                 fc = fc)
      
      # add peptides used to detect the proteins to count data frame
      ora.result[[2]]$Peptides <- rep(peptides[i], length = nrow(ora.result[[2]]))
      
      # add ora results to its list:
      ora_results_list[[i]] <- ora.result
      
      if (grepl(">", peptides[i]) == TRUE) {
        
        num.peptides[i] <- paste0(eval(parse(text = substr(peptides[i], 2, str_length(peptides[i]))))+1, 
                                  ".or.more.peptides")
        
      } else if (grepl("-", peptides[i]) == TRUE) {
        
        num.peptides[i] <- paste0(str_replace(peptides[i], "-", "_"), 
                                  ".peptides")
        
      } else {
        
        if (peptides[i] == 1) {
          num.peptides[i] <- paste0(peptides[i], ".peptide")
        } else {
          num.peptides[i] <- paste0(peptides[i], ".peptides")
        }
        
      }
      
      
      dir.create(file.path(num.peptides[i], "ORA.results"), showWarnings = FALSE)
      
      
      # 7.2. Fold change loop (j) 
      
      
      for (j in 1:length(ora_results_list[[i]][[1]])) {
        
        path <- paste0(num.peptides[i], "/ORA.results")
        dest.file <- names(ora_results_list[[i]][[1]][j])
        
        # create fold change folders
        # "≥ 4 peptides/ORA.results/fold_change=1.5"
        
        dir.create(file.path(path, dest.file), showWarnings = FALSE)
        
        # create differentially expressed genes folder
        # "≥ 4 peptides/ORA.results/fold_change=1.5/differentially.expressed.proteins"
        
        pathh <- paste0(path, "/", dest.file)
        dest.filee <- "differentially.expressed.proteins"
        dir.create(file.path(pathh, dest.filee), showWarnings = FALSE)
        
        
        # 7.3. Regulation loop (k)
        
        
        for (k in 1:length(names(ora_results_list[[i]][[1]][[j]]))) {
          
          # differentially expressed proteins for each regulation:
          differentially.expressed.prot <- data.frame(ora_results_list[[i]][[3]][[j]][k])
          differentially.expressed.prot$ENTREZID <- rownames(differentially.expressed.prot)
          differentially.expressed.prot <- differentially.expressed.prot[, c(2, 1)]
          colnames(differentially.expressed.prot)[2] <- "fold_change"
          
          # save this file with UNIPROT ids:
          differentially.expressed.prot.uniprot <- convert.ids.function(differentially.expressed.prot,
                                                                        fromType = "ENTREZID",
                                                                        toType = "UNIPROT") %>%
            arrange(desc(fold_change))
          
          
          
          # create DEPs folder for each folder:
          path.deps <- paste0(pathh, "/", dest.filee)
          dest.deps <- names(ora_results_list[[i]][[3]][[j]][k])
          dir.create(file.path(path.deps, dest.deps), showWarnings = FALSE)
          
          # save differentially expressed proteins for each regulation:
          write.table(differentially.expressed.prot.uniprot, paste0(path.deps, "/", dest.deps, "/", dest.deps,
                                                                    ".differentially.expressed.proteins.", dest.file, ".csv"), 
                      row.names=FALSE, quote=FALSE, sep = "\t")
          
          
          # 7.4. Ontology loop (l)
          
          
          for (l in 1:length(names(ora_results_list[[i]][[1]][[j]][[k]]))) {
            
            # create ontology folder:
            dest.ontology <- names(ora_results_list[[i]][[1]][[j]][[k]][l])
            dir.create(file.path(pathh, dest.ontology), showWarnings = FALSE)
            
            # in case it is GO, go through all specific ontologies:
            if (names(ora_results_list[[i]][[1]][[j]][[k]][l]) == "GO") {
              
              # 7.5. Specific Gene Ontology loop (m)
              
              for (m in 1:length(names(ora_results_list[[i]][[1]][[j]][[k]][[l]]))) {
                
                # create specific GO ontologies folders:
                path.go.ont <- paste0(pathh, "/", dest.ontology)
                dest.go.ont <- names(ora_results_list[[i]][[1]][[j]][[k]][[l]][m])
                dir.create(file.path(path.go.ont, dest.go.ont), showWarnings = FALSE)
                
                # create regulation folder
                path.go.reg <- paste0(path.go.ont, "/", dest.go.ont)
                dest.go.reg <- names(ora_results_list[[i]][[3]][[j]][k])
                dir.create(file.path(path.go.reg, dest.go.reg), showWarnings = FALSE)
                
                # create two folders inside: enriched.categories and plots
                path.go.reg.2 <- paste0(path.go.reg, "/", dest.go.reg)
                dest.go.reg.2 <- "enriched.categories"
                dest.go.reg.3 <- "plots"
                dir.create(file.path(path.go.reg.2, dest.go.reg.2), showWarnings = FALSE)
                dir.create(file.path(path.go.reg.2, dest.go.reg.3), showWarnings = FALSE)
                
                # define a general title for plots and file names:
                file.name.categories <- paste0(path.go.reg.2, "/", dest.go.reg.2, "/", 
                                               dest.ontology, "(", dest.go.ont, ").ora.", dest.deps) 
                

                
                file.name.plots <- paste0(path.go.reg.2, "/", dest.go.reg.3, "/",
                                          dest.ontology, "(", dest.go.ont, ")")
                title <- paste0(dest.ontology, "(", dest.go.ont, ") enriched categories")
                
                
                # 7.5.1. Specific Gene Ontology (m) categories
                
                
                enriched.categories <- as.data.frame(ora_results_list[[i]][[1]][[j]][[k]][[l]][[m]])
                colnames(enriched.categories)[which(names(enriched.categories) == "geneID")] <- "proteinID"
                enriched.categories <- core.enrichment.translation(enriched.categories,
                                                                   col.name = "proteinID")
                
                write.table(enriched.categories,
                            paste0(file.name.categories, ".enriched.categories.", 
                                   num.peptides[i], ".", dest.file, ".csv"), 
                            row.names=FALSE, quote=FALSE, sep = "\t")
                
                # Reduced GO terms using REVIGO (C=0.5)
                
                if (nrow(enriched.categories) > 10) {
                  
                  revigo.terms <- enriched.categories %>% arrange(p.adjust)
                  
                  simMatrix <- calculateSimMatrix(revigo.terms$ID,
                                                  orgdb = "org.Hs.eg.db",
                                                  ont = ont[m],
                                                  method = "Rel")
                  scores <- setNames(-log10(revigo.terms$qvalue), revigo.terms$ID)
                  reducedTerms <- reduceSimMatrix(simMatrix,
                                                  scores,
                                                  threshold = 0.5,
                                                  orgdb = "org.Hs.eg.db")
                  
                  # save reduced GO terms:
                  write.table(reducedTerms,
                              paste0(file.name.categories, ".reduced.enriched.categories.", 
                                     num.peptides[i], ".", dest.file, ".csv"), 
                              row.names = FALSE, quote = FALSE, sep = "\t")
                  
                  
                  if (!nrow(reducedTerms) < 50) {
                    
                    pdf(file=paste0(file.name.plots, ".reduced.terms.pdf"), width = 10, height = 6)
                    
                    treemap::treemap(reducedTerms[1:50,], title = "",
                                     index = c("parentTerm", "term"),
                                     fontsize.labels=c(12, 10),
                                     fontcolor.labels=c("white","black"),
                                     vSize = "score",
                                     fontface.labels=c(2,1),
                                     bg.labels=c("transparent"),
                                     type = "index",
                                     overlap.labels=1)
                    dev.off()
                    
                  } else {
                    
                    pdf(file=paste0(file.name.plots, ".reduced.terms.pdf"), width = 10, height = 6)
                    
                    treemap::treemap(reducedTerms, title = "",
                                     index = c("parentTerm", "term"),
                                     fontsize.labels=c(12, 10),
                                     fontcolor.labels=c("white","black"),
                                     vSize = "score",
                                     fontface.labels=c(2,1),
                                     bg.labels=c("transparent"),
                                     type = "index",
                                     overlap.labels=1)
                    dev.off()
                    
                  }
                  
                }
                
                
                
                
                # 7.5.2. Specific Gene Ontology (m) plots
                
                
                a <- ora_results_list[[i]][[1]][[j]][[k]][[l]][[m]]
                a.df <- data.frame(a)
                
                
                # in case a = 0, make it have 0 rows
                if (!nrow(data.frame(a)) == 0) {
                  if (data.frame(a)[1,1] == 0) {
                    a <- NULL
                    a.df <- data.frame()
                  }
                }
                
                # do not plot in case there are no enriched categories:
                if (!nrow(a.df) == 0) {
                  
                  edox <- setReadable(a, 'org.Hs.eg.db', 'ENTREZID')
                  
                  # barplot
                  ggplot(a, showCategory = showCategory.plots, 
                         aes(x = Count, y = Description, fill = p.adjust)) + 
                    ggtitle(title) +
                    geom_bar(stat = "identity", alpha = 1) + xlab("Number of proteins in each category") + 
                    ylab(NULL) + scale_fill_gradient(low = "red", high = "blue") +
                    ggtitle(title) +
                    theme_ridges() + theme(axis.title.x = element_text(hjust = 0.5),
                                           plot.title = element_text(hjust = 0.5))
                  ggsave(paste0(file.name.plots, ".barplot.png"), 
                         width = 10, height = 6)
                  
                  if ((length(a.df$Description) != 1) & (length(a@gene) != 1)) {
                    
                    # cnetplot
                    cnetplot(edox,
                             foldChange = ora_results_list[[i]][[3]][[j]][[m]],
                             circular = TRUE, colorEdge = TRUE, 
                             cex_label_gene = 0.6, cex_label_category = 0.7,
                             cex_category = 0.9,
                             showCategory = showCategory.plots) 
                    ggsave(paste0(file.name.plots, ".cnetplot.png"), 
                           width = 10, height = 8)
                    
                    #cnetplot2
                    cnetplot(edox, 
                             node_label="all",
                             foldChange = ora_results_list[[i]][[3]][[j]][[m]],
                             colorEdge = TRUE,
                             cex_label_category = 0.8,
                             cex_label_gene = 0.6,
                             showCategory = showCategory.plots) 
                    ggsave(paste0(file.name.plots, ".cnetplot2.png"), 
                           width = 10, height = 8)
                    
                  }
                  
                  # upsetplot: at least three categories to compare
                  if (nrow(a.df) > 3) {
                    upset <- upsetplot(a, showCategory = showCategory.plots) + 
                      ylab("Proteins") + 
                      xlab("Protein sets overlap")
                    ggsave(paste0(file.name.plots, ".upsetplot.png"), upset, 
                           width = 12, height = 8)
                  }
                  
                  # heatplot 
                  heatplot(edox, foldChange = ora_results_list[[i]][[3]][[j]][[m]], 
                           showCategory = showCategory.plots) + 
                    coord_equal()
                  ggsave(paste0(file.name.plots, ".heatplot.png"),
                         width = 10, height = 4)
                  
                  tryCatch({
                    # goplot
                    goplot(a, showCategory = showCategory.plots)
                    ggsave(paste0(file.name.plots, ".goplot.png"),
                           width = 13, height = 6)
                  }, error = function(e) print("goplot has not been able to be plotted"))
                  
                  
                } # end condition no enriched categories plot
                
              } # end specific GO ontology loop (m)
              
            } # end specific GO ontologies
            
            else {
              
              
              # 7.4. Ontology loop (l)
              
              
              # create regulation folder
              path.reg <- paste0(pathh, "/", dest.ontology)
              dest.reg <- names(ora_results_list[[i]][[3]][[j]][k])
              dir.create(file.path(path.reg, dest.reg), showWarnings = FALSE)
              
              # create two folders inside: enriched.categories and plots
              path.reg.2 <- paste0(path.reg, "/", dest.reg)
              dest.reg.2 <- "enriched.categories"
              dest.reg.3 <- "plots"
              dir.create(file.path(path.reg.2, dest.reg.2), showWarnings = FALSE)
              dir.create(file.path(path.reg.2, dest.reg.3), showWarnings = FALSE)
              
              # define a general title for plots and file names:
              file.name.categories <- paste0(path.reg.2, "/", dest.reg.2, "/",
                                             dest.ontology, ".ora.", dest.deps)
              file.name.plots <- paste0(path.reg.2, "/", dest.reg.3, "/",
                                        dest.ontology)
              title <- paste0(dest.ontology, " enriched categories")
              
              
              # 7.4. Ontology loop (l) categories
              
              
              enriched.categories <- data.frame(ora_results_list[[i]][[1]][[j]][[k]][[l]])
              colnames(enriched.categories)[which(names(enriched.categories) == "geneID")] <- "proteinID"
              enriched.categories <- core.enrichment.translation(enriched.categories,
                                                                 col.name = "proteinID")
              
              write.table(enriched.categories,
                          paste0(file.name.categories, ".enriched.categories.", 
                                 num.peptides[i], ".", dest.file, ".csv"), 
                          row.names=FALSE, quote=FALSE, sep = "\t")
              
              
              # 7.4. Ontology loop (l) plots
              
              
              a <- ora_results_list[[i]][[1]][[j]][[k]][[l]]
              a.df <- data.frame(a)
              
              
              # in case a = 0, make it have 0 rows
              if (!nrow(data.frame(a)) == 0) {
                if (data.frame(a)[1,1] == 0) {
                  a <- NULL
                  a.df <- data.frame()
                }
              }
              
              
              # do not plot in case there are no enriched categories:
              if (!nrow(a.df) == 0) {
                
                edox <- setReadable(a, 'org.Hs.eg.db', 'ENTREZID')
                
                # barplot
                ggplot(a, showCategory = showCategory.plots, 
                       aes(x = Count, y = Description, fill = p.adjust)) + 
                  ggtitle(title) +
                  geom_bar(stat = "identity", alpha = 1) + xlab("Number of proteins in each category") + 
                  ylab(NULL) + scale_fill_gradient(low = "red", high = "blue") +
                  ggtitle(title) +
                  theme_ridges() + theme(axis.title.x = element_text(hjust = 0.5),
                                         plot.title = element_text(hjust = 0.5))
                
                ggsave(paste0(file.name.plots, ".barplot.png"), width = 10, height = 6)
                
                
                if ((length(a.df$Description) != 1) & (length(a@gene) != 1)) {
                  
                  # cnetplot
                  cnetplot(edox,
                           foldChange = ora_results_list[[i]][[3]][[j]][[k]],
                           circular = TRUE, colorEdge = TRUE, 
                           cex_label_gene = 0.6, cex_label_category = 0.7,
                           cex_category = 0.9, 
                           showCategory = showCategory.plots) 
                  ggsave(paste0(file.name.plots, ".cnetplot.png"), 
                         width = 10, height = 8)
                  
                  #cnetplot2
                  cnetplot(edox, 
                           node_label="all",
                           foldChange = ora_results_list[[i]][[3]][[j]][[k]],
                           colorEdge = TRUE,
                           cex_label_category = 0.8,
                           cex_label_gene = 0.6,
                           showCategory = showCategory.plots) 
                  ggsave(paste0(file.name.plots, ".cnetplot2.png"), 
                         width = 10, height = 8)
                  
                }
                
                # upsetplot: at least three categories to compare
                if (nrow(a.df) > 3) {
                  upset <- upsetplot(a, showCategory = showCategory.plots) + 
                    ylab("Proteins") + 
                    xlab("Protein sets overlap")
                  ggsave(paste0(file.name.plots, ".upsetplot.png"), upset, 
                         width = 12, height = 8)
                }
                
                # heatplot:at least two categories
                if (nrow(a.df) > 2) {
                  heatplot(edox, foldChange = ora_results_list[[i]][[3]][[j]][[k]], 
                           showCategory = showCategory.plots) + 
                    coord_equal()
                  ggsave(paste0(file.name.plots, ".heatplot.png"), 
                         width = 10, height = 4)
                }
              } # end condition no enriched categories plot
              
            } # end no GO ontology condition
            
          } # end ontology loop (l)
          
        } # end regulation loop (k)
        
      } # end fold change loop (j)
      
      
      # 7.5. Peptides results
      
      
      # save ora enriched categories count for all peptides:
      count.ora.result <- rbind(count.ora.result, as.data.frame(ora_results_list[[i]][[2]]))
      
    } # end protein list loop (i)
    
    
    # 7.6. Summary results
    
    
    # create ora.enriched.categories folder in summary:
    dir.create(file.path("summary", "ora.enriched.categories"), showWarnings = FALSE)
    
    # create fold change folder:
    path.fc.reg <- paste0("summary/ora.enriched.categories")
    
    # fold change loop (i)
    for (i in 1:length(unique(count.ora.result$fold_change))) {
      fchange <- unique(count.ora.result$fold_change)[i]
      dest.sum.fc <- paste0("fold_change=", fchange)
      dir.create(file.path(path.fc.reg, dest.sum.fc), showWarnings = FALSE)
      
      # create regulation folder:
      path.sum.reg <- paste0(path.fc.reg, "/", dest.sum.fc)
      
      # regulation loop (j)
      for (j in 1:length(unique(count.ora.result$regulation))) {
        
        reg <- unique(count.ora.result$regulation)[j]
        dest.sum.reg <- paste0(reg, ".regulated")
        dir.create(file.path(path.sum.reg, dest.sum.reg), showWarnings = FALSE)
        file.pathh <- paste0(path.sum.reg, "/", dest.sum.reg)
        
        # save ora enriched categories count for each regulation:
        ora.res <- count.ora.result[count.ora.result$regulation==reg,]
        write.table(ora.res, paste0(file.pathh, "/", dest.sum.reg, ".ora.enriched.categories.csv"), 
                    row.names=FALSE, quote=FALSE, sep = "\t")
        
        ora.res$Peptides <- factor(ora.res$Peptides, levels = unique(ora.res$Peptides))
        
        # plot of enriched categories for each peptide and fold change:
        ggplot(ora.res, aes(x = factor(Peptides), y = Count, fill = Category)) +
          geom_bar(stat = 'identity', position="dodge", width = 0.5) +
          scale_y_continuous(breaks = seq(0, max(ora.res$Count), by = 1), 
                             minor_breaks = 1) + 
          xlab("Peptides") + ylab("Enriched categories") + 
          ggtitle(paste0(reg, "-regulated enriched categories"))
        ggsave(paste0(file.pathh, "/ora.enriched.categories.png"))
        
        # plot of differentially expressed proteins:
        ggplot(ora.res, aes(x = factor(Peptides), y = DEGs, fill = factor(fold_change))) +
          geom_bar(stat = 'identity', position="dodge", width = 0.5) +
          scale_y_continuous(breaks = seq(0, max(ora.res$DEGs), by = 1), 
                             minor_breaks = 1)+
          scale_fill_discrete(name = "fold change") + xlab("Peptides") +
          ylab("Differentially expressed proteins") +
          ggtitle(paste0(reg, "-regulated differentially expressed proteins"))
        ggsave(paste0(file.pathh, "/differentially.expressed.proteins.png"))
        
      } # end regulation loop (j) 
      
    } # end fold change loop (i)
    
    write.table(count.ora.result, paste0("summary/ora.enriched.categories/ora.enriched.categories.csv"), 
                row.names=FALSE, quote=FALSE, sep = "\t")
    
    
    
    # 7.7. Bubble chart
    
    # Save enriched categories for each ontology:
    
    # Ontology loop (i)
    
    if (ont == "all") {
      ont <- c("BP", "MF", "CC")
    }
    
    if (regulation == "all") {
      regulation <- c("UP", "DOWN", "UP.and.DOWN")
    }
    
    bubble.cols <- c("ID", "Description", "GeneRatio", "p.adjust", "Count")
    
    all.files <- list.files(getwd(), recursive = TRUE)
    
    if ("GO" %in% ora) {
      go.names <- paste0("GO\\(", ont, ")")
      ora.names <- c(go.names, ora[ora != "GO"])
    } else {
      ora.names <- ora
    }
    
    # ontology loop (i)
    
    for (i in 1:length(ora.names)) {
      
      
      # regulation loop (b)
      
      for (b in 1:length(regulation)) {
        
        if ("GO" %in% ora) {
          go.names <- paste0("GO\\(", ont, ")")
          ora.names <- c(go.names, ora[ora != "GO"])
        } else {
          ora.names <- ora
        }
        
        assign(paste0(regulation[b], ".files"), 
               grep(paste0(ora.names[i], ".ora.", regulation[b], ".regulated.enriched.categories.*peptide.*.csv"), 
                    all.files, value = TRUE))
        
        
        # Fold change loop (k): add a column containing this information in the 
        # merged enriched categories for each ontology
        
        for (k in 1:length(fc)) {
          
          assign(paste0(regulation[b], ".files.fc"), 
                 grep(paste0("fold_change=", fc[k]), get(paste0(regulation[b], ".files")), value = TRUE))
          
          
          assign(paste0(regulation[b], ".regulated.enriched.categories"), 
                 data.frame())
          
          # Peptide loop (j): we want to merge all categories for the different peptides
          
          
          for (j in 1:length(peptides)) {
            
            # extract selected files 
            # add a column indicating the peptide
            
            # "pept.DOWN.enriched.cat"
            assign(paste0("pept.", regulation[b], ".enriched.cat"),
                   read.csv(get(paste0(regulation[b], ".files.fc"))[j], header = TRUE, sep = "\t", dec = ".")[,bubble.cols])
            
            pept.enriched.cat <- get(paste0("pept.", regulation[b], ".enriched.cat"))
            
            
            if (nrow(pept.enriched.cat) == 0) {
              pept.enriched.cat[nrow(pept.enriched.cat)+1,] <- NA
            }
            
            pept.enriched.cat$peptides <- factor(peptides[j])
            
            assign(paste0(regulation[b], ".regulated.enriched.categories"),
                   rbind(get(paste0(regulation[b], ".regulated.enriched.categories")),
                         pept.enriched.cat))
            
          } # end peptide loop (j)
          
          
          # create bubble chart for each regulation and save in summary folder 
          enriched.cat <- paste0(regulation[b], ".regulated.enriched.categories")
          
          if (nrow(get(enriched.cat)) != 0) {
            
            # arrange in ascending order by p.adjust:
            enriched.matrix <- get(enriched.cat) %>% arrange(p.adjust) 
            
            # as some categories are enriched in different protein lists
            # we make sure we plot unique categories using showCategory.plots
            
            categories.subset.to.plot <- unique(enriched.matrix$Description)[1:showCategory.plots]
            
            if (paste0("GO\\(", ont[i], ")") %in% ora.names[i]) {
              
              ora.names[i] <- paste0(paste0("GO(", ont[i], ")"))
              
            }
            
            
            
            # do not plot in case of non existing enriched categories
            if (nrow(enriched.matrix[is.na(enriched.matrix[, 1]),]) != nrow(enriched.matrix)) {
              
              x_label_values = levels(enriched.matrix$peptides)
              
              enriched.matrix <- enriched.matrix[complete.cases(enriched.matrix), ]
              
              enriched.matrix[enriched.matrix$Description %in% categories.subset.to.plot,] %>%
                ggplot(aes(x = peptides, 
                           y = factor(Description, levels = rev(unique(Description))))) +
                geom_point(aes(size = Count, 
                               color = p.adjust), alpha = 0.7) +
                #Adjust the range of points size:
                scale_size(name = "Count", breaks = function(x) unique(floor(pretty(seq(1, (max(x) + 1) * 1.1)))),
                           range = c(2, 14)) +
                scale_color_gradient(low = "red", high = "blue") + 
                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=0)) +
                # drop makes sure all levels are plot despite there are no values for them
                scale_x_discrete(position = "top", drop = FALSE, 
                                 breaks = x_label_values,
                                 labels = x_label_values) +
                xlab(paste0("Detected proteins (fold_change=", fc[k], ")")) +
                ylab(paste(ora.names[i], regulation[b], "regulated enriched categories")) 
              
              dir.create(file.path("summary/ora.enriched.categories", "bubble.chart"), showWarnings = FALSE)
              
              ggsave(paste0("summary/ora.enriched.categories/bubble.chart/", ora.names[i], ".", 
                            regulation[b], ".regulated.fold_change=", fc[k] ,".png"))
              
            } # bubble plot function
            
            
          } # end condition
          
          
        } # end fold change loop (k)
        
        
      } # end bubble chart regulation loop (b)
      
      
    } # end ontology loop (i)
    
    
  } # end ORA function  
  
  
} # end main.function
