#!/usr/bin/Rscript
## Author: Taylor Falk
## tfalk@bu.edu
## BU BF591
## Assignment Week 2

#### Bioconductor ####
# it is standard among R packages to define libraries and packages at the 
# beginning of a script. Also note that a package should NOT be installed every 
# time a script runs.
# The bioconductor repository has installation instructions for biomaRt: 
# https://bioconductor.org/install/

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

if (!require("BiocManager", quietly = TRUE)){
  
}
if (!require("biomaRt", quietly = TRUE)){
  
}
# load tidyverse and your new bioconductor package
library(tidyverse)
library(BiocManager)
#setwd('/project/bf528/project_1/data/')

#### Loading and processing data ####
#' Load Expression Data
#'
#' @param filepath A text string of the full filepath to the file to load.
#'
#' @return A tibble containing the data loaded from the CSV in `filepath`. 
#' 
#' @details Note that not all CSVs are created equal, and there are often cases where 
#' the data will not load in correctly on the first try. You may want to write this functon to 
#' adjust the CSV file being loaded into this assignment so that it can be formed into a 
#' tibble correctly.
#'
#' @examples 
#' `data <- load_expression('/project/bf528/project_1/data/example_intensity_data.csv')`
load_expression <- function(filepath) {
  data <- readr::read_delim(filepath," ")
  tbl <- tibble(data) 
  return(tbl)
}

#' Filter 15% of the gene expression values.
#'
#' @param tibble A tibble of expression values, rows by probe and columns by sample.
#'
#' @return A tibble of affymetrix probe names from the input expression data tibble. 
#' These names match the rows with 15% or more of the expression values about log2(15).
#' 
#' @details This is similar to the filters being implemented in BF528's project 1. 
#' We do not necessarily want to capture all parts of the assay in our analysis, so 
#' filters like this serve to reduce the noise and amount of data to examine.
#'
#' @examples `samples <- filter_15(data_tib)`
#' `> str(samples)`
#' `tibble [40,158 Ã 1] (S3: tbl_df/tbl/data.frame)`
#' `$ probeids: chr [1:40158] "1007_s_at" "1053_at" "117_at" "121_at" ...`
filter_15 <- function(data){
  thr <- .15 * (ncol(data[-1]))
  filtered <- data[rowSums(data[-1] > log2(15)) > (thr),]
  probeids_vect <- filtered[1]
  probeids_vect <- as.data.frame(probeids_vect)
  return(probeids_vect)
}
#expressed_filter<-data[-1][rowSums(data[-1] > log2(15)) >= (0.15*ncol(data)),]

#' `> affy_to_hgnc(tibble(c('202860_at', '1553551_s_at')))`
#' `affy_hg_u133_plus_2 hgnc_symbol`
#' `1        1553551_s_at      MT-ND1`
#' `2        1553551_s_at       MT-TI`
#' `3        1553551_s_at       MT-TM`
#' `4        1553551_s_at      MT-ND2`
#' `5           202860_at     DENND4B`

affy_to_hgnc <- function(affy_vector){
  #library(bioMart)
  BiocManager::install("affy")
  ensembl = useMart("ensembl", dataset ="hsapiens_gene_ensembl")
  gene_map <- as_tibble(
    getBM(attributes=c("affy_hg_u133_plus_2", "hgnc_symbol"), 
          filters = "affy_hg_u133_plus_2",
          values= affy_vector,
          mart=ensembl)
  )
  return(gene_map)
}

#### ggplot ####

#' Reduce a tibble of expression data to only the rows in good_genes or bad_genes.
#'
#' @param expr_tibble A tibble holding the expression data, each row corresponding to
#' one affymetrix probe ID and each column to a sample.
#' @param names_ids A two column tibble that associates affy IDs with HGNC gene IDs. 
#' Generated `with affy_to_hgnc()`.
#' @param good_genes A list of gene names stored as a vector of strings.
#' @param bad_genes A list of gene names stored as a vector of strings.
#'
#' @return A tibble with two additional columns added:
#' 1. HGNC gene IDs 
#' 2. Does the gene is this row fall into "good" or "bad" genes?
#' This tibble should be reduced to only rows present in good or bad genes. All
#' other rows can be discarded.
#' 
#' @details In order to plot only our genes of interest, we need to rearrange our 
#' data to include only the elements we want to see. We also want to add to columns, 
#' one that associates the probeids with the HGNC gene name, and one that says if 
#' that gene is in the good or bad sets of genes.
#'
#' @examples 
#' `plot_tibble <- reduce_data(expr_tibble = expr, names_ids = sample_names,`
#' `                           goodGenes, badGenes)`
#' `> head(plot_tibble)`
#' `A tibble: 6 Ã 38`
#' `  probeids    hgnc    gene_set    GSM972389 ...`
#' `  <chr>       <chr>   <chr>       <dbl>     ...`
#' `1 202860_at   DENND4B good        7.16      ...`
#' `2 204340_at   TMEM187 good        6.40      ...`
reduce_data <- function(expr_tibble, names_ids, good_genes, bad_genes){
  reduce_sample_names <- as_tibble(filter(names_ids,
                                          (hgnc_symbol %in% good_genes) | (hgnc_symbol %in% bad_genes)))
  
  reduced_expr <- as_tibble(filter(expr_tibble, probeids %in% reduce_sample_names$affy_hg_u133_plus_2))
  combined_expr <- as_tibble(merge(x=reduced_expr,y=reduce_sample_names[,c("affy_hg_u133_plus_2","hgnc_symbol")],
                                   by.x = "probeids", by.y = "affy_hg_u133_plus_2"))
  combined_expr <- combined_expr %>% relocate(hgnc_symbol, .after = probeids)
  combined_expr <- combined_expr %>% add_column(gene_set = "good", .after = "hgnc_symbol")
  combined_expr <- combined_expr %>% mutate_at(vars(c('gene_set')),~ifelse(hgnc_symbol %in% bad_genes, "bad", .))
  return(combined_expr)
}

#' Plot a boxplot of good and bad genes.
#'
#' @param tibble A reduced tibble of expression data, with information about
#' good and bad genes and gene names.
#'
#' @return A ggplot object which contains a boxplot of the genes and samples we 
#' are interested in.
#' 
#' @details This function performs one additional step before using `ggplot()`: 
#' converting the _wide_ format of the input tibble to a _long_ format.
#'
#' @examples `p <- plot_ggplot(plot_tibble)`
plot_ggplot <- function(data_tibble) {
  
  goodGenes <- c("PKD1", "NOS3", "AGTR1", "COL4A5", "ECE1", "MEN1", "OLR1", "F7")
  badGenes <- c("TP53", "EGFR", "BRAF", "KRAS", "PIK3CA", "ERBB2", "MAPK1", "NRAS")
  
  long_tibble <- data_tibble%>%
    pivot_longer(cols = starts_with("GSM"), names_to = "GSM", values_to = "Expression_levels")
  
  # group long_tibble by good and bad genes
  long_tibble$hgnc_symbol <- factor(long_tibble$hgnc_symbol, levels=c(goodGenes, badGenes))
  #view(grouped_data)
  
  # Default plot
  box_plot <- ggplot(long_tibble, aes(x = hgnc_symbol, y = Expression_levels, 
                                      fill = gene_set)) + geom_boxplot() + theme_bw() +
    scale_fill_manual(values=c("#DF2935", "#71B48D")) +
    labs(title = "Boxplot of 16 somewhat randomly chosen genes",
         subtitle = "like really I just kinda picked some random cancer ones")
  # Set the theme for box plot
  bplot <- box_plot + xlab("Gene") + ylab("Expression levels") 
  bplot <- bplot + theme(legend.position = "bottom", 
                         axis.text.x = element_text(angle = 45, vjust = 0.20, hjust = 0.40),
                         plot.title = element_text(size = 14),
                         plot.subtitle = element_text(size = 11))
  return(bplot)
}


