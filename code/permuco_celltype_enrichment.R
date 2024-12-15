#!/usr/bin/env Rscript

# Code to run cell type enrichment analyses  
# This example uses the human cell landscape expression specificity matrix (ctd) but the zeisel one in data/zeisel_2018.rds could also be used
# This example also uses the PAML branch results but could also be adapted to use the branch-site results as seen in the hashed out code


library(permuco)


# Read in ETE3 and HyPhy res 
branch_res <- read.csv("results/S1_PAML_branch_res_Multiz30wayPrimate.csv")

# Read in CDS and GC content 
cds_gc <- read.csv("data/hg38_cds_len_gc_content.csv")

# Load zeisel ctd
ctd <- readRDS(file = "data/hcl2020_ctd.rds")

permuco_evo_res <- list()

for (i in seq_along(ctd)) {
  current_level <- ctd[[i]]
  
  # Iterate over each branch in evo_res
  for (branch_name in unique(evo_res$branch)) {
      
      # Filter hits based on omega and/or p value 
      hits <- evo_res$gene[evo_res$p < 0.05 & evo_res$branch == branch_name]
     #hits <- evo_res$gene[evo_res$p < 0.05 & evo_res$branch == branch_name] 
      
      in_genelist <- rep(0, nrow(current_level$specificity))
      in_genelist[rownames(current_level$specificity) %in% hits] <- 1
      
      permuco_results <- data.frame(celltype = colnames(current_level$specificity), p = 0, coefs = 0)
      
      for (ct in colnames(current_level$specificity)) {
        newDat <- data.frame(in_genelist = as.factor(in_genelist), spc = current_level$specificity[, ct],gene=rownames(current_level$specificity))

        # ... add the transcript_length data alongside the geneset binary representation
        newDat2=merge(newDat,cds_gc,by="gene",all.x=FALSE,all.y=FALSE)
        newDat2$percentage_gene_gc_content <- scale(newDat2$percentage_gene_gc_content, scale=FALSE)
        newDat2$cds_length <- scale(newDat2$cds_length, scale=FALSE)
        
        # The below line is the one that actually uses permuco
        mod <- lmperm(formula = spc ~ in_genelist * cds_length * percentage_gene_gc_content, data = newDat2, np = 10000)
        
        permuco_results[permuco_results$celltype == ct, ]$p <- mod$table$`resampled Pr(>|t|)`[2]
        permuco_results[permuco_results$celltype == ct, ]$coefs <- mod$coefficients[2]
      }
      
      permuco_evo_res[[branch_name]]$permuco_results <- permuco_results
      permuco_evo_res[[branch_name]]$permuco_results$ctd_level <- i
      permuco_evo_res[[branch_name]]$permuco_results$branch <- branch_name

    }
}
    


