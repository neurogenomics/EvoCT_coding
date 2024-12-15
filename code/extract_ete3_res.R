#!/usr/bin/env Rscript

# Code for extracting results from ete3 analysis

library(dplyr)
library(stringr)

# read in tree marks (found in data folder of the repo) 
#tree_marks <- read.csv("/rds/general/user/cm1118/home/EvolutionOfCelltypes/data/tree_marks.csv")

#######################################
# 1. Extract data from taxa.txt files
#######################################

taxaFiles = list.files(pattern=".taxa.txt", recursive=TRUE, full.names=TRUE)

count=0
for(tf in taxaFiles){
  taxa = readLines(tf)
  if(length(taxa)>0){
    geneName = gsub(".*/", "", tf)
    geneName = gsub(".taxa.txt", "", geneName)
    tmpTaxa = data.frame(gene=geneName)
    taxaList = as.list(taxa)
    for(i in 1:length(taxaList)){
      species = taxaList[[i]]
      tmpTaxa[[species]] = 1
    }
    
    ###############################################################################################################
    # 2. Generate a table to record how many taxa (and which ones) were included in the analysis of each gene
    ###############################################################################################################
    
    count=count+1
    if(count==1){
      taxaAnalysed=tmpTaxa
    }else{
      taxaAnalysed=bind_rows(taxaAnalysed, tmpTaxa)
    }
  }
}

taxaAnalysed[is.na(taxaAnalysed)] <- 0
taxaAnalysed <- taxaAnalysed %>% mutate(total=rowSums(select_if(., is.numeric)))

######################################
# 3. Get LRT p-values for all models
######################################

logFiles <- list.files(pattern=".log", recursive=TRUE, full.names=TRUE)

for (lf in logFiles){
  if(length(grep("\\.branch.log", lf) > 0)){
    geneName_b = gsub(".*/", "", lf)
    geneName_b = gsub(".branch.log", "", geneName_b)

    if (geneName_b %in% allBranchRes$gene) {
      next  # Skip only the if statement and continue with the next else if statement
    }

    dat_b = readLines(lf)
    whichStart_b = grep("LRT",dat_b) # Block containing the p values for the LRTs starts with this term
    if(length(whichStart_b)>0){
      whichEnd_b = grep("SUMMARY BY MODEL",dat_b)
      tableLines_b = dat_b[(whichStart_b+4):(whichEnd_b-2)]
      res_b = strsplit(tableLines_b," \\| ")
      for(i in 1:length(res_b)){
        null_b = trimws(res_b[[i]][1])
        alt_b = trimws(res_b[[i]][2])
        p_b = as.numeric(gsub("\\*","",trimws(res_b[[i]][3])))
        null_species_b = gsub("\\~.*|.*\\.|\\-.*", "\\1", alt_b)
        alt_species_b = gsub("\\~.*|.*\\.|\\-.*", "\\1", alt_b)
        null_model_b = gsub("(.*)\\..*|\\~.*","\\1",null_b)
        alt_model_b = gsub("(.*)\\..*|\\~.*","\\1",alt_b)
        tmp_branch = data.frame(gene=geneName_b,null_model=null_model_b,alt_model=alt_model_b,species=alt_species_b,p=p_b)
        if(i==1){
          BranchRes=tmp_branch
        }else{
          BranchRes=rbind(BranchRes,tmp_branch)
        }
        
      }
      
      #########################################
      # 4. Get omega values for branch model
      #########################################
      
      omega_start_lines = grep("- Model b_free.",dat_b) # Each block containing an omego value starts with this term
      BranchRes$omega = NA
      for(jj in omega_start_lines){
        subBlock_omega = dat_b[jj:(jj+7)]
        whichLine = grep("#1  =>",subBlock_omega)
        subBlock_omega[whichLine]
        omega = as.numeric(gsub("(.*#1\\s*=>\\s*)(\\d*)","\\2",subBlock_omega[whichLine]))
        whichModel = gsub("(.*Model[ ]*)(b_free.*\\..*)-(.*)","\\2",subBlock_omega[1])
        whichModel = gsub("b_free.", "", whichModel)
        BranchRes[BranchRes$null_model=="b_neut" & BranchRes$alt_model=="b_free" & grepl(whichModel,BranchRes$species),]$omega =  omega
        
      }
      
      ###############################################
      # 5. Get tree marks for branch model    
      ###############################################
      
      tree_mark_start_lines = grep("marking branches",dat_b)
      BranchRes$tree = NA
      for(kk in tree_mark_start_lines){
        subBlock_tree = dat_b[kk+2]
        tree_mark = gsub("#1.*","#1",subBlock_tree)
        tree_mark_2 = gsub(".*,|\\(","",tree_mark)
        branch_number = gsub(".*?([0-9]+).*","\\1", dat_b[kk])
        BranchRes[BranchRes$species==branch_number,]$tree = tree_mark_2
        BranchRes$tree = trimws(BranchRes$tree, which = c("both"))
      }
      
      BranchRes <- merge(BranchRes, tree_marks, by="tree")
      BranchRes$tree <- NULL
      
      ################################################
      # 6. Construct results table for branch model
      ################################################
      
      if(!exists('allBranchRes')){
        allBranchRes=BranchRes
      }else{
        allBranchRes=rbind(allBranchRes,BranchRes)
      }
    }
    
  } else if(length(grep("\\.site.log", lf) > 0)){
    dat_s = readLines(lf)
    whichStart_s = grep("LRT",dat_s) # Block containing the p values for the LRTs starts with this term
    if(length(whichStart_s)>0){
      whichEnd_s = grep("SUMMARY BY MODEL",dat_s)
      tableLines_s = dat_s[(whichStart_s+4):(whichEnd_s-2)]
      res_s = strsplit(tableLines_s," \\| ")
      for(j in 1:length(res_s)){
        geneName_s = gsub(".*/", "", lf)
        geneName_s = gsub(".site.log", "", geneName_s)
        null_s = trimws(res_s[[j]][1])
        alt_s = trimws(res_s[[j]][2])
        p_s = as.numeric(gsub("\\*","",trimws(res_s[[j]][3])))
        null_model_s = gsub("(.*)\\..*|\\~.*","\\1",null_s)
        alt_model_s = gsub("(.*)\\..*|\\~.*","\\1",alt_s)
        tmp_site = data.frame(gene=geneName_s,null_model=null_model_s,alt_model=alt_model_s,p=p_s)

        
        ##############################################
        # 7. Construct results table for site model
        ##############################################
        
        if(!exists('allSiteRes')){
          allSiteRes=tmp_site
        }else{
          allSiteRes=rbind(allSiteRes,tmp_site)
        }
      }
    }
    
  } else{

    geneName_bs = gsub(".*/", "", lf)
    geneName_bs = gsub(".branch-site.log", "", geneName_bs)

    if (geneName_bs %in% allBranchSiteRes$gene) {
      next  # Skip only the if statement and continue with the next else if statement
    }

    dat_bs = readLines(lf)
    whichStart_bs = grep("LRT", dat_bs) # Block containing the p values for the LRTs starts with this term
    if(length(whichStart_bs)>0){
      whichEnd_bs = grep("SUMMARY BY MODEL",dat_bs)
      tableLines_bs = dat_bs[(whichStart_bs+4):(whichEnd_bs-2)]
      res_bs = strsplit(tableLines_bs," \\| ")
      for(k in 1:length(res_bs)){
        null_bs = trimws(res_bs[[k]][1])
        alt_bs = trimws(res_bs[[k]][2])
        p_bs = as.numeric(gsub("\\*","",trimws(res_bs[[k]][3])))
        null_species_bs = gsub("\\~.*|.*\\.|\\-.*", "\\1", null_bs)
        null_species_bs = gsub("^([^_]*_[^_]*)_.*$", "\\1", null_species_bs)
        alt_species_bs = gsub("\\~.*|.*\\.|\\-.*", "\\1", alt_bs)
        null_model_bs = gsub("(.*)\\..*|\\~.*","\\1",null_bs)
        alt_model_bs = gsub("(.*)\\..*|\\~.*","\\1",alt_bs)
        tmp_branchsite = data.frame(gene=geneName_bs,null_model=null_model_bs,alt_model=alt_model_bs,species=alt_species_bs,p=p_bs)
        if(k==1){
          BranchSiteRes=tmp_branchsite
        }else{
          BranchSiteRes=rbind(BranchSiteRes,tmp_branchsite)
        }
      }
      
      tree_mark_start_lines = grep("marking branches",dat_bs)
      BranchSiteRes$tree = NA
      for(kk in tree_mark_start_lines){
        subBlock_tree = dat_bs[kk+2]
        tree_mark = gsub("#1.*","#1",subBlock_tree)
        tree_mark_2 = gsub(".*,|\\(","",tree_mark)
        branch_number = gsub(".*?([0-9]+).*","\\1", dat_bs[kk])
        BranchSiteRes[BranchSiteRes$species==branch_number,]$tree = tree_mark_2
        BranchSiteRes$tree = trimws(BranchSiteRes$tree, which = c("both"))
      }
      
      BranchSiteRes <- merge(BranchSiteRes, tree_marks, by="tree")
      BranchSiteRes$tree <- NULL
      
      ####################################################
      # 8. Construct results table for branch-site model
      ####################################################
      
      if(!exists('allBranchSiteRes')){
        allBranchSiteRes=BranchSiteRes
      }else{
        allBranchSiteRes=rbind(allBranchSiteRes, BranchSiteRes)
      }
      
    }
  }
}

save.image("/rds/general/user/cm1118/home/EvolutionOfCelltypes/results_Jan22/allRes_ete3_nostopcodon_nobranchlen2.RData")

############################################################################################################
# 9. Get codon positions for positively selected, relaxed, and conserved sites (branch-site and site model)
############################################################################################################

for (lf in logFiles){
  if(length(grep("\\.branch-site.log", lf) > 0)){
    dat_bs_sites = readLines(lf)
    dat_bs_sites = append(dat_bs_sites, "") # append empty line to end of each log file
    sites_start_lines_bs = grep("* Sites significantly caracterized",dat_bs_sites)+1
    sites_end_lines_bs =  which(!nzchar(dat_bs_sites))
    for(ss in sites_start_lines_bs){
      lastLine_bs = sites_end_lines_bs[sites_end_lines_bs>ss][1]-1
      sitesTable_bs = dat_bs_sites[(ss-4):lastLine_bs]
      codonPositions_bs = dat_bs_sites[(ss+2):lastLine_bs]
      geneName_bs = gsub(".*/", "", lf)
      geneName_bs = gsub('\\..*', "", geneName_bs)
      model_species_bs = gsub(" - Model ", "", sitesTable_bs[1])
      model_bs = gsub('\\..*', "", model_species_bs)
      species_bs = gsub('\\-.*', "", model_species_bs)
      species_bs = gsub("^.*\\.","", species_bs)
      res_bs_sites = strsplit(codonPositions_bs," \\| ")
      for(l in 1:length(res_bs_sites)){
        position_bs=trimws(res_bs_sites[[l]][1])
        category_bs=trimws(res_bs_sites[[l]][2])
        tmp_bs_sites = data.frame(gene=geneName_bs,model=model_bs, species=species_bs, codon_position=position_bs, category_probability=category_bs)
        tmp_bs_sites$category = sapply(strsplit(as.character(tmp_bs_sites$category),' '), "[", 1)
        tmp_bs_sites$posterior_probability = sapply(strsplit(as.character(tmp_bs_sites$category_probability),'>'), "[", 2)
        tmp_bs_sites$posterior_probability = gsub(")", "",tmp_bs_sites$posterior_probability)
        tmp_bs_sites$category_probability = NULL
        tmp_bs_sites$model <- sub("^$", "M1", tmp_bs_sites$model)
        tmp_bs_sites$species <- sub("^$", "NA", tmp_bs_sites$species)
        
        if(!exists('allBranchSiteRes_codon_pos')){
          allBranchSiteRes_codon_pos=tmp_bs_sites
        }else{
          allBranchSiteRes_codon_pos=rbind(allBranchSiteRes_codon_pos, tmp_bs_sites)
        }
      }
      
    }
    
  } else if(length(grep("\\.site.log", lf) > 0)){
    dat_s_sites = readLines(lf)
    dat_s_sites = append(dat_s_sites, "") # append empty line to end of each log file
    sites_start_lines_s = grep("* Sites significantly caracterized",dat_s_sites)+1
    sites_end_lines_s =  which(!nzchar(dat_s_sites))
    for(tt in sites_start_lines_s){
      lastLine_s = sites_end_lines_s[sites_end_lines_s>tt][1]-1
      sitesTable_s = dat_s_sites[(tt-3):lastLine_s]
      codonPositions_s = dat_s_sites[(tt+2):lastLine_s]
      geneName_s = gsub(".*/", "", lf)
      geneName_s = gsub('\\..*', "", geneName_s)
      model_species_s = gsub(" - Model ", "", sitesTable_s[1])
      model_s = gsub('\\..*', "", model_species_s)
      species_s = gsub('\\-.*', "", model_species_s)
      species_s = gsub("^.*\\.","", species_s)
      res_s_sites = strsplit(codonPositions_s," \\| ")
      for(m in 1:length(res_s_sites)){
        position_s=trimws(res_s_sites[[m]][1])
        category_s=trimws(res_s_sites[[m]][2])
        tmp_s_sites = data.frame(gene=geneName_s,model=model_s, model=model_bs, codon_position=position_s, category_probability=category_s)
        tmp_s_sites$category = sapply(strsplit(as.character(tmp_s_sites$category),' '), "[", 1)
        tmp_s_sites$posterior_probability = sapply(strsplit(as.character(tmp_s_sites$category_probability),'>'), "[", 2)
        tmp_s_sites$posterior_probability = gsub(")", "",tmp_s_sites$posterior_probability)
        tmp_s_sites$category_probability = NULL
        tmp_s_sites$model <- sub("^$", "M1", tmp_s_sites$model)
        tmp_s_sites$species <- sub("^$", "NA", tmp_s_sites$species)
        if(!exists('allSiteRes_codon_pos')){
          allSiteRes_codon_pos=tmp_s_sites
        }else{
          allSiteRes_codon_pos=rbind(allSiteRes_codon_pos,tmp_s_sites)
        }
      }
    }
  }
}

