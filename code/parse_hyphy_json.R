#!/usr/bin/env Rscript

rm(list = ls())

library(jsonlite)
library(purrr)
library(dplyr)

# code to parse hyphy absrel json results files

setwd("/rds/general/user/cm1118/projects/cellevolution/live/genes")

files <- list.files(full.names = TRUE, pattern = '*_ABSREL.json', recursive = TRUE)

# function to extract species from the left side of the node
get_string_left_of_node <- function(node_name) {
  if (startsWith(node_name, "Node")) {
    # Find the position of the node_name in the tree_string
    node_position <- regexpr(paste0("\\b", node_name, "\\b"), tree_string)
    
    # If node_name is found in the tree_string
    if (node_position != -1) {
      # Extract the substring before the node_name
      left_string <- trimws(substr(tree_string, 1, node_position - 1))
      return(left_string)
    }
  }
  # If node_name is not found, return NA
  return(node_name)
}

results <- list()

for (file in files) {
  tryCatch({
    json <- read_json(file)

    # Get gene name from file name
    base_name <- basename(file)
    gene <- gsub("_ABSREL.json", "", base_name)
    
    # extract the contents of the 'trees' key
    tree_string <- json$input$trees$`0`

    # Initialise an empty list for the file results
    file_results <- list()
    
    # get node names that were tested 
    nodes_tested <- names(json$tested$`0`[json$tested$`0` == "test"])

    # Iterate over the nodes
    for (node_name in nodes_tested) {
      
      species_left <- get_string_left_of_node(node_name)
      
      # Initialize branch as "Homo-Pan"
      branch <- "Homo-Pan"
      
      # Classify branches based on the extracted species
      if (node_name == "hg38") {
        branch <- "Human"
      } else if (grepl("otoGar3", species_left)) {
        branch <- "Primate"
      } else if (grepl("tarSyr2", species_left)) {
        branch <- "Haplorhine"
      } else if (grepl("ponAbe2", species_left)) {
        branch <- "Great Ape"
      }

      # Check if the node exists in the json
      if (!is.null(json[["branch attributes"]][["0"]][[node_name]])) {
        file_results[[node_name]] <- data.frame(gene = gene, node_id = node_name, branch = branch, p = json[["branch attributes"]][["0"]][[node_name]][["Uncorrected P-value"]])
      } else {
        file_results[[node_name]] <- data.frame(gene = gene, node_id = node_name, branch = branch, p = "NULL")
      }
    }

    # Append file results to the overall results list
    results[[file]] <- file_results

  }, error = function(e) {
    # Error occurred during parsing, continue to the next file
    return()
  })
}

results_list <- lapply(results, function(element) {
  do.call(rbind, element)
})

# Unlist the combined list into a single data frame
results_df <- do.call(rbind, results_list)


write.csv(results_df, file="../hyphy_absrel_res.csv", quote=F, row.names=F)
save.image("../hyphy_absrel_res.RData")