# code to split MAF (.fa) according to gene symbol and rename so that you only have ensembl transcript ID and species ID e.g. ENST00000000233.10_hg38
#  instead of e.g. >ENST00000367772.8_hg38 2229 chr1:169853713-169888840-

# Load packages
library(Biostrings) 
library(dplyr)
library(tidyr)
library(parallel)
library(data.table)
library(rtracklayer)
library(biomaRt)

# Download fasta file from: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz30way/alignments/knownCanonical.protNuc.fa
# fasta="/rds/general/user/cm1118/projects/cellevolution/live/EvoCT_coding/data/knownCanonical.protNuc.fa"

# Load FASTA file
genes_fasta <- readDNAStringSet(fasta)

# Convert to data.table for speed
fasta_dt <- data.table(
  seq.name = names(genes_fasta),
  sequence = as.character(genes_fasta)
)

# Remove duplicates
fasta_dt <- unique(fasta_dt)

# Split seq.name into components: 
fasta_dt[, c("transcript_id", "null", "chr") := tstrsplit(seq.name, " ", fixed = TRUE)]
fasta_dt$null <- NULL
fasta_dt[, c("transcript_id", "species_id") := tstrsplit(transcript_id, "_", fixed = TRUE)]
fasta_dt$seq.name <- NULL

# Filter sequences mapping to alt contigs or fix patches 
fasta_dt <- fasta_dt[!grepl("_alt|fix", chr)]

# Create a UCSC table query for hg38 knownCanonical
session <- browserSession("UCSC")
genome(session) <- "hg38"

# Define the query 
query <- ucscTableQuery(session, "GENCODE V36")
data <- getTable(query)

# Subset to transcript ID and gene symbol
gene_symbols <- data[,c("name", "geneName")]
gene_symbols <- unique(gene_symbols)
setnames(gene_symbols, c("transcript_id", "gene_symbol"))

# Merge FASTA and gene-symbol table
merged_dt <- merge(fasta_dt, gene_symbols, by = "transcript_id", all.x = TRUE)
merged_dt <- merged_dt[!is.na(gene_symbol)]

# Find which genes have more than 1 transcript and subset to canonical transcript using ensembl data 
transcript_counts <- merged_dt[, .N, by = gene_symbol]
multiple_transcripts <- transcript_counts[N > 30]

mart <- useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
canonical_tx <- getBM(attributes=c("ensembl_gene_id",
                                   "ensembl_transcript_id_version",
                                   "hgnc_symbol",
                                   "transcript_is_canonical"),
                      mart=mart, 
                      filters = "hgnc_symbol", 
                      values = multiple_transcripts$gene_symbol)

canonical_tx_subset <- subset(canonical_tx, is.na(transcript_is_canonical))

merged_dt <- merged_dt[!merged_dt$transcript_id %in% canonical_tx_subset$ensembl_transcript_id_version,]

# The remaining genes with 2 transcripts are the same transcript annotated to X and Y chromosomes so we keep 1 
# apart from MATR3 which I manually filtered 
merged_dt <- unique(merged_dt)
merged_dt <- subset(merged_dt, transcript_id != "ENST00000361059.7")

# Create seq.name column to include transcript and species IDs
#merged_dt <- unite(merged_dt, seq.name, c(transcript_id, species_id), sep="_")


# Group sequences by gene_symbol and save
save_split_fasta <- function(gene_data, save_path) {
  # Get the unique gene name
  gene <- unique(gene_data$gene_symbol)
  
  # Create a subdirectory for the gene
  gene_dir <- file.path(save_path, gene)
  if (!dir.exists(gene_dir)) {
    dir.create(gene_dir, recursive = TRUE)
  }
  
  # Prepare the DNAStringSet
  seqs <- Biostrings::DNAStringSet(setNames(gene_data$sequence, gene_data$species_id))
  
  # Save the sequences to a FASTA file in the gene-specific directory
  Biostrings::writeXStringSet(seqs, filepath = file.path(gene_dir, paste0(gene, ".fasta")))
}


#save_path = "/rds/general/user/cm1118/projects/cellevolution/live/EvoCT_coding/data"

# Parallelized saving
no_cores <- max(detectCores() - 1, 1)
cl <- makeCluster(no_cores)
clusterExport(cl, varlist = c("save_split_fasta", "merged_dt"))
parLapply(cl, split(merged_dt, merged_dt$gene_symbol), save_split_fasta, save_path = save_path)
stopCluster(cl)
