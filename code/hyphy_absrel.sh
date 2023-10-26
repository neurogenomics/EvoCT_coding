#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=64gb
#PBS -J 1-10000
#PBS -j oe
#PBS -o hyphy.log

# activate conda environment
module load anaconda3/personal
source activate hyphy

# Define environment variables
export repo_dir=/rds/general/user/cm1118/projects/cellevolution/live/EvoCT_coding
export analysis_dir=/rds/general/user/cm1118/projects/cellevolution/live/genes
export temp_dir=/rds/general/user/cm1118/projects/cellevolution/ephemeral


# Determine the gene name to be processed based on the PBS_ARRAY_INDEX (this works for Imperial College London's HPC array jobs)
gene=$(head -"$PBS_ARRAY_INDEX" $repo_dir/data/hyphy_genes.txt | tail -1)

# Define the path to the JSON result file
json_file="$analysis_dir/$gene/${gene}_ABSREL.json"

# check if JSON results file already exists
if [ -e "$json_file" ]; then
    echo "JSON file $json_file exists. Skipping analysis for $gene."
else
    echo "Running analysis for $gene"

    # Remove any existing log file for this gene.
    rm "$analysis_dir/$gene/${gene}_ABSREL.log"
    
    # count number of sequences in fasta 
    count=$(grep -o ">" "$analysis_dir/$gene/${gene}_nostop.fasta" | wc -l)

# hyphy requires full path names 
# If there are 30 sequences in the fasta, run HyPhy using the full species tree 
# If there aren't 30 sequences in the fasta then you need to re-label the species tree using label-mrca.bf (this function was written by HyPhy developers when I asked about labelling the most recent common ancestor of two species)
if [ "$count" -eq 30 ]; then
    hyphy absrel --alignment "/rds/general/user/cm1118/projects/cellevolution/live/genes/$gene/${gene}_nostop.fasta" --tree "/rds/general/user/cm1118/home/EvolutionOfCelltypes/data/hg38.30way.nh.hyphy.nobranchlen.txt" --branches Foreground --multiple-hits Double+Triple --srv Yes --output "/rds/general/user/cm1118/projects/cellevolution/live/genes/$gene/${gene}_ABSREL.json" >> "/rds/general/user/cm1118/projects/cellevolution/live/genes/$gene/${gene}_ABSREL.log"
else
    hyphy /rds/general/user/cm1118/projects/cellevolution/live/genes-analyses/LabelTrees/label-tree.bf --tree /rds/general/user/cm1118/projects/cellevolution/live/genes/"$gene"/species_tree_nobranchlen.nw --list /rds/general/user/cm1118/projects/cellevolution/ephemeral/hyphy_temp/list.txt --output "/rds/general/user/cm1118/projects/cellevolution/live/genes/$gene/hyphy_hg38.nwk"
    
    # Generate a list of taxa from the fasta file
    awk 'NR==FNR{A[$1];next}$2 && $3 in A' /rds/general/user/cm1118/projects/cellevolution/live/genes/"$gene"/"$gene.taxa.txt" /rds/general/user/cm1118/home/EvolutionOfCelltypes/data/branch_mark_dictionary_hyphy.txt | awk '!seen[$1]++ {print $4}' | awk 'BEGIN { ORS = ";" } { print substr($0, 1, length($0)-1) }' | sed 's/;$//' > "/rds/general/user/cm1118/projects/cellevolution/live/genes/$gene/$gene_branch_marks.txt"
    # Read the list of taxa
    taxa=$(cat "/rds/general/user/cm1118/projects/cellevolution/live/genes/$gene/$gene_branch_marks.txt")

    # Label most recent common ancestor
    hyphy /rds/general/user/cm1118/projects/cellevolution/live/genes-analyses/LabelTrees/label-mrca.bf --tree "/rds/general/user/cm1118/projects/cellevolution/live/genes/$gene/hyphy_hg38.nwk" --taxa "${taxa[*]}" --output "/rds/general/user/cm1118/projects/cellevolution/live/genes/$gene/hyphy_labelled.nwk"

    # Run HyPhy aBSREL
    hyphy absrel --alignment "/rds/general/user/cm1118/projects/cellevolution/live/genes/$gene/${gene}_nostop.fasta" --tree "/rds/general/user/cm1118/projects/cellevolution/live/genes/$gene/hyphy_labelled.nwk" --branches Foreground --multiple-hits Double+Triple --srv Yes --output "/rds/general/user/cm1118/projects/cellevolution/live/genes/$gene/${gene}_ABSREL.json" >> "/rds/general/user/cm1118/projects/cellevolution/live/genes/$gene/${gene}_ABSREL.log"
fi

