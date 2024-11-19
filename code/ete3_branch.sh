#!/bin/bash

# Load necessary modules and activate the environment
#module load anaconda3/personal
#source activate ete3_env

# Define important environment variables
export repo_dir=/rds/general/user/cm1118/projects/cellevolution/live/EvoCT_coding
export analysis_dir=/rds/general/user/cm1118/projects/cellevolution/live/genes
export temp_dir=/rds/general/user/cm1118/projects/cellevolution/ephemeral

# Get the gene from the list indexed by PBS_ARRAY_INDEX
gene=$(head -$PBS_ARRAY_INDEX $repo_dir/data/branch_genes.txt | tail -1)

# create temp directory for gene
if [ ! -d "$temp_dir/$gene" ]; then
    mkdir -p "$temp_dir/$gene"
fi 

# Check if the log file already exists and exit if it does
if [ -e "${gene}_nobranchlen.branch.log" ]; then
    echo "Log file already exists for $gene. Exiting."
    exit 0
fi

# In the fasta file, find lines that only have gaps (---...) and the name line before it
sed -n '/^[-]*$/{s/.*//;x;d;};x;p;${x;p;}' $analysis_dir/$gene/$gene.fasta | sed '/^$/d' > $analysis_dir/$gene/$gene-clean.fasta

# Remove stop codons
if [ -f "$analysis_dir/$gene/mixed_types/clustalo_default-none-none-fasttree_full/$gene-clean-aa.fasta.final_tree.used_alg.fa" ]; then
    java -jar ~/anaconda3/envs/hyphy/share/macse-2.07-0/macse_v2.07.jar -prog exportAlignment -align "$analysis_dir/$gene/mixed_types/clustalo_default-none-none-fasttree_full/$gene-clean-aa.fasta.final_tree.used_alg.fa" -codonForFinalStop NNN -codonForInternalStop NNN -out_NT "$analysis_dir/$gene/${gene}_nostop.fasta"
else
    java -jar ~/anaconda3/envs/hyphy/share/macse-2.07-0/macse_v2.07.jar -prog exportAlignment -align "$analysis_dir/$gene/$gene-clean.fasta" -codonForFinalStop NNN -codonForInternalStop NNN -out_NT "$analysis_dir/$gene/${gene}_nostop.fasta"
fi

# Double-check that the last line isn't a species name
tail -1 $analysis_dir/$gene/$gene-clean.fasta | if grep -q ">" - ; then sed -i '$ d' $analysis_dir/$gene/$gene-clean.fasta ; fi

# Assign remaining taxa to a variable and save them to a file
taxa=$(grep ">" $analysis_dir/$gene/$gene-clean.fasta | sed ':a;N;$!ba;s/\n/ /g' | sed 's/>//g')
echo $taxa | sed 's/ /\n/g' > $analysis_dir/$gene/$gene.taxa.txt

# Build a new species tree only containing the post-trimming taxa
ete3 mod --prune $taxa -t $repo_dir/data/hg38.30way.nh.nobranchlen.txt  > $analysis_dir/$gene/species_tree.nw

# Remove branch lengths (recommended by developers of codeML)
perl -ne '$_=~s/:[\d\.]+//g; print $_;' $analysis_dir/$gene/species_tree.nw > $analysis_dir/$gene/species_tree_nobranchlen.nw

# Check if the species tree is not empty, if empty, create an empty taxa file
if [ -s $analysis_dir/$gene/species_tree.nw ]
then
    echo "Species tree is not empty"
else
    rm $analysis_dir/$gene/$gene.taxa.txt
    touch $analysis_dir/$gene/$gene.taxa.txt
fi

# Find which species are available for branch marks
awk 'NR==FNR{A[$1];next}$2 && $3 in A' $analysis_dir/$gene/$gene.taxa.txt $repo_dir/data/branch_mark_dictionary.txt | awk '!seen[$1]++ {print $4}' | awk 'BEGIN { ORS = " " } { print }' > $analysis_dir/$gene/branch_marks.txt

# Select branch marks for testing
marks=$(cat $analysis_dir/$gene/branch_marks.txt)

# Select which positive selection models to run
models=(b_free b_neut)

# Specify which models to run log-likelihood ratio tests (LRT) between
tests=(b_free b_neut)

# Run positive selection detection
echo "****** Running branch model for $gene ******"
ete3 evol -t $analysis_dir/$gene/species_tree_nobranchlen.nw --alg $analysis_dir/$gene/${gene}_nostop.fasta -o $temp_dir/$gene/$gene-branch-new --models ${models[*]} --tests ${tests[*]} --cpu 32 --mark ${marks[*]} >> $analysis_dir/$gene/${gene}.branch.log
echo "Done with $gene"
