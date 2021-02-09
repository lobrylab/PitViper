directory=$(find ../../../PitViper_development/data/essentials_genes_sampling/)
for file in $directory
do
    if test -f "$file"; then
        regex='.+\/([0-9]+)'
        [[ $file =~ $regex ]]
        sample_id=${BASH_REMATCH[1]}
        echo $sample_id
        snakemake --configfile config/configfile_shRNA.yaml --use-conda --cores 1 --config essential=../../../PitViper_development/data/essentials_genes_sampling/${sample_id}_essentials_genes.txt token=${sample_id}
    fi
done
