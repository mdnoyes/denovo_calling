#! /bin/bash


working_dir="/net/eichler/vol28/projects/denovo_variation/nobackups/combined_batches"

while read sample; do 
	file=$(ls "$working_dir"/snv_calling/results/"$sample".snv_calls.tsv);
	echo "sample chr pos ref alt id" | sed 's/ /\t/g' > $sample.denovo.txt; 
	awk -v s="$sample" '{print s, $1, $2, $4, $5, $3}' $file | tail -n+2 | sed 's/ /\t/g' >> $sample.denovo.txt; 
done < "$working_dir"/manifests/all_children.txt
