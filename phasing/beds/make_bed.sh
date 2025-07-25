#! /bin/bash

working_dir="/net/eichler/vol28/projects/denovo_variation/nobackups/combined_batches"

while read sample; do 
	file=$(ls "$working_dir"/snv_calling/results/"$sample".snv_calls.tsv); 
	tail -n+2 $file | awk '{print $1, $2-100000, $2+100000}' | sed 's/ /\t/g' > $sample.bed; 
done < "$working_dir"/manifests/all_children.txt
