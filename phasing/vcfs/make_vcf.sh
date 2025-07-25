#! /bin/bash

module load tabix/0.2.6

working_dir="/net/eichler/vol28/projects/denovo_variation/nobackups/combined_batches"

while read sample; do 
	dnm_file=$(ls ../denovo_files/"$sample".denovo.txt);
	vcf=$(ls "$working_dir"/snv_calling/raw_calls/"$sample".hifi_denovo_calls.tsv);
	family=$(echo "$sample" | cut -d "_" -f1);
	header=$(echo "#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT" $family"_fa" $family"_mo" $sample);
	cat vcf_header.txt > $sample.vcf;
	echo "$header" | sed 's/ /\t/g' >> $sample.vcf;
	cut -f6 $dnm_file | tail -n+2 > dnm_temp.txt;
	grep -f dnm_temp.txt $vcf >> $sample.vcf;
	bgzip $sample.vcf;
	tabix $sample.vcf.gz;
done < "$working_dir"/manifests/all_children.txt 

rm dnm_temp.txt
