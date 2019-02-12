#!/bin/bash
# parse GTF
# end with table that has chrom,start,end, gene_id, gene_name
geneFile="/home/shraddhapai/genome_annotation/grch38/gencode.v26.annotation.gtf.genesonly.gtf"

cat $geneFile | cut -f 1,4,5,7 > first.txt
cat $geneFile | cut -f 9 | awk -F";" '{
	sub(/gene_id /,"",$1);
	gsub(/\"/,"",$1)
	sub(/gene_name/,"",$3); 
	gsub(/\"/,"",$3)
	print $1"\t"$3
}' > second.txt
paste first.txt second.txt > gencode.v26.annotation.gtf.genesonly.txt
