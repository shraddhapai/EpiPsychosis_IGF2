plink=/home/shraddhapai/software/plink/plink
$plink --bfile ../SUS19399.hg19 --exclude Exclude-SUS19399.hg19-1000G.txt --make-bed --out TEMP1
$plink --bfile TEMP1 --update-map Chromosome-SUS19399.hg19-1000G.txt --update-chr --make-bed --out TEMP2
$plink --bfile TEMP2 --update-map Position-SUS19399.hg19-1000G.txt --make-bed --out TEMP3
$plink --bfile TEMP3 --flip Strand-Flip-SUS19399.hg19-1000G.txt --make-bed --out TEMP4
$plink --bfile TEMP4 --reference-allele Force-Allele1-SUS19399.hg19-1000G.txt --make-bed --out SUS19399.hg19-updated
$plink --bfile SUS19399.hg19-updated --reference-allele Force-Allele1-SUS19399.hg19-1000G.txt --chr 1 --recode vcf --out SUS19399.hg19-updated-chr1 
$plink --bfile SUS19399.hg19-updated --reference-allele Force-Allele1-SUS19399.hg19-1000G.txt --chr 2 --recode vcf --out SUS19399.hg19-updated-chr2 
$plink --bfile SUS19399.hg19-updated --reference-allele Force-Allele1-SUS19399.hg19-1000G.txt --chr 3 --recode vcf --out SUS19399.hg19-updated-chr3 
$plink --bfile SUS19399.hg19-updated --reference-allele Force-Allele1-SUS19399.hg19-1000G.txt --chr 4 --recode vcf --out SUS19399.hg19-updated-chr4 
$plink --bfile SUS19399.hg19-updated --reference-allele Force-Allele1-SUS19399.hg19-1000G.txt --chr 5 --recode vcf --out SUS19399.hg19-updated-chr5 
$plink --bfile SUS19399.hg19-updated --reference-allele Force-Allele1-SUS19399.hg19-1000G.txt --chr 6 --recode vcf --out SUS19399.hg19-updated-chr6 
$plink --bfile SUS19399.hg19-updated --reference-allele Force-Allele1-SUS19399.hg19-1000G.txt --chr 7 --recode vcf --out SUS19399.hg19-updated-chr7 
$plink --bfile SUS19399.hg19-updated --reference-allele Force-Allele1-SUS19399.hg19-1000G.txt --chr 8 --recode vcf --out SUS19399.hg19-updated-chr8 
$plink --bfile SUS19399.hg19-updated --reference-allele Force-Allele1-SUS19399.hg19-1000G.txt --chr 9 --recode vcf --out SUS19399.hg19-updated-chr9 
$plink --bfile SUS19399.hg19-updated --reference-allele Force-Allele1-SUS19399.hg19-1000G.txt --chr 10 --recode vcf --out SUS19399.hg19-updated-chr10 
$plink --bfile SUS19399.hg19-updated --reference-allele Force-Allele1-SUS19399.hg19-1000G.txt --chr 11 --recode vcf --out SUS19399.hg19-updated-chr11 
$plink --bfile SUS19399.hg19-updated --reference-allele Force-Allele1-SUS19399.hg19-1000G.txt --chr 12 --recode vcf --out SUS19399.hg19-updated-chr12 
$plink --bfile SUS19399.hg19-updated --reference-allele Force-Allele1-SUS19399.hg19-1000G.txt --chr 13 --recode vcf --out SUS19399.hg19-updated-chr13 
$plink --bfile SUS19399.hg19-updated --reference-allele Force-Allele1-SUS19399.hg19-1000G.txt --chr 14 --recode vcf --out SUS19399.hg19-updated-chr14 
$plink --bfile SUS19399.hg19-updated --reference-allele Force-Allele1-SUS19399.hg19-1000G.txt --chr 15 --recode vcf --out SUS19399.hg19-updated-chr15 
$plink --bfile SUS19399.hg19-updated --reference-allele Force-Allele1-SUS19399.hg19-1000G.txt --chr 16 --recode vcf --out SUS19399.hg19-updated-chr16 
$plink --bfile SUS19399.hg19-updated --reference-allele Force-Allele1-SUS19399.hg19-1000G.txt --chr 17 --recode vcf --out SUS19399.hg19-updated-chr17 
$plink --bfile SUS19399.hg19-updated --reference-allele Force-Allele1-SUS19399.hg19-1000G.txt --chr 18 --recode vcf --out SUS19399.hg19-updated-chr18 
$plink --bfile SUS19399.hg19-updated --reference-allele Force-Allele1-SUS19399.hg19-1000G.txt --chr 19 --recode vcf --out SUS19399.hg19-updated-chr19 
$plink --bfile SUS19399.hg19-updated --reference-allele Force-Allele1-SUS19399.hg19-1000G.txt --chr 20 --recode vcf --out SUS19399.hg19-updated-chr20 
$plink --bfile SUS19399.hg19-updated --reference-allele Force-Allele1-SUS19399.hg19-1000G.txt --chr 21 --recode vcf --out SUS19399.hg19-updated-chr21 
$plink --bfile SUS19399.hg19-updated --reference-allele Force-Allele1-SUS19399.hg19-1000G.txt --chr 22 --recode vcf --out SUS19399.hg19-updated-chr22 

echo "bgzipping"
for chr in {1..22};do
	echo $chr
	bgzip SUS19399.hg19-updated-chr${chr}.vcf
done

rm TEMP*
