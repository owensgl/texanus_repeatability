for pop in `ls ../meta/ | grep samples | sed s/.txt//g`
do
	vcftools --vcf texanus.freebayes.parentagesnps.deb.vcf --keep ../meta/$pop.txt --geno-r2 --out $pop.deb
	vcftools --vcf texanus.freebayes.parentagesnps.ann1.vcf --keep ../meta/$pop.txt --geno-r2 --out $pop.ann1
	vcftools --vcf texanus.freebayes.parentagesnps.ann2.vcf --keep ../meta/$pop.txt --geno-r2 --out $pop.ann2
	vcftools --vcf texanus.freebayes.parentagesnps.ann3.vcf --keep ../meta/$pop.txt --geno-r2 --out $pop.ann3
done
