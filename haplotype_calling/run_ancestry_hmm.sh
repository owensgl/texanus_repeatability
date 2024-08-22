#cat ../meta/samplelists.txt | parallel -j 10 bash run_ancestry_hmm.sh
i=$1
pop=$(echo $i | cut -f 2 -d '.');
gen=$(echo $i | cut -f 3 -d '.');
let gen_plus=$gen+2
prefix=$(echo $i | sed 's/.samples.txt//g')
for parent in deb ann1 ann2 ann3
do
  cd $parent
  bcftools view -S ../../meta/$i ../../vcf/texanus.freebayes.parentagesnps.vcf | perl ../../bin/vcf2ancestryhmm_tex.pl ../../vcf/texanus.freebayes.BC1.10AF.50missing.5dp.bi.15-35MAF.snpparentage.$parent.txt ../../vcf/texanus.freebayes.5dp.parentagesnps.cmlocations.txt $prefix ../../meta/$i
  proportion_ancestry=$(cat $prefix.ancestry_proportion.txt | cut -f 2 -d " ")
  other_proportion=$(echo "1-$proportion_ancestry" | bc)
  ../../bin/Ancestry_HMM-1.0/src/ancestry_hmm -i $prefix.ancestry_hmminput.txt -s ../../meta/ancestry_hmm_ploidy/$i -a 2 $proportion_ancestry $other_proportion -p 1 10000 0.75 -p 0 -$gen_plus 0.25 -min $gen_plus 2> $prefix.ancestryhmm.log
  cd ..
done
