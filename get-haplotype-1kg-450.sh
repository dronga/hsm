#!/bin/bash
#!/bin/bash
TABIX=/apps/well/tabix/0.2.6/tabix
VCFTOOLS=/apps/well/vcftools/16102015/bin/vcftools
PLINK=/apps/well/plink/1.07/plink
CHAOS=/well/lindgren/alexd/chaos/xchaos
HAPLOVIEW=~/bin/Haploview.jar

chr=`sed -n ${1}p ${2} | cut -f2 -d' '`
pos=`sed -n ${1}p ${2} | cut -f3 -d' '`
cg=`sed -n ${1}p ${2} | cut -f1 -d' '`


#$TABIX -fh ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr${1}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ${1}:${2}-${3} > temp.$5.vcf
$TABIX -fh /well/lindgren/alexd/1000genomes/ALL.chr${chr}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ${chr}:$((${pos} - 100000))-$((${pos}+100000)) > temp.$1.vcf

grep -E 'SNP|#' temp.$1.vcf > temp.$1.snp.vcf

$VCFTOOLS --vcf temp.$1.snp.vcf --keep EUR.sample.txt --plink --out temp.$1

$PLINK --ped temp.$1.ped --map temp.$1.map --recodeHV --out temp.$1 --noweb

java -Xmx13g -jar $HAPLOVIEW -nogui -pedfile temp.$1.ped -info temp.$1.info -out temp.$1 -blockoutput 

$CHAOS motif vcf2fasta --i temp.$1.snp.vcf --o temp.$1.fasta -lf 1 -rf 1
#cut -f2 temp.$1.map > temp.$1.snp
#grep -f temp.$1.snp chr${1}.fasta | sed 's/    /\n/' > temp.$1.fasta
#grep -A1 -f temp.2.snp chr1.fasta | sed '/^--$/d' > temp.$1.fasta

./find-methyl-alex.pl temp.$1.fasta > temp.$1.methyl.list

./score-blocks-pos2.pl temp.$1.info temp.$1.methyl.list temp.$1.GABRIELblocks ${pos} ${1} ${cg}

## Calculate Dosage

$PLINK --ped temp.$1.ped --map temp.$1.map  --recodeA --out temp.$1 --noweb

R --vanilla "--args mfile=\"temp.$1.methyl.list\" rfile=\"temp.$1.raw\" sfile=\"hap.markers/$1.$pos.hap.markers\"  cg=\"${cg}\"" < EUR.score.R

rm temp.$1*
