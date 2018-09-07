#!/bin/bash

#usage: ./check-snp-hap.sh line snpfile outfile

chr=`sed -n ${1}p ${2} | cut -f1`
pos=`sed -n ${1}p ${2} | cut -f2`
rs=`sed -n ${1}p ${2} | cut -f3`

bash -xv get-haplotype-1kg.sh ${chr} $(($pos - 50000)) $(($pos + 50000)) $pos $1

if [ -f "meth.haps/${1}.${pos}.meth.haps" ]
then
  R --vanilla  "--args cis.gibbs.num=${1};outfile=\"${3}\";psid=${1};posid=${pos}" < score.R
fi


#rm temp.${1}.*
