#!/bin/bash
#$ -cwd -V
#$ -N hsm -j y
#$ -P lindgren.prja -q long.qa
#$ -t 100001-100200:100
### -pe shmem 2

echo "************************************************************"
echo "SGE job ID: "$JOB_ID
echo "SGE task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "************************************************************"

## for testing:
## SGE_TASK_ID=100

#bash -xv ./check-snp-hap-t2d-random.sh $SGE_TASK_ID snp.sample wmcg.random
#bash -xv ./check-snp-hap-t2d-random.sh 50005 snp.sample wmcgs
##bash -xv ./check-snp-hap-t2d-random.sh ${SGE_TASK_ID} snp.sample.05.random snp.sample.05/${SGE_TASK_ID}.wmcg
#obesity gwas
#bash -xv ./check-snp-hap-t2d-random.sh ${SGE_TASK_ID} ob.af ob.results/${SGE_TASK_ID}.wmcg

cd /well/drong/hsm2014

for i in `seq ${SGE_TASK_ID} $((${SGE_TASK_ID}+100))`
do 
bash -xv get-haplotype-1kg-450.sh ${SGE_TASK_ID} 450k.auto.chrpos.txt
done

cd script_out

echo "************************************************************"
echo "Finished at: "`date`
echo "************************************************************"
exit 0
