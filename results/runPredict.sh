## file:
##  runPredict.sh
## description:
##  example shell script to submit parallelized tasks to SGE cluster

## job name
#$ -N runPredict
## number of cores per job
#$ -pe smp 1
## memory allocation
#$ -l h_vmem=6G
## Make sure that the .e and .o file arrive in the working directory
#$ -cwd
## Merge the standard out and standard error to one file in one folder
#$ -j y
## Disable log output
#$ -o /dev/null
#$ -e /dev/null
## set up distributed jobs for nclust range, equal to nrow of param file
#$ -t 1-32400
## limit the number of simultaneous jobs
#$ -tc 100

SEEDFILE=$PWD/runPredict.param.txt
FIELD1=$(grep "^$SGE_TASK_ID " $SEEDFILE | cut -d' ' -f1)
FIELD2=$(grep "^$SGE_TASK_ID " $SEEDFILE | cut -d' ' -f2)
FIELD3=$(grep "^$SGE_TASK_ID " $SEEDFILE | cut -d' ' -f3)
FIELD4=$(grep "^$SGE_TASK_ID " $SEEDFILE | cut -d' ' -f4)
FIELD5=$(grep "^$SGE_TASK_ID " $SEEDFILE | cut -d' ' -f5)
FIELD6=$(grep "^$SGE_TASK_ID " $SEEDFILE | cut -d' ' -f6)
FIELD7=$(grep "^$SGE_TASK_ID " $SEEDFILE | cut -d' ' -f7)
FIELD8=$(grep "^$SGE_TASK_ID " $SEEDFILE | cut -d' ' -f8)
FIELD9=$(grep "^$SGE_TASK_ID " $SEEDFILE | cut -d' ' -f9)
FIELD10=$(grep "^$SGE_TASK_ID " $SEEDFILE | cut -d' ' -f10)
Rscript $PWD/runPredict.R $FIELD1 $FIELD2 $FIELD3 $FIELD4 $FIELD5 $FIELD6 $FIELD7 $FIELD8 $FIELD9 $FIELD10
