#!/bin/bash

#SBATCH --job-name=array
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=96:00:00
#SBATCH --mem=10G
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --partition=standard 

##############################################
## Submit an array
##############################################
WORKDIR=/home/rojas/septoriasis_model/02_comp_genomics/02_TEs/REPET2.5

cd ${WORKDIR}
#Make a list of jobs
ls *REPET_pipeline.sh> list.jobs.txt

#Create submitting jobs
# Create the submission script
echo -e '#!/bin/bash' > array_batch.sh
echo -e '# get sample' >> array_batch.sh
echo -e 'printf -v SLURM_ARRAY_TASK_ID_three "%03d" $SLURM_ARRAY_TASK_ID' >> array_batch.sh # $SLURM_ARRAY_TASK_ID is a variable that takes a number among the array interval the user specify. But I need to make it 3-digits to use grep and avoid multiple matches, so this step.
echo -e 'sample_ID=$(cat list.jobs.txt | grep "$SLURM_ARRAY_TASK_ID.REPET_pipeline.sh")' >> array_batch.sh    #CHANGE NAME OF THE JOB
echo -e '# run the command' >> array_batch.sh
echo -e 'bash $sample_ID' >> array_batch.sh


#Submit job
sbatch --job-name=REPET2.5 --array=110-140 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=96:00:00 --mem=36G --error=job_call_%A_%a.err --output=job_call_%A_%a.out --mail-type=FAIL --mail-user=rojas@evolbio.mpg.de --partition=global array_batch.sh
