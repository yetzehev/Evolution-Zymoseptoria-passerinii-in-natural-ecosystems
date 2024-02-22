#!/bin/bash
 
#SBATCH --job-name=call #Give your job a name. 
#SBATCH --nodes=1#Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4 #Multithreading.
#SBATCH --time=96:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=32G #Total Memory per node to use for the job
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=rojas@evolbio.mpg.de #Email for notifications from previous line
#SBATCH --partition=standard #Request a specific partition for the resource allocation.


####################################################################
#
#      Alignment and sort using  Zpa796 reference genome
# 
####################################################################
##............................................
##      PROJECT ORGANIZATION
## Define your variables and substitute values 
##............................................

## Variables are defined in capital letters
PROJECT=SSLB_model
WORKDIR=/home/rojas/septoriasis_model/00_SNVcalling
META=${WORKDIR}/meta
BIN=${WORKDIR}/bin # This directory contains this script

## INPUTS
FASTQ_TRIM=/home/rojas/septoriasis_model/00_SNVcalling/trimmed ## input reads
REF=/home/rojas/ref_genomes/Zpa796/IPAnp/assembly/Z796IPnp.fasta ## reference genome


#................................................
## SOFTWARE
#Paths to the software
#ngm  and java are installed and accesible to all user
#................................................
NGM=/data/biosoftware/ngm/ngm/ngm
PICARDTOOLS=/data/biosoftware/Picard/v2.23.5
NCORES="1"
RAM=36 # RAM in gigabases

#................................................
#
#	Create ouput directories
#.................................................
## OUTPUTS
#mkdir ${WORKDIR}/bam
BAM=${WORKDIR}/bam

mkdir ${WORKDIR}/gvcf
GVCF=${WORKDIR}/gvcf


## Create reference index
#samtools faidx $REF

## Create a dictionary for the reference genome
samtools dict ${REF} > /home/rojas/ref_genomes/Zpa796/IPAnp/assembly/Z796IPnp.dict

## Create a list of files
ls ${FASTQ_TRIM} |  grep fastq | egrep 'R1|R2'|sed 's/.bbduk_R1.fastq.gz//g' |sed 's/.bbduk_R2.fastq.gz//g' | sort -r |uniq > ${META}/samplelist_${PROJECT}.trim.txt

############################################################
##	Create a script for calling 
##	aligning and sorting each sample
############################################################
##Set array values
cd ${WORKDIR}
NAME=($(cat ${META}/samplelist_${PROJECT}.trim.txt))
NJOB=($(seq 401 512)) ## As many paired samples: 112

## Open loop for creating individual scripts
cd ${BIN}
for i in "${!NJOB[@]}"; do

echo -e '#!/bin/bash' > ${NJOB[$i]}.snvCall.sh

#.........................................................
#
#	First module
#	Aligning reads against reference
#.........................................................
## Aligning paired reads
echo "#Aligning paired ${NAME[$i]}" >> ${NJOB[$i]}.snvCall.sh
echo "${NGM} -r ${REF} -1 ${FASTQ_TRIM}/${NAME[$i]}.bbduk_R1.fastq.gz  -2 ${FASTQ_TRIM}/${NAME[$i]}.bbduk_R2.fastq.gz -o ${BAM}/${NAME[$i]}.paired.bam -t ${NCORES} -b --rg-id ${NAME[$i]} --rg-sm ${NAME[$i]} --rg-pl illumina --rg-pu ${PROJECT}"  >> ${NJOB[$i]}.snvCall.sh

## Aligning unpaired reads
echo "#Aligning unpaired ${NAME[$i]}"  >> ${NJOB[$i]}.snvCall.sh
echo    "${NGM} -r ${REF} -q ${FASTQ_TRIM}/${NAME[$i]}.bbduk_R1.fastq.gz -o ${BAM}/${NAME[$i]}.R1.unpaired.bam -t ${NCORES} -b"  >> ${NJOB[$i]}.snvCall.sh
echo    "${NGM} -r ${REF} -q ${FASTQ_TRIM}/${NAME[$i]}.bbduk_R2.fastq.gz -o ${BAM}/${NAME[$i]}.R2.unpaired.bam -t ${NCORES} -b" >> ${NJOB[$i]}.snvCall.sh

#.........................................................
#	Test second module: 
#	Procesing reads: merge, sort, add groups and clean
#.........................................................
## Merging bam files
echo "#Processing for ${NAME[$i]}"  >> ${NJOB[$i]}.snvCall.sh
echo "java -jar ${PICARDTOOLS}/picard.jar MergeSamFiles I=${BAM}/${NAME[$i]}.paired.bam I=${BAM}/${NAME[$i]}.R1.unpaired.bam I=${BAM}/${NAME[$i]}.R2.unpaired.bam O=${BAM}/${NAME[$i]}.merged.bam VALIDATION_STRINGENCY=LENIENT COMPRESSION_LEVEL=0 " >> ${NJOB[$i]}.snvCall.sh
echo "samtools view -h ${BAM}/${NAME[$i]}.merged.bam > ${BAM}/${NAME[$i]}.merged.sam"  >> ${NJOB[$i]}.snvCall.sh

## Sort sam files
echo "#sort sample ${NAME[$i]}" >> ${NJOB[$i]}.snvCall.sh
echo "samtools view -hSb ${BAM}/${NAME[$i]}.merged.sam > ${BAM}/${NAME[$i]}.merged.v0.bam"  >> ${NJOB[$i]}.snvCall.sh
echo "java -jar ${PICARDTOOLS}/picard.jar SortSam I=${BAM}/${NAME[$i]}.merged.v0.bam O=${BAM}/${NAME[$i]}.merged.v1.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT COMPRESSION_LEVEL=0"  >> ${NJOB[$i]}.snvCall.sh

#.........................................................
## Module 3. Add or replace groups 
#.........................................................
echo "#AddOrReplaceGroups for ${NAME[$i]}" >> ${NJOB[$i]}.snvCall.sh
echo "java -jar ${PICARDTOOLS}/picard.jar AddOrReplaceReadGroups I=${BAM}/${NAME[$i]}.merged.v1.bam O=${BAM}/${NAME[$i]}.merged.v2.bam RGID=${NAME[$i]} RGLB=${PROJECT} RGPL=ILLUMINA RGPU=${PROJECT} RGSM=${NAME[$i]} SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT COMPRESSION_LEVEL=0" >> ${NJOB[$i]}.snvCall.sh
echo "java -jar ${PICARDTOOLS}/picard.jar CleanSam I=${BAM}/${NAME[$i]}.merged.v2.bam O=${BAM}/${NAME[$i]}.merged.v3.bam VALIDATION_STRINGENCY=LENIENT" >> ${NJOB[$i]}.snvCall.sh

#.........................................................
## Optional. Mark Duplicates
#.........................................................
echo	"java -jar ${PICARDTOOLS}/picard.jar MarkDuplicates I=${BAM}/${NAME[$i]}.merged.v3.bam O=${BAM}/${NAME[$i]}.bam  M=${BAM}/${NAME[$i]}.marked_dup_metrics.txt" >> ${NJOB[$i]}.snvCall.sh


#................................................
## Module 4.	Run GATK/ Haplotype Caller
#................................................
echo "#Run haplotype caller for ${NAME[$i]}" >> ${NJOB[$i]}.snvCall.sh
echo "samtools index -b ${BAM}/${NAME[$i]}.bam" >> ${NJOB[$i]}.snvCall.sh
echo "#Call GATK HaplotypeCaller for sample ${NAME}" >> ${NJOB[$i]}.snvCall.sh
echo 'gatk --java-options "-Xmx'${RAM}'g" HaplotypeCaller \' >> ${NJOB[$i]}.snvCall.sh
echo	'-R '${REF}'  \' >> ${NJOB[$i]}.snvCall.sh
echo	'-I '${BAM}/${NAME[$i]}'.bam \' >> ${NJOB[$i]}.snvCall.sh
echo	'-O '${GVCF}/${NAME[$i]}'.gvcf.vcf \' >> ${NJOB[$i]}.snvCall.sh
echo	'-ploidy 1 \' >> ${NJOB[$i]}.snvCall.sh
echo	'--do-not-run-physical-phasing \' >> ${NJOB[$i]}.snvCall.sh
echo	'--native-pair-hmm-threads '${NCORES}' \' >> ${NJOB[$i]}.snvCall.sh
echo	"-ERC GVCF  #Emmitting reference confidence scores, gvcf format"  >> ${NJOB[$i]}.snvCall.sh

#.........................................................
## Remove intermediate files
#.........................................................
echo "Remove intermediate files" >> ${NJOB[$i]}.snvCall.sh
echo	"rm ${BAM}/${NAME[$i]}.paired.bam" >> ${NJOB[$i]}.snvCall.sh
echo	"rm ${BAM}/${NAME[$i]}.merged.sam" >> ${NJOB[$i]}.snvCall.sh
echo	"rm ${BAM}/${NAME[$i]}.merged.bam" >> ${NJOB[$i]}.snvCall.sh
echo	"rm ${BAM}/${NAME[$i]}.merged.v0.bam" >> ${NJOB[$i]}.snvCall.sh
echo	"rm ${BAM}/${NAME[$i]}.merged.v1.bam" >> ${NJOB[$i]}.snvCall.sh
echo	"rm ${BAM}/${NAME[$i]}.merged.v2.bam" >> ${NJOB[$i]}.snvCall.sh
echo	"rm ${BAM}/${NAME[$i]}.merged.v3.bam" >> ${NJOB[$i]}.snvCall.sh
echo	"rm ${BAM}/${NAME[$i]}.R1.unpaired.bam" >> ${NJOB[$i]}.snvCall.sh
echo	"rm ${BAM}/${NAME[$i]}.R2.unpaired.bam" >> ${NJOB[$i]}.snvCall.sh

done < ${META}/samplelist_${PROJECT}.trim.txt


##############################################
## Submit an array
##############################################
#Make a list of jobs
ls *.snvCall.sh > list.jobs.txt

#Create submitting jobs
# Create the submission script
echo -e '#!/bin/bash' > array_batch.sh
echo -e '# get sample' >> array_batch.sh
echo -e 'printf -v SLURM_ARRAY_TASK_ID_three "%03d" $SLURM_ARRAY_TASK_ID' >> array_batch.sh # $SLURM_ARRAY_TASK_ID is a variable that takes a number among the array interval the user specify. But I need to make it 3-digits to use grep and avoid multiple matches, so this step.
echo -e 'sample_ID=$(cat list.jobs.txt | grep "$SLURM_ARRAY_TASK_ID.snvCall.sh")' >> array_batch.sh    #CHANGE NAME OF THE JOB
echo -e '# run the command' >> array_batch.sh
echo -e 'bash $sample_ID' >> array_batch.sh


#Submit job
sbatch --job-name=snvCall --array=401-512 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=96:00:00 --mem=36G --error=job_call_%A_%a.err --output=job_call_%A_%a.out --mail-type=FAIL --mail-user=rojas@evolbio.mpg.de --partition=global array_batch.sh
