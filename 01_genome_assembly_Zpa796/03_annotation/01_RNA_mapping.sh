#!/bin/bash
 
#SBATCH --job-name=hisat2.bbduk2 #Give your job a name.
#SBATCH --nodes=1#Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4 #Multithreading.
#SBATCH --time=24:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=24G #Total Memory per node to use for the job
#sbatch --gridOptions="--time=24:00:00 --partition norm,b1"
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=rojas@evolbio.mpg.de #Email for notifications from previous line
#SBATCH --partition=fast #Request a specific partition for the resource allocation.


####################################################################
#
#       BRAKER (DB + RNA)
# 
#       Trimmed & Filtered Zpa796 reads
#
####################################################################


#............................................
#PROJECT ORGANIZATION
#Define your variables
#............................................
WORKDIR=/home/rojas/Zpasserinii/1_PacBio_assembly/6_braker_RNA/trimming/bbduk.test2 # Trimmed reads 
BIN=${WORKDIR}/bin2.bbduk # Contains this script
READS=${WORKDIR}/Zpa796.reads #Merged &trimmed reads
REFERENCE=/home/rojas/Zpasserinii/1_PacBio_assembly/1_assemblies_edited/IPAnp/assembly/Z796IPnp.fna
GTF=/home/rojas/Zpasserinii/1_PacBio_assembly/6_braker/BRAKER.IPnp/braker/braker.gtf
PROJECT=Z796bbduk2v3
EXTENSION=bbduk


#............................................
#SET YOUR ENVIRONMENT
# THIS SCRIPT INCLUDES MODIFICATIONS TO THE bbdul trimming. I added quality trimming
#............................................
#SOFTWARE
#Hisat is in the PATH
source /data/modules/python/python-anaconda3/etc/profile.d/conda.sh
hisat2=/data/biosoftware/hisat2/hisat2

#..........................
#Make directories
#..........................
cd ${WORKDIR}
#mkdir reference
mkdir bam


#............................
#Make indexes and references in the REFERENCE directory
#............................
#Copy reference genome and GTF file
cd reference 
#cp ${REFERENCE} .
#cp ${GTF} .

#Rename references
mv *.fa ${PROJECT}.fa
mv *.gtf ${PROJECT}.gtf

#Make index
samtools faidx ${PROJECT}.fa
#Create splice site list
hisat2_extract_splice_sites.py ${PROJECT}.gtf> ${PROJECT}.gtf.splice_sites
#Create exon list
hisat2_extract_exons.py ${PROJECT}.gtf > ${PROJECT}.gtf.exons
#Create index
${hisat2}/hisat2-build -p 4 --ss ${PROJECT}.gtf.splice_sites --exon ${PROJECT}.gtf.exons ${PROJECT}.fa ${PROJECT}
#Move to WORKDIR
cd ${WORKDIR}


#............................
#Make an script for aligning the filtered reads
#............................
cd ${WORKDIR}/${BIN}
#I have  20 fastq files that only contain Zpa796 reads, the reads for Hordeum were removed with fasq_screen
#********************
# No white spaces among the list of fastq reads
#********************
#Read mapping against the reference genome
echo -e '#!/bin/bash' > ${PROJECT}.hisat2.sh
echo "#Align all reads from Zpa796 and generate a single bam file" >> ${PROJECT}.hisat2.sh
echo "${hisat2}/hisat2 -p 4 --very-sensitive --dta -x ${WORKDIR}/reference/${PROJECT} -1 ${READS}/Z796Ax72h01_ME_${EXTENSION}_1_001.fastq.gz,${READS}/Z796Ax72h02_ME_${EXTENSION}_1_001.fastq.gz,${READS}/Z796Ax72h04_ME_${EXTENSION}_1_001.fastq.gz -2 ${READS}/Z796Ax72h01_ME_${EXTENSION}_2_001.fastq.gz,${READS}/Z796Ax72h02_ME_${EXTENSION}_2_001.fastq.gz,${READS}/Z796Ax72h04_ME_${EXTENSION}_2_001.fastq.gz -S ${WORKDIR}/bam/Zpa796.allReads.${EXTENSION}.hisat2.sam"  >> ${PROJECT}.hisat2.sh
echo "#Sort and indexing alignments (samtools)" >> ${PROJECT}.hisat2.sh
echo "samtools sort -o ${WORKDIR}/bam/Zpa796.allReads.${EXTENSION}.hisat2.bam ${WORKDIR}/bam/Zpa796.allReads.${EXTENSION}.hisat2.sam" >> ${PROJECT}.hisat2.sh
echo "samtools index ${WORKDIR}/bam/Zpa796.allReads.${EXTENSION}.hisat2.bam"  >> ${PROJECT}.hisat2.sh

#done

sbatch --job-name=hisat2 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=4 --time=96:00:00 --mem=36G --error=job_hisat2.err --output=job_hisat2.out --mail-type=FAIL --mail-user=rojas@evolbio.mpg.de --partition=global ${PROJECT}.hisat2.sh

