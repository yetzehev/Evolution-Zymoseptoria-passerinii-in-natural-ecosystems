#!/bin/bash
 
#SBATCH --job-name=vcftools #Give your job a name.
#SBATCH --nodes=1#Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4 #Multithreading.
#SBATCH --time=96:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=48G #Total Memory per node to use for the job
#sbatch --gridOptions="--time=24:00:00 --partition norm,b1"
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=rojas@evolbio.mpg.de #Email for notifications from previous line
#SBATCH --partition=global #Request a specific partition for the resource allocation.


####################################################################
##
##       SNV calling for Z. passerinii
## 	Filtering
##  	Check https://www.ddocent.com//filtering/	
##	Check https://speciationgenomics.github.io/filtering_vcfs/
####################################################################


##............................................
##      PROJECT ORGANIZATION
## Define your variables
##............................................
## Variables are defined in capital letters
PROJECT=SSLB_model

## SOFTWARE

VCFpath=/home/rojas/septoriasis_model/00_SNVcalling/out

cd ${VCFpath}

## Apply a three step filter. We are going to only keep variants that have been successfully genotyped in 50% of individuals, a minimum quality score of 30, and a minor allele count of 3.

vcftools --gzvcf ${VCFpath}/SSLB_model_v4.2.5.vcf.gz --max-missing 0.5 --mac 3 --minQ 30 --recode --recode-INFO-all --out ${VCFpath}/SSLB.g5mac3

##kept 2701967 out of a possible 6549866 Sites

##  Missingness on a per-individual basis
vcftools --gzvcf ${VCFpath}/SSLB_model_v4.2.5.vcf.gz --missing-indv  --out ${PROJECT}

## Check missingnes data per individual

#```
# awk '!/IN/' SSLB_model.imiss| cut -f5 > totalmissing
# gnuplot << \EOF 
# set terminal dumb size 120, 30
# set autoscale 
# unset label
# set title "Histogram of % missing data per individual"
# set ylabel "Number of Occurrences"
# set xlabel "% of missing data"
# set yr [0:100000]
# binwidth=0.01
# bin(x,width)=width*floor(x/width) + binwidth/2.0
# plot 'totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
# pause -1
# EOF
#```

## Several samples have more than 0.8 missing data



##.......................................
##	Remove indivudals with more tha 0.5 missing data
##.......................................

#awk '$5 > 0.5' SSLB_model.imiss | cut -f1 > lowDP.indv

## In total 26 samples of Z. passerinii, 4 samples of Z. ardabiliae, 2 samples of Z. pseudiotriotici and 20 of C. beticola 
vcftools --vcf ${VCFpath}/SSLB.g5mac3.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out ${VCFpath}/SSLB_model.g5mac3dplm

## Check Mean depth per site averaged across all individuals
vcftools --gzvcf  ${VCFpath}/SSLB_model.g5mac3dplm.recode.vcf --site-mean-depth --out ${VCFpath}/${PROJECT}.g5mac3dplm
vcftools --gzvcf  ${VCFpath}/SSLB_model.g5mac3dplm.recode.vcf --min-meanDP 10 --max-meanDP 25 --recode --recode-INFO-all --out ${VCFpath}/SSLB_model.g5mac3dplm.10-25x
##After filtering, kept 2247524 out of a possible 2701967 Sites


## Check missingness per site
vcftools --gzvcf ${VCFpath}/SSLB_model.g5mac3dplm.recode.vcf --missing-site --recode --out ${VCFpath}/${PROJECT}.g5mac3dplm



## Check missing-site per species
#awk '$2 == "Zar"' popmap > Zar.keep && awk '$2 == "Zbr"' popmap > Zbr.keep && awk '$2 == "Zpa"' popmap > Zpa.keep &&  awk '$2 == "Zps"' popmap > Zps.keep && awk '$2 == "Ztr"' popmap > Ztr.keep && awk '$2 == "Cb"' popmap > Cb.keep 

## Estimate missing data for loci in each population
#vcftools --vcf SSLB_model.g5mac3dplm.10-25x.recode.vcf --keep ../meta/Zar.keep --missing-site --out Zar
#vcftools --vcf SSLB_model.g5mac3dplm.10-25x.recode.vcf --keep ../meta/Zbr.keep --missing-site --out Zbr
#vcftools --vcf SSLB_model.g5mac3dplm.10-25x.recode.vcf --keep ../meta/Zpa.keep --missing-site --out Zpa
#vcftools --vcf SSLB_model.g5mac3dplm.10-25x.recode.vcf --keep ../meta/Zps.keep --missing-site --out Zps 
#vcftools --vcf SSLB_model.g5mac3dplm.10-25x.recode.vcf --keep ../meta/Ztr.keep --missing-site --out Ztr
#vcftools --vcf SSLB_model.g5mac3dplm.10-25x.recode.vcf --keep ../meta/Cb.keep --missing-site --out Cb

## We can combine the two files and make a list of loci about the threshold of 10% missing data to remove. Note this is double the overall rate of missing data.

cat Zar.lmiss  Zbr.lmiss  Zpa.lmiss  Zps.lmiss  Ztr.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci.cat

## Remove badloci.cat
vcftools --vcf SSLB_model.g5mac3dplm.10-25x.recode.vcf --exclude-positions badloci.cat --recode --recode-INFO-all --out SSLB_model.g5mac3dplm.10-25x.mis0.1
## After filtering, kept 1167977 out of a possible 2247524 Sites


## Check missingnes per site again
vcftools --vcf SSLB_model.g5mac3dplm.10-25x.mis0.1.recode.vcf --missing-site --out SSLB_model.g5mac3dplm.10-25x.mis0.1



