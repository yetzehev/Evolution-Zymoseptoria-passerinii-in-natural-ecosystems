#!/bin/bash

#SBATCH --job-name=REPET
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --time=336:00:00
#SBATCH --mem=64G
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --partition=global
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=rojas@evolbio.mpg.de #Email for notifications from previous line
#SBATCH --partition=global #Request a specific partition for the resource allocation.

# Reference:
# Cécile Lorrain, Alice Feurtey, Mareike Möller, Janine Haueisen, Eva Stukenbrock
# Dynamics of transposable elements in recently diverged fungal pathogens: lineage-specific transposable element content and efficiency of genome defenses
# G3 Genes|Genomes|Genetics, 
# Volume 11, Issue 4, April 2021, jkab068
# https://doi.org/10.1093/g3journal/jkab068

# Marco Guerreiro
# 07 March 2022
# REPET v3.0


#############################################
#### Configuration - Path to directories #### 

# Path to genome assembly directory and TE prediction directory (where config folder, this bash file and outputs will be)
# Genome asssembly files must be .fna
GENOMESFOLDER=/home/rojas/Zpasserinii/1_PacBio_assembly/5_Repet/REPET3.IPnp/assembly
TE_PREDICTION=/home/rojas/Zpasserinii/1_PacBio_assembly/5_Repet/REPET3.IPnp
#### Configuration - Genome assembly file name and project name #### 
# Genome assembly file name without path or extension (must be .fna)
# Project name must be maximum 8 characters
# For each project name, REPET creates an SQL table that cannot be overwritten. Repeating calculations with a previous project name will fail. Tables can be removed from the database.
SPECIES=Z796IPnp    
PROJECT=Z796IPnp    
        ########
#############################################


# Load modules
. config_files/setEnv.sh 


### REPET Pipeline
cd $TE_PREDICTION

mkdir $PROJECT
mkdir $PROJECT/TEdenovo
mkdir $PROJECT/TEdenovo/temp

# Renames fasta headers -> spaces and some characters might be a problem. Avoid space (" ") or symbols such as "=", ";", ":", "|" in headers.
# The genome fasta file must be “project_name.fa”
awk '/^>/{print ">'$PROJECT'_" ++i; next}{print}' < $GENOMESFOLDER/$SPECIES.fna > $TE_PREDICTION/$PROJECT/TEdenovo/$PROJECT.fa

# Make a link (ln -s) to access the databanks used in similarity based classification
ln -s /data/biosoftware/repet/hmm/ProfilesBankForREPET_Pfam32.0_GypsyDB.hmm $PROJECT/TEdenovo/
ln -s /data/biosoftware/repet/RepBase/RepBase20.05_REPET.embl/repbase20.05_aaSeq_cleaned_TE.fa $PROJECT/TEdenovo/
ln -s /data/biosoftware/repet/RepBase/RepBase20.05_REPET.embl/repbase20.05_ntSeq_cleaned_TE.fa $PROJECT/TEdenovo/


### Run TEdenovo pipeline

# get config file and update it with the directories path
cp config_files/TEdenovo.cfg $PROJECT/TEdenovo/

sed -i 's/project_name: /project_name: '$PROJECT'/g' $PROJECT/TEdenovo/TEdenovo.cfg
sed -i 's|project_dir: |project_dir: '$TE_PREDICTION'\/'$PROJECT'\/TEdenovo|g'  $PROJECT/TEdenovo/TEdenovo.cfg
sed -i 's|tmpDir: |tmpDir: '$TE_PREDICTION'\/'$PROJECT'\/TEdenovo\/temp|g'  $PROJECT/TEdenovo/TEdenovo.cfg


## STEP 1
# Move to project directory
cd $TE_PREDICTION/$PROJECT/TEdenovo/

# run TEdenovo pipeline
$REPET/launch_TEdenovo.py -P $PROJECT -f MCL > denovo.txt

# Important output files:
# Number of consensus before filtering:
# PROJECT_Blaster_Grouper_Map_consensus.fa
# PROJECT_Blaster_Piler_Map_consensus.fa
# PROJECT_Blaster_Recon_Map_consensus.fa
# TEdenovo metrics: 
# PROJECT_Blaster_GrpRecPil_Map_TEclassif_Filtered/PROJECT_sim_denovoLibTEs_PC_filtered.classif_stats.txt

# Move to project main folder
cd $TE_PREDICTION/


### Run TEannot pipeline
mkdir $PROJECT/TEannot
mkdir $PROJECT/TEannot/temp

# get config file and update it with the directories path
cp $TE_PREDICTION/config_files/TEannot.cfg $TE_PREDICTION/$PROJECT/TEannot/

sed -i 's/project_name: /project_name: '$PROJECT'_annot/g' $TE_PREDICTION/$PROJECT/TEannot/TEannot.cfg
sed -i 's|project_dir: |project_dir: '$TE_PREDICTION'\/'$PROJECT'\/TEannot|g' $TE_PREDICTION/$PROJECT/TEannot/TEannot.cfg
sed -i 's|tmpDir: |tmpDir: '$TE_PREDICTION'\/'$PROJECT'\/TEannot\/temp|g' $TE_PREDICTION/$PROJECT/TEannot/TEannot.cfg
sed -i 's/classif_table_name:/classif_table_name: '$PROJECT'/g' $TE_PREDICTION/$PROJECT/TEannot/TEannot.cfg

# Move to TEannot directory
cd $TE_PREDICTION/$PROJECT/TEannot/

# Link outputs from TEdenovo
ln -s $TE_PREDICTION/$PROJECT/TEdenovo/"$PROJECT"_Blaster_GrpRecPil_Map_TEclassif_Filtered/"$PROJECT"_sim_denovoLibTEs_filtered.fa "$PROJECT"_annot_refTEs.fa
ln -s $TE_PREDICTION/$PROJECT/TEdenovo/"$PROJECT".fa "$PROJECT"_annot.fa

## Step 1: The first step prepares all the data banks required in the next steps
$REPET/TEannot.py -P "$PROJECT"_annot -C TEannot.cfg -S 1

# Tables on MySQL: PROJECT_annot_chr_seq ; PROJECT_annot_chk_seq ; PROJECT_annot_chk_map ; PROJECT_annot_refTEs_seq ; PROJECT_annot_refTEs_map
# Files created: PROJECT_annot_db

## Step 2: aligns the reference TE sequences on each genomic chunk via BLASTER (high sensitivity,ls  followed by MATCHER) AND/OR REPEATMASKER (cutoff at 200) AND/OR CENSOR (high sensitivity)
$REPET/TEannot.py -P "$PROJECT"_annot -C TEannot.cfg -S 2 -a BLR -v 2
$REPET/TEannot.py -P "$PROJECT"_annot -C TEannot.cfg -S 2 -a RM  -v 2
$REPET/TEannot.py -P "$PROJECT"_annot -C TEannot.cfg -S 2 -a CEN -v 2

## Step 2 bis: idem to step 2 on randomized sequences to generate filter threshold
$REPET/TEannot.py -P "$PROJECT"_annot -C TEannot.cfg -S 2 -a BLR -r -v 2
$REPET/TEannot.py -P "$PROJECT"_annot -C TEannot.cfg -S 2 -a RM  -r -v 2
$REPET/TEannot.py -P "$PROJECT"_annot -C TEannot.cfg -S 2 -a CEN -r -v 2

#Files: PROJECT_annot_Tedetect ; PROJECT_annot_TEdetect_rnd

#Step3: filters and combines the high-scoring segment pairs (HSPs) obtained at step 2, i.e. the TE annotations
$REPET/TEannot.py -P "$PROJECT"_annot -C TEannot.cfg -S 3 -c BLR+RM+CEN -v 2 

#Tables MySQL: PROJECT_annot_chk_allTEs_path ; PROJECT_annot_chr_allTEs_path
#Files: PROJECT_annot_Tedetect/Comb

#Step 7: performs successive procedures on the MySQL tables such as removal of redundant TE, removal of SSR annotations included into TE annotations and "long join procedure" 
$REPET/TEannot.py -P "$PROJECT"_annot -C TEannot.cfg -S 7 -v 2
 
#Tables MySQL: PROJECT_annot_chk_allTEs_nr_path ; PROJECT_annot_chk_allTEs_nr_noSSR_path ; PROJECT_annot_chk_allTEs_nr_noSSR_join_path

# STEP 3
# Get the metrics of the genome
$REPET/GiveInfoFasta.py -i "$PROJECT"_annot.fa

# Get genome length
GENOME_LENGTH=$(grep -o -P '(?<=lengths\): ).*(?= bp)' "$PROJECT"_annot.fa.stats)

# STEP 4
# Extract information from the MySQL datasets and parse them into tables and text files (-g for genome size in pb)
$REPET/PostAnalyzeTELib.py -a 3 -g "$GENOME_LENGTH" -p "$PROJECT"_annot_chr_allTEs_nr_join_path -s "$PROJECT"_annot_refTEs_seq 

# STEP 5
$REPET/GetSpecificTELibAccordingToAnnotation.py -i "$PROJECT"_annot_chr_allTEs_nr_join_path.annotStatsPerTE.tab -t "$PROJECT"_annot_refTEs_seq


# STEP 6
# Run an other complete TEannot (steps 1 to 8) on your original genome using this validated TEs library (here the full length fragment library which is similar to full length copies but just a bit more stringent).  


cd $TE_PREDICTION/
mkdir $PROJECT/FLFannot ## full length fragment library
mkdir $PROJECT/FLFannot/temp

cp $TE_PREDICTION/$PROJECT/TEdenovo/repbase20.05_aaSeq_cleaned_TE.fa  $TE_PREDICTION/$PROJECT/FLFannot/
cp $TE_PREDICTION/$PROJECT/TEdenovo/repbase20.05_ntSeq_cleaned_TE.fa  $TE_PREDICTION/$PROJECT/FLFannot/

# get config file and update it with the directories path
cp $TE_PREDICTION/config_files/TEannot.cfg $TE_PREDICTION/$PROJECT/FLFannot/FLFannot.cfg

sed -i 's/project_name: /project_name: '$PROJECT'_FLF/g' $TE_PREDICTION/$PROJECT/FLFannot/FLFannot.cfg
sed -i 's|project_dir: |project_dir: '$TE_PREDICTION'\/'$PROJECT'\/FLFannot|g' $TE_PREDICTION/$PROJECT/FLFannot/FLFannot.cfg
sed -i 's|tmpDir: |tmpDir: '$TE_PREDICTION'\/'$PROJECT'\/FLFannot\/temp|g' $TE_PREDICTION/$PROJECT/FLFannot/FLFannot.cfg
sed -i 's/classif_table_name:/classif_table_name: '$PROJECT'/g' $TE_PREDICTION/$PROJECT/FLFannot/FLFannot.cfg


cd $TE_PREDICTION/$PROJECT/FLFannot/

# Consensus library link (it is important to give new names to avoid the previous files to be replaced and lost)
ln -s $TE_PREDICTION/$PROJECT/TEannot/"$PROJECT"_annot_chr_allTEs_nr_join_path.annotStatsPerTE_FullLengthFrag.fa  $TE_PREDICTION/$PROJECT/FLFannot/"$PROJECT"_FLF_refTEs.fa
# Genome link
ln -s $TE_PREDICTION/$PROJECT/TEannot/"$PROJECT"_annot.fa $TE_PREDICTION/$PROJECT/FLFannot/"$PROJECT"_FLF.fa

# And launch all steps in one command line (TEannot.cfg is the same as previously)
$REPET/launch_TEannot.py -P "$PROJECT"_FLF -C FLFannot.cfg > FLF.txt 


# PostAnalysis
$REPET/PostAnalyzeTELib.py -a 3 -g "$GENOME_LENGTH" -p "$PROJECT"_FLF_chr_allTEs_nr_noSSR_join_path -s "$PROJECT"_FLF_refTEs_seq 

# concatenate the .gff3 files 
cat "$PROJECT"_FLF_GFF3chr/*.gff3 |grep -v "##" > "$PROJECT"_refTEs.gff

cd $TE_PREDICTION

