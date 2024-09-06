#!/bin/bash

#SBATCH --job-name=telom
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --time=70:00:00
#SBATCH --mem=20G
#SBATCH --error=telom.sterr
#SBATCH --output=telom.stout
#SBATCH --mail-type=NONE
#SBATCH --mail-user=rojas@evolbio.mpg.de
#SBATCH --partition=standard

#   ----------------
# |  INPUTS         |
#   ----------------


home="/home/rojas/Zpasserinii/1_PacBio_assembly/"
fasta_file="/home/rojas/Zpasserinii/1_PacBio_assembly/1_assemblies_edited/IPAnp/assembly/Z796IPnp.fasta"
work_dir="$home/4_evaluateTelomeres/"
project=Z796

## Source of FindTelomers.py 
## https://github.com/JanaSperschneider/FindTelomeres

## SOFTWARE, SCRIPTS
FindTelomers=/home/rojas/Zpasserinii/1_PacBio_assembly/1_assembly/bin/FindTelomeres.py

cd ${work_dir}
python3 ${FindTelomers} ${fasta_file} > ${project}.txt
