#!/bin/bash
 
#SBATCH --job-name=braker #Give your job a name.
#SBATCH --nodes=1#Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4 #Multithreading.
#SBATCH --time=96:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=32G #Total Memory per node to use for the job
#sbatch --gridOptions="--time=24:00:00 --partition norm,b1"
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=rojas@evolbio.mpg.de #Email for notifications from previous line
#SBATCH --partition=standard #Request a specific partition for the resource allocation.

# Note: activate environment evaluate_myAssembly
## $ conda activate evaluate_myAssembly

#Project organization
project="Zpa796_assembly"
home="/home/rojas/Zpasserinii/1_PacBio_assembly/assemblies"
ipa="$home/ipa.fasta"

cd $home

#Extract contig names
assembly="ipa"
grep '>' $ipa > scaffolds.$assembly.txt
cat scaffolds.$assembly.txt | wc -l

#Change scaffolds name
for num in {1..31}
do
Newheader=">contig.$num"
echo $Newheader >> scaffolds.numbers.$assembly.txt
done

# merge names and new header.
paste scaffolds.$assembly.txt scaffolds.numbers.$assembly.txt > scaffolds.name.numbers.$assembly.txt

# replace with sed in the loop
while read line
do
linearray=( $line )
OldName=${linearray[0]}
NewName=${linearray[1]}
sed -i "s|$OldName|$NewName|g" $ipa 
done < scaffolds.name.numbers.$assembly.txt


