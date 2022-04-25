#!/bin/bash

#READ FIRST
# Run the script from the home folder which contains the experiments folder 
# make sure you create a config file specific for the experiment which is name starts with "aux", and finishes with the "exp_name" - this file will hold the path to the initial seeds and the path to the relevant samples to the experiment
## fill in the experiment name (should be the same as appeares in the name of the config file
exp_name="gradual"

config=$(find "$(pwd)" -maxdepth 1 -name "aux_*"$exp_name)
#load config file - MAKE SURE TO FILL IT CORRECTLY
source $config

cd experiments
mkdir $exp_name 

##create an environment to contain all relevant packages (multiqc, fastp, cutadapt)

#load packages
module load bcl2fastq/2.20.0.422
module load fastqc/0.11.9
module load BBMap/38.45
module load R/4.0.0  

##load miniconda first, change 'analysis' to whatever your conda environment name is
source activate analysis ## multiqc, fastp, cutadapt - packages to install on the conda environment




##11 - clustering seeds
echo "Clustering seeds"
#creating a fastq file to contatin all relevant files
lines=$(cat $initial_Seeds)
echo -n cat > conc_initial_Seeds
echo -n zcat > conc_umis_Seeds

for i in $lines
do
  file=$(ls | grep $i)
  echo -n " "$file >> conc_initial_Seeds 
done
echo -n " > "$exp_name"_initial_Seeds.fastq" >> conc_initial_Seeds 
     
source conc_initial_Seeds

#creating a fastq file to contatin all relevant umi files
cd umis  
for i in $lines
do
  file=$(ls | grep $i)
  echo -n " "$file >> conc_umis_Seeds 
done
echo -n " > "$exp_name"_umis_Seeds.fastq" >> conc_umis_Seeds 
     
source conc_umis_Seeds
cd ..


##move all relevant files to the new folder of the experiment
touch files_to_move
cat $samples_to_cluster | awk '{system("ls | grep "$1" >> files_to_move")}'
folder="gradual"

sample_lines=$(cat files_to_move)
for i in $sample_lines
  do
    mv $i $exp_name"/"$i
done
mv $exp_name"_initial_Seeds.fastq" $exp_name"/"$exp_name"_initial_Seeds.fastq" 
rm files_to_move

cd $exp_name
bsub -K -q gpu-long -n 16 -R "rusage[mem=16000]" -R "rusage[ngpus_physical=1]" -J cd_hit_seeds -o cd_hit_seeds.out -e cd_hit_seeds.err \
cd-hit-est -i $exp_name"_initial_Seeds.fastq" -o $exp_name"_initial_seeds.clst" -c 0.95 -n 5 -d 0 -M 0 -r 0
echo "finished job1"
comp=$(find . -name "cd_hit_seeds.out" | awk '{system("cat "$1)}' | grep -c Successfully)
if (( $comp == 1 )); then echo "Successfully finished clustering initial seeds"; else echo "Something went wrong with initial seeds clustering, please check log file" && exit ; fi


##12cluster_2d
echo "Clustering samples against seeds"
sample_lines=$(cat $samples_to_cluster)
for i in $sample_lines
  do
    file_name=$(ls | grep $i)
    bsub -K -q gpu-short -n 4 -R "rusage[mem=10000]" -R "rusage[ngpus_physical=1]" -J $i"_cluster" -o $i"_cluster.out" -e $i"_cluster.err" \
    cd-hit-est-2d -i $exp_name"_initial_seeds.clst" -i2 $file_name -o $file_name"_VSall.clst" -c 0.95 -n 5 -d 0 -M 0 -r 0 -T 0 &
done
wait

num_samples=$(cat $samples_to_cluster | wc -l)
comp=$(find . -name "*_cluster.out" | awk '{system("cat "$1)}' | grep -c Successfully)

if (( $comp == $num_samples )); then echo "Successfully finished clustering_2d"; else echo "Something went wrong with 2d clustering, please check log file" && exit ; fi

for file in $(ls | grep VSall.*clstr$  | awk '{print $1}')
do
  clstr2txt.pl $file > $file"_df.txt" &
done
wait
comp=$(find . -name "*df.txt" | wc -l)
if (( $comp == $num_samples )); then echo "Successfully finished creating cluster_2d text files"; else echo "Something went wrong with creating the cluster2d text files, please check log file" && exit ; fi

##run R_code

for file in $(ls *clstr_df.txt | awk '{print $1}')
  do 
    bsub -q gpu-short -n 1 -R "rusage[mem=120000]" -R "rusage[ngpus_physical=1]" -J "$file"_R -o "$file"_R.out -e "$file"_R.err \
    Rscript $r_Script $file
  done



