#!/bin/bash

##create an environment to contain all relevant packages (multiqc, fastp, cutadapt)

#load packages
module load bcl2fastq/2.20.0.422
module load fastqc/0.11.9
module load BBMap/38.45
module load R/3.6.0  
##load miniconda first, change 'analysis' to whatever your conda environment name is
source activate analysis ## multiqc, fastp, cutadapt - packages to install on the conda environment

#load config file - MAKE SURE TO FILL IT CORRECTLY
source ./aux_01_config_file


##9 creating internal indices files and internal demultiplexing
sed 's/ /_/g' $path_to_int_index > int_index_1
sed 's/,/_/' int_index_1 > int_index_2
sed 's/,/\n^/' int_index_2 > int_index_3
sed 's/Gr/\>Gr/' int_index_3 > int_index_final
ls | grep "int_index_.*[1|2|3]$" | awk '{system("rm "$1)}'

for folder in $(ls -l |grep ^d | grep -v -E 'Reports|Stats' | awk '{print $9}')
  do
    cd $folder
    cat ../int_index_final | grep -A 1 $folder"_" > "int_index_"$folder 
    cd ..
  done 

for folder in $(ls -l |grep ^d | grep -v -E 'Reports|Stats' | awk '{print $9}')
  do 
    cd $folder
    barcode=$(ls int_index*)
    index=$(ls *only_index.fastq.gz)
    read1=$(ls *R1_repaired.fastq.gz)
    
    bsub -K -q gpu-short -n 1 -R "rusage[mem=32000]" -R "rusage[ngpus_physical=1]" -J $folder"_int_index" -o $folder"_int_index.out" -e $folder"_int_index.err" \
    cutadapt -e 0.15 --no-indels -g file:$barcode -o int-{name}.2.fastq -p int-{name}.1.fastq $index $read1 &
    cd ..
  done
wait

comp=$(find . -name "*_int_index.out" | awk '{system("cat "$1)}' | grep -c Successfully)
if (( $comp == $samples )); then echo "Successfully finished demultimplexing intenral indices"; else echo "Something went wrong with internal de-multiplexing, please check log files" && exit ; fi 

##10 - finding matching umis for each sample

for folder in $(ls -l |grep ^d | grep -v -E 'Reports|Stats' | awk '{print $9}')
  do
    cd $folder 
    for file in $(ls int-*.1.fastq | awk '{print $1}')
    do
      group=$(ls Gr*only_umi.fastq.gz)
      bsub -K -q gpu-long -n 1 -R "rusage[mem=64000]" -R "rusage[ngpus_physical=1]" -J $file"_umi_int_index" -o $file"_umi_int_index.out" -e $file"_umi_int_index.err" \
      repair.sh in1=$group in2=$file out1=ni.fastq.gz out2=$file"_umis.fastq.gz" repair overwrite=true &
    done
    rm ni.fastq.gz
    cd ..  
  done
wait

comp=$(find . -name "*_umi_int_index.out" | awk '{system("cat "$1)}' | grep -c Successfully)
samples=$(find . -name "int-*.1.fastq" | wc -l)
if (( $comp == $samples )); then echo "Successfully finished demultimplexing umis by internal indices"; else echo "Something went wrong with internal de-multiplexing of umis, please check log files" && exit ; fi 

##10 - arranging files for later work

mkdir ../experiments
mkdir ../experiments/umis
find . -name "int*.1.fastq" -exec cp {} "../experiments/" \;
find . -name "*_umis.fastq.gz" -exec cp {} "../experiments/umis/" \;
cd ../experiments
