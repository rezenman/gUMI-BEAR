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

##1 - bcl2fastq

folder=$(ls -l | grep ^d | awk '{print $9}') 
bsub -K -q gpu-medium -R "rusage[mem=16000]" -R "rusage[ngpus_physical=1]" -n 8 -J bcl2fastq -o bcl2fastq.out -e bcl2fastq.err \
bcl2fastq \
    -p 8 \
    --barcode-mismatches 1 --input-dir $folder"/Data/Intensities/BaseCalls/" \
    --runfolder-dir $folder \
    --intensities-dir $folder"/Data/Intensities" \
    --sample-sheet SampleSheet.csv \
    --output-dir Multiplexed \
    --ignore-missing-controls \
    --ignore-missing-bcls \
    --no-bgzf-compression \
    --no-lane-splitting \
    --use-bases-mask Y*,I*,N*,Y*

comp=$(find . -name "bcl2fastq.out" | awk '{system("cat "$1)}' | grep -c Successfully)
if (( $comp == 1 )); then echo "Successfully finished de-multiplexing and base-calling"; else echo "Something went wrong with the bcl2fastq script, please check log file" && exit ; fi


cd Multiplexed

##2 - fastqc

for folder in $(ls -l |grep ^d | grep -v -E 'Reports|Stats' | awk '{print $9}')
  do 
    cd $folder
    read1=$(ls|grep R1)
    read2=$(ls|grep R2)
    
    bsub -K -q gpu-short -n 1 -R "rusage[mem=8000]" -R "rusage[ngpus_physical=1]" -J $folder"_fastqc" -o $folder"_fastqc.out" -e $folder"_fastqc.err" \
    fastqc $read1 $read2 &
    cd ..
  done
 
wait

comp=$(find . -name "*fastqc.out" | awk '{system("cat "$1)}' | grep -c Successfully)
if (( $comp == $samples )); then echo "Successfully finished creating quality reports"; else echo "Something went wrong with the quality reports, please check log files" && exit ; fi
cd ..

##3 - MultiQC  

bsub -K -q gpu-short -n 1 -R "rusage[mem=16000]" -R "rusage[ngpus_physical=1]" -J MultiQC -o MultiQC.out -e MultiQC.err \
multiqc .

comp=$(find . -name "MultiQC.out" | awk '{system("cat "$1)}' | grep -c Successfully)
if (( $comp == 1 )); then echo "Successfully finished aggregating fastqc reports"; else echo "Something went wrong with the MultiQC script, please check log file" && exit ; fi

cd Multiplexed
##4 - trimming uneccessary basepairs  

for folder in $(ls -l |grep ^d | grep -v -E 'Reports|Stats' | awk '{print $9}')
  do 
    cd $folder
    R1=$(ls *R1_001.fastq.gz)
    R2=$(ls *R2_001.fastq.gz)
    
    bsub -K -q gpu-short -n 4 -R "rusage[mem=8000]" -R "rusage[ngpus_physical=1]" -J $folder"_fastp_trim" -o $folder"_fastp_trim.out" -e $folder"_fastp_trim.err" \
    fastp -i $R1 -I $R2 -o $folder"_only_barcode.fastq.gz" -O $folder"_index_umi.fastq.gz" -f 5 -F 5 -T 3 -w 10 -h $folder"_fastp_trim_report.html" &
    cd ..
  done
wait

comp=$(find . -name "*fastp_trim.out" | awk '{system("cat "$1)}' | grep -c Successfully)
if (( $comp == $samples )); then echo "Successfully finished trimming uneccessary base-pairs"; else echo "Something went wrong with the fasp trim script, please check log files" && exit ; fi



##5 - Quality filtering for R1 and R2

for folder in $(ls -l |grep ^d | grep -v -E 'Reports|Stats' | awk '{print $9}')
  do 
    cd $folder
    R1=$(ls *barcode.fastq.gz)
    R2=$(ls *index_umi.fastq.gz)
    
    bsub -K -q gpu-short -n 4 -R "rusage[mem=8000]" -R "rusage[ngpus_physical=1]" -J $folder"_fastp_quality_R1" -o $folder"_fastp_quality_R1.out" -e $folder"Gr1_fastp_quality_R1.err" \
    fastp -i $R1 -o $folder"_R1_filtered.fastq.gz" -q 25 -u 9 -n 0 -l 24 -h $folder"_fastp_R1_quality_report.html" &
    
    bsub -K -q gpu-short -n 4 -R "rusage[mem=8000]" -R "rusage[ngpus_physical=1]" -J $folder"_fastp_quality_R2" -o $folder"_fastp_quality_R2.out" -e $folder"_fastp_quality_R2.err" \
    fastp -i $R2 -o $folder"_R2_filtered.fastq.gz" -n 2 -l 18 -h $folder"_fastp_R2_quality_report.html" &
    
    cd ..
  done
wait

comp_1=$(find . -name "*fastp_quality_R1.out" | awk '{system("cat "$1)}' | grep -c Successfully)
comp_2=$(find . -name "*fastp_quality_R2.out" | awk '{system("cat "$1)}' | grep -c Successfully)

if (( $comp_1 == $samples )); then echo "Successfully filtered R1 reads"; else echo "Something went wrong with the R1 reads filtering, please check log files" && exit ; fi
if (( $comp_2 == $samples )); then echo "Successfully filtered R2 reads"; else echo "Something went wrong with the R2 reads filtering, please check log files" && exit ; fi


##6 - Repairing reads
for folder in $(ls -l |grep ^d | grep -v -E 'Reports|Stats' | awk '{print $9}')
  do 
    cd $folder
    R1=$(ls *R1_filtered.fastq.gz)
    R2=$(ls *R2_filtered.fastq.gz)
    
    bsub -K -q gpu-short -n 1 -R "rusage[mem=128000]" -R "rusage[ngpus_physical=1]" -J $folder"_repair" -o $folder"_repair.out" -e $folder"_repair.err" \
    repair.sh in1=$R1 in2=$R2 out1=$folder"_R1_repaired.fastq.gz" out2=$folder"_R2_repaired.fastq.gz" outs=$folder"_singles.fastq.gz" repair &
    
    cd ..
  done
wait

comp=$(find . -name "*_repair.out" | awk '{system("cat "$1)}' | grep -c Successfully)
if (( $comp == $samples )); then echo "Successfully finished repairing reads"; else echo "Something went wrong with repairing the reads, please check log files" && exit ; fi

##7 - seperating R2 reads to umi reads and internal indices reads

for folder in $(ls -l |grep ^d | grep -v -E 'Reports|Stats' | awk '{print $9}')
  do 
    cd $folder
    R2=$(ls *R2_repaired.fastq.gz)
    
    bsub -K -q gpu-short -n 10 -R "rusage[mem=32000]" -R "rusage[ngpus_physical=1]" -J $folder"_trim_umi" -o $folder"_trim_umi.out" -e $folder"_trim_umi.err" \
    cutadapt -u 8 -o $folder"_only_umi.fastq.gz" $R2 -j 0  &
    
    bsub -K -q gpu-short -n 10 -R "rusage[mem=32000]" -R "rusage[ngpus_physical=1]" -J $folder"_trim_index" -o $folder"_trim_index.out" -e $folder"_trim_index.err" \
    cutadapt -u -10 -o $folder"_only_index.fastq.gz" $R2 -j 0 &
    cd ..
  done
wait


comp_1=$(find . -name "*_trim_umi.out" | awk '{system("cat "$1)}' | grep -c Successfully)
comp_2=$(find . -name "*_trim_index.out" | awk '{system("cat "$1)}' | grep -c Successfully)

if (( $comp_1 == $samples )); then echo "Successfully seperated umi reads"; else echo "Something went wrong with the umi seperation, please check log files" && exit ; fi
if (( $comp_2 == $samples )); then echo "Successfully seperated index reads"; else echo "Something went wrong with the index seperation, please check log files" && exit ; fi


##8 - creating report

touch final_report
echo "sample_name" "barcode" "index" "umi" >> final_report
for folder in $(ls -l |grep ^d | grep -v -E 'Reports|Stats' | awk '{print $9}')
  do
     cd $folder
     index=$(zcat *only_index.fastq.gz | grep -c "^@")
     read1=$(zcat *R1_repaired.fastq.gz |grep -c "^@")
     umi=$(zcat *only_umi.fastq.gz |grep -c "^@")
     echo $folder $read1 $index $umi >> ../final_report
     cd ..
  done
echo "finished creating report"

