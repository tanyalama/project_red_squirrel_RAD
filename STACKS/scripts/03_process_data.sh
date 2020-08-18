#!/bin/bash
#run process_radtags and clone_filter

#BSUB -q long
#BSUB -W 48:00
#BSUB -R rusage[mem=16000]
#BSUB -n 4
#BSUB -R span[hosts=1]
#BSUB -e ./logs/03_process_data.err
#BSUB -oo ./logs/03_process_data.log

##Need to specify barcode file
barcode_file="./stacks_barcodes.txt"

#Run atropos to trim first two bases
for file in ./00_RAW_data/*.fq.gz
do
sample=${file##./00_RAW_data/}
atropos trim --cut 2 -se $file -o 01_trimmed_RAW_data/${sample}
done

wait

#Run process_radtags
module load stacks/2.3d

for file in ./01_trimmed_RAW_data/*1.fq.gz
do
read_name=${file%%1.fq.gz}
process_radtags -1 $file  -2 ${read_name}2.fq.gz -b $barcode_file -o 02_process_out -c -q -r --inline_null -e sbfI --bestrad
done

wait

#Move "remainder reads" into separate folder
mkdir ./02_process_out/remainder_reads
mv ./02_process_out/*.rem.* ./02_process_out/remainder_reads

#Run clone_filter
for file in ./02_process_out/*.1.fq.gz
do
sample=${file%%.1.fq.gz}
clone_filter -1 ${sample}.1.fq.gz -2 ${sample}.2.fq.gz -i gzfastq -o ./03_clone_filter_out
done

wait

#Rename clone_filter output to remove second number from end - Read 1
for file in ./03_clone_filter_out/*.1.1.fq.gz
do
sample=${file%%.1.1.fq.gz}
mv $file ${sample}.1.fq.gz
done

#Read 2
for file in ./03_clone_filter_out/*.2.2.fq.gz
do
sample=${file%%.2.2.fq.gz}
mv $file ${sample}.2.fq.gz
done
