#!/bin/bash
#SBATCH -J relernn
#SBATCH -o relernn."%j".out
#SBATCH -e relernn."%j".err
#SBATCH --mail-user marine.salson@etu.umontpellier.fr
#SBATCH --mail-type=ALL
#SBATCH -p gpu
#SBATCH -A gpu_group
#SBATCH -w node26
#SBATCH -c 10

ReLERNN_SIMULATE -v spatial_dataset.vcf -g all_chr.bed -d output/

ReLERNN_TRAIN -d output/

ReLERNN_PREDICT -v spatial_dataset.vcf -d output/

ReLERNN_BSCORRECT -d output/ -t 10


