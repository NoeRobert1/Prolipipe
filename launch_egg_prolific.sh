#!/bin/bash
#SBATCH --job-name=egg_prolific
#SBATCH --chdir=/scratch/norobert/prolific_project/
#SBATCH -p genouest
#SBATCH -o log_%x.log            # output file
#SBATCH -e log_%x.log            # error file
#SBATCH --cpus-per-task 8

. /local/env/envconda.sh 
#conda activate ~/miniconda3/envs/prolific
conda activate /home/genouest/dyliss/ytirlet/my_env

# using eggnog-mapper
source /local/env/enveggnog-mapper-2.1.9.sh 
/home/genouest/dyliss/ytirlet/my_env/bin/python prolific/pipeline.py --input prolific/toy_example/genomes/ \
--padmet_ref /scratch/norobert/metacyc_26.0.padmet \
--output /scratch/norobert/ \
--ptsc /scratch/norobert/ \
--ptsi mpwt_26.sif \
--tax prolific/toy_example/all_taxons.tsv \
--pwy prolific/toy_example/pathways/ \
--annot eggnog \
--egg_path /db/eggnog/5.0.2/ -k -v 
