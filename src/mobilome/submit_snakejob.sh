#!/bin/bash
#SBATCH -J ControlJob
#SBATCH -t 20-00:00:00
#SBATCH -N 1
#SBATCH --output /cluster/CBIO/data1/fmassip/HGT/ProjectMisha/HGTnew/multi_comparisons/ThreeFold/src/logs/%x-%j.out
#SBATCH --error /cluster/CBIO/data1/fmassip/HGT/ProjectMisha/HGTnew/multi_comparisons/ThreeFold/src/logs/%x-%j.err
eval "$(conda shell.bash hook)"

# activate conda 
conda activate snakemake

# make things fail on errors
set -o nounset
set -o errexit
set -x

export LOGDIR=/cluster/CBIO/data1/fmassip/HGT/ProjectMisha/HGTnew/multi_comparisons/ThreeFold/src/logs/${SLURM_JOB_NAME}-${SLURM_JOB_ID}
### run your commands here!

mkdir -p $LOGDIR

snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/mobilome/1.construct_mobilomes.smk -n --unlock

snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/mobilome/1.construct_mobilomes.smk \
           --use-conda \
	   --cluster-config ~/HGTnew/multi_comparisons/ThreeFold/src/config_sge.yml \
	   --cluster "sbatch -N 1 -c 1 -J Mobilome  -o $LOGDIR/%j.log -t {cluster.time} --mem {cluster.mem}" \
	   --jobs 20 \
           --configfile ~/HGTnew/multi_comparisons/ThreeFold/src/mobilome/config_mobilome.yml \
	   --rerun-incomplete 

#snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/2.compare_pairwise.smk -n --unlock
#
#snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/2.compare_pairwise.smk \
#           --use-conda \
#	   --cluster-config config_sge.yml \
#	   --cluster "sbatch -N 1 -c 1 -J PntLum2  -o $LOGDIR/%j.log -t {cluster.time} --mem {cluster.mem}" \
#	   --jobs 30 \
#           --configfile config.yml \
#	   --rerun-incomplete 
#
#snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/3.compare_three.smk -n --unlock
#
#snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/3.compare_three.smk \
#           --use-conda \
#	   --cluster-config config_sge.yml \
#	   --cluster "sbatch -N 1 -c 1 -J PntLum3  -o $LOGDIR/%j.log -t {cluster.time} --mem {cluster.mem}" \
#	   --jobs 30 \
#           --configfile config.yml \
#	   --rerun-incomplete 
#
#
#snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/4.compare.smk -n --unlock
#
#snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/4.compare.smk \
#           --use-conda \
#	   --cluster-config config_sge.yml \
#	   --cluster "sbatch -N 1 -c 1 -J PntLum4  -o $LOGDIR/%j.log -t {cluster.time} --mem {cluster.mem}" \
#	   --jobs 30 \
#           --configfile config.yml \
#	   --rerun-incomplete 
#
#snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/5.compare.smk -n --unlock
#
#snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/5.compare.smk \
#           --use-conda \
#	   --cluster-config config_sge.yml \
#	   --cluster "sbatch -N 1 -c 1 -J PntLum5  -o $LOGDIR/%j.log -t {cluster.time} --mem {cluster.mem}" \
#	   --jobs 30 \
#           --configfile config.yml \
#	   --rerun-incomplete 
#
#snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/6.compare.smk -n --unlock
#
#snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/6.compare.smk \
#           --use-conda \
#	   --cluster-config config_sge.yml \
#	   --cluster "sbatch -N 1 -c 1 -J PntLum6  -o $LOGDIR/%j.log -t {cluster.time} --mem {cluster.mem}" \
#	   --jobs 20 \
#           --configfile config.yml \
#	   --rerun-incomplete 
#
#snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/7.compare.smk -n --unlock
#
#snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/7.compare.smk \
#           --use-conda \
#	   --cluster-config config_sge.yml \
#	   --cluster "sbatch -N 1 -c 1 -J PntLum7 -o $LOGDIR/%j.log -t {cluster.time} --mem {cluster.mem}" \
#	   --jobs 10 \
#           --configfile config.yml \
#	   --rerun-incomplete 
#
#
#
#snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/8.compare.smk -n --unlock
#
#snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/8.compare.smk \
#           --use-conda \
#	   --cluster-config config_sge.yml \
#	   --cluster "sbatch -N 1 -c 1 -J PntLum8  -o $LOGDIR/%j.log -t {cluster.time} --mem {cluster.mem}" \
#	   --jobs 2 \
#           --configfile config.yml \
#	   --rerun-incomplete 
#


