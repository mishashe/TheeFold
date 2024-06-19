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

snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/1.download.smk -n --unlock

snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/1.download.smk \
           --use-conda \
	   --cluster-config config_sge.yml \
	   --cluster "sbatch -N 1 -c 1 -J PntLum1  -o $LOGDIR/%j.log -t {cluster.time} --mem {cluster.mem}" \
	   --jobs 30 \
           --configfile config_dl.yml \
	   --rerun-incomplete 

snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/2.compare_pairwise.smk -n --unlock

snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/2.compare_pairwise.smk \
           --use-conda \
	   --cluster-config config_sge.yml \
	   --cluster "sbatch -N 1 -c 1 -J PntLum2  -o $LOGDIR/%j.log -t {cluster.time} --mem {cluster.mem}" \
	   --jobs 30 \
           --configfile config.yml \
	   --rerun-incomplete 

snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/3.compare_three.smk -n --unlock

snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/3.compare_three.smk \
           --use-conda \
	   --cluster-config config_sge.yml \
	   --cluster "sbatch -N 1 -c 1 -J PntLum3  -o $LOGDIR/%j.log -t {cluster.time} --mem {cluster.mem}" \
	   --jobs 30 \
           --configfile config.yml \
	   --rerun-incomplete 


snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/4.compare.smk -n --unlock

snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/4.compare.smk \
           --use-conda \
	   --cluster-config config_sge.yml \
	   --cluster "sbatch -N 1 -c 1 -J PntLum4  -o $LOGDIR/%j.log -t {cluster.time} --mem {cluster.mem}" \
	   --jobs 30 \
           --configfile config.yml \
	   --rerun-incomplete 

snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/5.compare.smk -n --unlock

snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/5.compare.smk \
           --use-conda \
	   --cluster-config config_sge.yml \
	   --cluster "sbatch -N 1 -c 1 -J PntLum5  -o $LOGDIR/%j.log -t {cluster.time} --mem {cluster.mem}" \
	   --jobs 30 \
           --configfile config.yml \
	   --rerun-incomplete 

snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/6.compare.smk -n --unlock

snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/6.compare.smk \
           --use-conda \
	   --cluster-config config_sge.yml \
	   --cluster "sbatch -N 1 -c 1 -J PntLum6  -o $LOGDIR/%j.log -t {cluster.time} --mem {cluster.mem}" \
	   --jobs 20 \
           --configfile config.yml \
	   --rerun-incomplete 

snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/7.compare.smk -n --unlock

snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/7.compare.smk \
           --use-conda \
	   --cluster-config config_sge.yml \
	   --cluster "sbatch -N 1 -c 1 -J PntLum7 -o $LOGDIR/%j.log -t {cluster.time} --mem {cluster.mem}" \
	   --jobs 10 \
           --configfile config.yml \
	   --rerun-incomplete 



snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/8.compare.smk -n --unlock

snakemake -s ~/HGTnew/multi_comparisons/ThreeFold/src/8.compare.smk \
           --use-conda \
	   --cluster-config config_sge.yml \
	   --cluster "sbatch -N 1 -c 1 -J PntLum8  -o $LOGDIR/%j.log -t {cluster.time} --mem {cluster.mem}" \
	   --jobs 2 \
           --configfile config.yml \
	   --rerun-incomplete 


#for i in `head -n 2  Species_list`
#do	
#	for j in `cat Species_list`
#	do
#		if [ $i != $j ]
#			then
#			my_list=$(echo $i $j | xargs -n1 | sort | xargs)
#			cat config_min.yml >config.yml
#			echo $my_list |sed -r 's/(.*) (.*)/SPECIES1 : \1\nSPECIES2 : \2/' >>config.yml
#			echo $my_list |sed -r 's/(.*) (.*)/SPECIES1 : \1\nSPECIES2 : \2/' >>test
#
#
#			snakemake -s ~/HGTnew/dating_long_distance/scripts/2b.lastZ.smk -n --unlock
#				
#			snakemake -s ~/HGTnew/dating_long_distance/scripts/2b.lastZ.smk \
#		           --use-conda \
#			   --cluster-config config_sge.yml \
#			   --cluster "sbatch -N 1 -c 1 -J lastZ -o $LOGDIR/%j.log -t {cluster.time} --mem {cluster.mem}" \
#			   --jobs 92 \
#			   --rerun-incomplete \
#			   --configfile config.yml \
#			   --latency-wait 60
##			   --resources cp_cores=10 \
#	
#	
#		fi
#	done
#done
#
#snakemake -s ~/HGTnew/dating_long_distance/scripts/3.measure_divergence.smk --unlock -n \
# --configfile config_dl.yml 
#
#
#snakemake -s ~/HGTnew/dating_long_distance/scripts/3.measure_divergence.smk \
#   --use-conda \
#   --cluster-config config_sge.yml \
#   --cluster "sbatch -N 1 -c 1 -J Diverg  -o $LOGDIR/%j.log -t {cluster.time} --mem {cluster.mem}" \
#   --jobs 92 \
#   --rerun-incomplete \
#   --configfile config_dl.yml \
#   --latency-wait 60
##			   --resources cp_cores=10 \
#
