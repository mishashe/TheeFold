##########
#
# Author: F. MASSIP
# date: 24/02/2023
#
# Pipeline to download genomes of different bacteria,
# and split multifasta to single fasta file

configfile:"config_mobilome.yml"

#SPECIES=[config["SPECIES1"],config["SPECIES2"]]


def get_combos(Genus_list):
	all_combos = []
	all_combos2 = []
	G = Genus_list
	n = len(G)
	for i1 in range(0,n-1):
		for i2 in range(i1+1,n):
			combo = G[i1] + " " + G[i2] 
			combo2 = G[i1] + "_" + G[i2] 

			all_combos.append(combo)
			all_combos2.append(combo2)
	return (all_combos,all_combos2)


(COMBOS,COMBOS2) = get_combos(config["GENUS"])

print(COMBOS)


rule all:
	input:
		#expand("{DATA_DIR}/list_samples/{GENUS}_csv",DATA_DIR=config["DATA_DIR"],GENUS=config["GENUS"]),
		#expand("{DATA_DIR}/fasta/{GENUS}.fna",DATA_DIR=config["DATA_DIR"],GENUS=config["GENUS"])
		expand("{DATA_DIR}/mobilome/{GENUS}_seq.fa",DATA_DIR=config["DATA_DIR"],GENUS=config["GENUS"]),
		expand("{DATA_DIR}/mobilome/{GENUS}_listFiles.txt",DATA_DIR=config["DATA_DIR"],GENUS=config["GENUS"]),
		expand("{DATA_DIR}/mobilome/alignments/pairs/{COMBOS2}.tsv",DATA_DIR=config["DATA_DIR"],COMBOS2=COMBOS2)
		expand("{DATA_DIR}/mobilome/Matches/pairs/{COMBOS2}.tsv",DATA_DIR = config["DATA_DIR"],COMBOS2=COMBOS2)

rule Extract:
	input:
		listFile=expand("{DATA_DIR}/external/list_samples/{{GENUS}}_csv",DATA_DIR=config["DATA_DIR"]),
		db = expand("{DATA_DIR}/external/database/nuc_db.fa",DATA_DIR=config["DATA_DIR"]),
		pdb = expand("{DATA_DIR}/external/database/prot_db.fa",DATA_DIR=config["DATA_DIR"])
	output:
		fasta=expand("{DATA_DIR}/mobilome/{{GENUS}}_seq.fa",DATA_DIR=config["DATA_DIR"]),
		listFiles = expand("{DATA_DIR}/mobilome/{{GENUS}}_listFiles.txt",DATA_DIR=config["DATA_DIR"])
	conda:
		config["CONDA_FILE2"]
	shell:
		"""python3 {config[CODE_DIR]}/extract_mobilome.py -c {config[DATA_DIR]}/external/fasta/{wildcards.GENUS}\
		 -n 20 \
		-o {config[DATA_DIR]}/mobilome/{wildcards.GENUS} \
		-d {config[DATA_DIR]}/external/database/nuc_db.fa \
		-p {config[DATA_DIR]}/external/database/prot_db.fa"""

rule align:
	input:
		fasta1=expand("{DATA_DIR}/mobilome/{{SP1}}_seq.fa",DATA_DIR=config["DATA_DIR"]),
		fasta2=expand("{DATA_DIR}/mobilome/{{SP2}}_seq.fa",DATA_DIR=config["DATA_DIR"])
	output:
		out = expand("{DATA_DIR}/mobilome/alignments/pairs/{{SP2}}_{{SP1}}.tsv",DATA_DIR = config["DATA_DIR"])
	conda:
		config["CONDA_FILE2"]
	shell:
		"""python3 {config[CODE_DIR]}/align_and_extract_coord_matches.py --st1 {input.fasta1} \
		--st2 {input.fasta2} \
		-o {output.out}"""

rule getMatchesNoAlign:
	input:
		fasta1=expand("{DATA_DIR}/mobilome/{{SP1}}_seq.fa",DATA_DIR=config["DATA_DIR"]),
		fasta2=expand("{DATA_DIR}/mobilome/{{SP2}}_seq.fa",DATA_DIR=config["DATA_DIR"])
	output:
		out = expand("{DATA_DIR}/mobilome/Matches/pairs/{{SP2}}_{{SP1}}.tsv",DATA_DIR = config["DATA_DIR"])
	conda:
		config["CONDA_FILE2"]
	shell:
		"""mkdir -p {config[DATA_DIR]}/mobilome//Matches/{wildcards.SP2}_{wildcards.SP1}/
 		cd {config[RES_DIR]}/mobilome//Matches/{wildcards.SP2}_{wildcards.SP1}/
		{config[CODE_DIR]}/../bfmem/bfmem -s b -t 1 -l 20 -k 20  -o {output.out} -r {input.fasta1} -q {input.fasta2}"""


#rule Intersect:
#	input:
#		 expand("{DATA_DIR}/mobilome/alignemnts/pairs/{{SP2}}_{{SP1}}.tsv",DATA_DIR=config["DATA_DIR"]),
#	output:
#		h5 = expand("{DATA_DIR}/mobilome/results/{{SP2}}_{{SP1}}/{{SP2}}_{{SP1}}.h5",RES_DIR=config["DATA_DIR"]),
#		tsv = expand("{DATA_DIR}/mobilome/results/{{SP2}}_{{SP1}}/{{SP2}}_{{SP1}}.tsv",RES_DIR=config["DATA_DIR"])
#
#	conda:
#		config["CONDA_FILE2"]
#	shell:
#		"""Rscript {config[CODE_DIR]}/Intersect2.R {wildcards.SP1} {wildcards.SP2} {output.h5} """
#
#
#
#


