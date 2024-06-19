##########
#
# Author: F. MASSIP
# date: 12/14/2021
#

configfile:"config.yml"


def get_countries(s,mem):
	count_list=[]
	if mem == "t" :
		myfileHM = config["DATA_DIR"]+"/list_samples/"+s+".txt_unique"
		with open(myfileHM) as fHM:
			for i in fHM.readlines():
				count_list.append(i.rstrip())
	return count_list

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
		expand("{RES_DIR}/pairs/{COMBOS2}.tsv",RES_DIR=config["RES_DIR"],COMBOS2=COMBOS2),
		expand("{RES_DIR}/{COMBOS2}/{COMBOS2}.tsv",RES_DIR=config["RES_DIR"],COMBOS2 = COMBOS2),
		expand("{RES_DIR}/{COMBOS2}/{COMBOS2}.h5",RES_DIR=config["RES_DIR"],COMBOS2 = COMBOS2)

rule Compare:
	input:
		fasta1=expand("{DATA_DIR}/fasta/{{SP1}}.fna",DATA_DIR=config["DATA_DIR"]),
		fasta2=expand("{DATA_DIR}/fasta/{{SP2}}.fna",DATA_DIR=config["DATA_DIR"])

	output:
		out = expand("{RES_DIR}/pairs/{{SP2}}_{{SP1}}.tsv",RES_DIR=config["RES_DIR"]),

	conda:
		config["CONDA_FILE2"]
	shell:
		"""mkdir -p {config[RES_DIR]}/{wildcards.SP2}_{wildcards.SP1}/
		cd {config[RES_DIR]}/{wildcards.SP2}_{wildcards.SP1}/
		{config[CODE_DIR]}/bfmem/bfmem -s b -t 1 -l 100  -o {output.out} -r {input.fasta1} -q {input.fasta2}"""


rule Intersect:
	input:
		 expand("{RES_DIR}/pairs/{{SP2}}_{{SP1}}.tsv",RES_DIR=config["RES_DIR"]),
	output:
		h5 = expand("{RES_DIR}/{{SP2}}_{{SP1}}/{{SP2}}_{{SP1}}.h5",RES_DIR=config["RES_DIR"]),
		tsv = expand("{RES_DIR}/{{SP2}}_{{SP1}}/{{SP2}}_{{SP1}}.tsv",RES_DIR=config["RES_DIR"])

	conda:
		config["CONDA_FILE2"]
	shell:
		"""Rscript {config[CODE_DIR]}/Intersect2.R {wildcards.SP1} {wildcards.SP2} {output.h5} """
