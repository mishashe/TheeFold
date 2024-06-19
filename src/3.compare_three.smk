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

def get_trios(Genus_list):
	all_trios = []
	all_trios2 = []
	n = len(Genus_list)
	for i in range(0,n-2):
		for j in range(i+1,n-1):
			for k in range(j+1,n):
				trio = Genus_list[i] + " " + Genus_list[j] + " " + Genus_list[k]
				trio2 = Genus_list[i] + "_" + Genus_list[j] + "_" + Genus_list[k]
				all_trios.append(trio)
				all_trios2.append(trio2)
	return (all_trios,all_trios2)


(TRIOS,TRIOS2) = get_trios(config["GENUS"])

rule all:	
	input:
		expand("{RES_DIR}/{TRIOS2}/{TRIOS2}.tsv",RES_DIR=config["RES_DIR"],TRIOS2 = TRIOS2)


rule Intersect:
	input:
		 expand("{RES_DIR}/pairs/{{SP3}}_{{SP2}}.tsv",RES_DIR=config["RES_DIR"]),
		 expand("{RES_DIR}/pairs/{{SP3}}_{{SP1}}.tsv",RES_DIR=config["RES_DIR"]),
		 expand("{RES_DIR}/pairs/{{SP2}}_{{SP1}}.tsv",RES_DIR=config["RES_DIR"])

	output:
		tsv = expand("{RES_DIR}/{{SP3}}_{{SP2}}_{{SP1}}/{{SP3}}_{{SP2}}_{{SP1}}.tsv",RES_DIR=config["RES_DIR"])

	conda:
		config["CONDA_FILE2"]
	shell:
		"""Rscript {config[CODE_DIR]}/Intersect3.R {wildcards.SP1} {wildcards.SP2} {wildcards.SP3} {output.tsv} """
