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
	for i1 in range(0,n-7):
		for i2 in range(i1+1,n-6):
			for i3 in range(i2+1,n-5):
				for i4 in range(i3+1,n-4):
					for i5 in range(i4+1,n-3):
						for i6 in range(i5+1,n-2):
							for i7 in range(i6+1,n-1):
								for i8 in range(i7+1,n):

									combo = G[i1] + " " + G[i2] + " " + G[i3] + " " + G[i4] + " " + G[i5] + " " + G[i6] + " " + G[i7] + " " + G[i8]
									combo2 = G[i1] + "_" + G[i2] + "_" + G[i3] + "_" + G[i4] + "_" + G[i5] + "_" + G[i6] + "_" + G[i7] + "_" + G[i8]

									all_combos.append(combo)
									all_combos2.append(combo2)
	return (all_combos,all_combos2)


(COMBOS,COMBOS2) = get_combos(config["GENUS"])
print (COMBOS)
print (COMBOS2)

rule all:	
	input:
		expand("{RES_DIR}/{COMBOS2}/{COMBOS2}.tsv",RES_DIR=config["RES_DIR"],COMBOS2 = COMBOS2)

rule Intersect:
	input:
		expand("{RES_DIR}/pairs/{{SP8}}_{{SP4}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP8}}_{{SP3}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP8}}_{{SP2}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP8}}_{{SP1}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP8}}_{{SP5}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP8}}_{{SP6}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP8}}_{{SP7}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP7}}_{{SP4}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP7}}_{{SP3}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP7}}_{{SP2}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP7}}_{{SP1}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP7}}_{{SP5}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP7}}_{{SP6}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP6}}_{{SP4}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP6}}_{{SP3}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP6}}_{{SP2}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP6}}_{{SP1}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP6}}_{{SP5}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP5}}_{{SP4}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP5}}_{{SP3}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP5}}_{{SP2}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP5}}_{{SP1}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP4}}_{{SP3}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP4}}_{{SP2}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP4}}_{{SP1}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP3}}_{{SP2}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP3}}_{{SP1}}.tsv",RES_DIR=config["RES_DIR"]),
		expand("{RES_DIR}/pairs/{{SP2}}_{{SP1}}.tsv",RES_DIR=config["RES_DIR"])

	output:
		tsv = expand("{RES_DIR}/{{SP8}}_{{SP7}}_{{SP6}}_{{SP5}}_{{SP4}}_{{SP3}}_{{SP2}}_{{SP1}}/{{SP8}}_{{SP7}}_{{SP6}}_{{SP5}}_{{SP4}}_{{SP3}}_{{SP2}}_{{SP1}}.tsv",RES_DIR=config["RES_DIR"])

	conda:
		config["CONDA_FILE2"]
	shell:
		"""Rscript {config[CODE_DIR]}/Intersect8.R {wildcards.SP1} {wildcards.SP2} {wildcards.SP3} {wildcards.SP4} {wildcards.SP5} {wildcards.SP6} {wildcards.SP7} {wildcards.SP8} {output.tsv} """
