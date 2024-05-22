##########
#
# Author: F. MASSIP
# date: 24/02/2023
#
# Pipeline to download genomes of different bacteria,
# and split multifasta to single fasta file

configfile:"config_dl.yml"

#SPECIES=[config["SPECIES1"],config["SPECIES2"]]

rule all:
	input:
#		expand("{DATA_DIR}/list_samples/{SPECIES}.txt",DATA_DIR=config["DATA_DIR"],SPECIES=config["SPECIES"]),
		expand("{DATA_DIR}/list_samples/{GENUS}_csv",DATA_DIR=config["DATA_DIR"],GENUS=config["GENUS"])


rule DownloadSP:
	input:
		meta=expand("{DATA_DIR}/bac120_metadata_r220.tsv",DATA_DIR=config["DATA_DIR"])
	output:
		listFile=expand("{DATA_DIR}/list_samples/{{GENUS}}_csv",DATA_DIR=config["DATA_DIR"]),
		fasta=directory(config["DATA_DIR"]+"/fasta/{GENUS}/")
	conda:
		config["CONDA_FILE"]
	shell:
		"""python3 {config[CODE_DIR]}/select_genomes.py -m {input.meta} \
                   -o {config[DATA_DIR]}/list_samples/{wildcards.GENUS}_to_dl.txt \
                   -l genus -t {wildcards.GENUS}

                   python3 {config[CODE_DIR]}/download_api.py {config[DATA_DIR]}/list_samples/{wildcards.GENUS}_to_dl.txt \
                   {config[DATA_DIR]}/new_taxdump/rankedlineage.dmp {config[DATA_DIR]}/bac120_taxonomy_r220.tsv.gz \
                   {config[DATA_DIR]}/fasta/{wildcards.GENUS}  {output.listFile}"""



