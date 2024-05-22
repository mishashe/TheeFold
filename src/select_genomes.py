# -*-coding:utf8 -*

import os
import sys
#import numpy as np
import random
import re
from collections import defaultdict
from optparse import OptionParser
import pandas as pd
################################



parser = OptionParser()
parser.add_option("-m","--metadata_file",dest="metadata_file",action="store",
help="file containing metadata about gtdb data")
parser.add_option("-o","--output_file",dest="output_file",action="store",
help="output file name")
parser.add_option("-l", "--taxa_level",dest="taxa_level",action="store",
help="level of the taxa we want to dl. should be one of [domain, phylum, class, order, family, genus, species]")
parser.add_option("-t", "--taxa_to_dl",dest="taxa_to_dl",action="store",
help="name of the taxa we want to dl")


(options, args) = parser.parse_args()
patt=re.compile("start=([0-9]+)")

metadata_file = options.metadata_file
taxa_level = options.taxa_level
taxa_to_dl = options.taxa_to_dl

output_file=options.output_file

#ori=defaultdict(list)


def clean_gtdb(gtdb_taxo):
    """cleans the gtdb 120 taxonomy"""
    taxonomic_levels = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    taxo_values = gtdb_taxo.split(";")
    gtdb_df = {}
    
    for key in taxonomic_levels:
        for value in taxo_values:
            gtdb_df[key] = value[3:]
            taxo_values.remove(value)
            break

#    gtdb = gtdb_taxo.split(";")
#    gtdb_df[taxonomic_levels] = gtdb_taxo.split(';')
#    for col in range(0:len(gtdb_df)):
#        gtdb_df[col] = gtdb_df[col].apply(lambda x : x[3:])
##    gtdb_df["assembly_accession"] = gtdb_df["assembly_accession"].apply(lambda x : x[3:])
##    gtdb_df = gtdb_df.drop(columns=['taxonomy'])
    return gtdb_df

f1 = open(output_file,"w")
f1.write("assembly_accession\n")

#count = 0
with open (metadata_file,'r') as meta: 
    meta.readline()
    for line in meta.readlines():
        line=line.rstrip()
        row=line.split('\t')
#        print(row[5],"\n",row[6],"\n",row[19])
#        gtdb_tax = pd.read_csv(row[19], sep="\t", names=["assembly_accession", "taxonomy"])
        gtdb = clean_gtdb(row[19])
#       print(gtdb)
        if float(row[5]) > 90 and float(row[6]) < 5:  #keep only genomes wich checkmcompletness >90% and contamination <5%
            if row[58] != "derived from metagenome" : #remove genomes extracted from metagenome analysis
                if gtdb [ taxa_level ] == taxa_to_dl:
                    f1.write(row[0][3:]+"\n")
#                    count = count+1
#                    if count >10:
#                        break

f1.close()

