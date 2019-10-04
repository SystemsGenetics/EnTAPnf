#!/usr/bin/env python3


import glob
import pandas as pd
import numpy as np
import re

filenames = sorted(glob.glob('*.tsv'))

################## IPR_mappings
file = open('IPR_mappings.txt','w')
file.close()

cols = np.arange(15)
for f in filenames:
    table = pd.read_table(f,
                          sep='\t',
                          header=None,
                          names=cols,          # to deal with some rows missing
                          usecols=[0,11,12])


    table.columns = ["Gene", "IP", "Description"]

    # remove rows with NaN
    table = table.dropna(subset= ["IP"])  # drop only the columns that has NaN in IPR column ???????????

    # remove duplicates
    table = table.drop_duplicates(keep='first')

    # removes _[[:digit]] at the end of the gene/transcript names
    table["Gene"] = table["Gene"].str.replace(r'^(.*?)_\d+', r'\1')

#     for i in table[0].index:
#         table[0][i] = re.sub(r'^(.*?)_\d+',r'\1',table[0][i])  ### rsplit option?? -

    #break

    #print(table)
    with open('IPR_mappings.txt', 'a') as txtfile:
        GO_terms = table.to_csv(txtfile,sep='\t',header = False,index=False)




############### GO_mappings
################## GO_mappings

cols = np.arange(15)
for f in filenames:
    table = pd.read_table(f,sep='\t',header=None,names=cols,usecols=[0,13])
    table.columns = ["Gene","GO"]

    # remove rows with NaN
    table = table.dropna()

    # remove duplicates
    table = table.drop_duplicates(keep='first')

    # removes _[[:digit]] at the end of the genenames
    for i in table[0].index:
         = re.sub(r'^(.*?)_\d+',r'\1',table[0][i])
        table[13][i].split('|')
    # separating the different GO terms and assigning them to the gene/transcript



    #print(table)
    #with open('IPR_mappings.txt', 'a') as txtfile:
    #    table.to_csv(txtfile,sep='\t',header = False,index=False)
