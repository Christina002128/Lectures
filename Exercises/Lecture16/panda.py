# exercise 16 #
#!/usr/local/bin/python3
import os, sys, re, numpy as np
import pandas as pd
df = pd.read_csv('eukaryotes.tsv', sep="\t", na_values=['-'])
df.index=df.apply(lambda x : "{} ({})".format(x['#Organism/Name'], x['BioSample Accession']), axis=1)
fun_100=df[(df['Size (Mb)']>100) & (df['Group']=='Fungi')].loc[:,['Group']].value_counts()[0]
fun_100
# 178
df[df['Group'].isin(['Plants','Animals','Fungi','Protists'])]['Group'].value_counts()

df[df.apply(lambda x : x['#Organism/Name'].startswith('Heliconius'),axis=1)][['#Organism/Name', 'Scaffolds']]

(df['Group'].str.lower()=='plants').value_counts() # number of plants

df[df['Group'].str.lower()=='plants'].sort_values(['Genes'],ascending=[False])[['Group','Center','Genes']].head()
# str.lower in re
df['proteins per gene']= df['Proteins'] / df['Genes']
df[df['proteins per gene']>1.1][['#Organism/Name','proteins per gene']].sort_values(['proteins per gene'],ascending=False)
len(df[df['proteins per gene']>1.1]) # number of lines 

