#!/usr/local/bin/python3

# user choose data, and ask user to do what
protein=input("protein name:\n")
organism=input("taxonomy name:\n")
accfetch='esearch -db Protein -query "('+protein+'[Protein Name]) AND '+organism+'[Organism]" | efetch -format acc > proteins.acc'
import os
os.system(accfetch)
fastafetch='esearch -db Protein -query "('+protein+'[Protein Name]) AND '+organism+'[Organism]" | efetch -format acc > proteins.fasta'
os.system(fastafetch)
