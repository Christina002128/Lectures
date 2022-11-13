#!/usr/local/bin/python3

# user input protein and taxinomic group
def pro_tx_seq():
    import os
    while True:
        protein=input("Please type in protein name:\n")
        taxonomy=input("Please type in taxonomy ID. For example, txid4890.\n")
        taxonomy=taxonomy.replace(' ','').lower() # delete spaces and change to lower case
        if taxonomy.startswith("txid"):
            print("searching proteins sequences from ncbi...") # start program
            # count the number of sequences searched from NCBI
            ecount='esearch -db Protein -query "('+protein+'[Protein Name]) AND '+taxonomy+'[Organism]" | grep "Count" > count.txt'
            os.system(ecount)
            count=open("count.txt").read()
            count=int(count.replace('<Count>','').replace('</Count>\n','').strip())
            # too many sequences, ask user 
            if count>1000:
                test=input("number of sequences more than 1000:",count,". Continue downloading? yes/no:\n")
                if test=='no':
                    print("Bye!\n")
                    quit()
            #download sequences
            fastafetch='esearch -db Protein -query "('+protein+'[Protein Name]) AND '+taxonomy+'[Organism]" | efetch -format fasta > '+taxonomy+'.fa'
            os.system(fastafetch) 
            print("protein sequences inside "+taxonomy+".fa file.")
            # create a file containing headers of all sequences
            os.system("grep '>'"+taxonomy+".fa >> "+taxonomy+"_title")
            break
        else:
            print("invalid taxonomy ID, try again.\n") # invalid, continue retyping loop
    return [protein,taxonomy]
# run function
dataset=pro_tx_seq()

# check if the file is empty
def test(dataset):
    while True:
        if os.stat(dataset[0]+'.fa').st_size==0:
            print("no proteins are caught. try again.\n")
            dataset=pro_tx_seq() # re-input the protein group
        else:
            break
    return dataset
# run function
dataset=test(dataset)

# Check how many species contained in the dataset
def species(taxonomy):
    fas=open(taxonomy+"_title").read()
    fas=fas.rstrip().replace(taxonomy+".fa:>",'').split('\n')
    i=0
    import os,sys,subprocess,pandas as pd
    dic={} # dictionary { accession : species }
    for pro in fas:
        key=pro.split(' ')[0].replace(">",'')
        value=pro.split('[')[-1].replace("]",'')
        dic[key]=value
        i+=1
    protein_numbers=len(dic.values())
    species_numbers=len(set(dic.values()))
    print("There are "+protein_numbers+"protein sequences from "+species_numbers+" kinds of species.\n")

species(dataset[1]) # run function

# ask user to continue or not
condition=input("Do you want to continue analysing? yes/no:\n")
if condition.lower()=="no":
    print("Bye!\n")
    quit()

protein=dataset[0]
taxonomy=dataset[1]
# protein sequence alignment
os.system('clustalo -i '+taxonomy+'.fa -o '+taxonomy+'.align.msf --outfmt=msf --threads=100')
# basic information about sequences in an input multiple sequence alignment
os.system('infoalign '+taxonomy+'.align.msf -out '+taxonomy+'.infoalign')
# sort by %Change from low to high, indecating conservation from high to low
os.system('sort -k10 '+taxonomy+'.infoalign >'+taxonomy+'.infoalign_sort')
# determine the level of conservation and display on the screen
os.system('cat '+taxonomy+'.infoalign_sort')
# plot the the level of conservation
plot="plotcon -sformat msf "+taxonomy+".align.msf -winsize 20 -gsubtitle '"+protein+" sequences in "+taxonomy+"' -graph ps -goutfile "+taxonomy
os.system(plot)
# display the plot
taxonomy.ps
