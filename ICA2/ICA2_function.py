#!/usr/local/bin/python3
import os
import pandas as pd
# user input protein name and taxinomic group ID
def pro_tx_seq():
    while True: # a loop till correct input
        protein=input("Please type in protein name:\n")
        taxonomy=input("Please type in taxonomy ID. For example, txid4890\n")
        taxonomy=taxonomy.replace(' ','').lower() # delete spaces and change to lower case
        # test if taxonomy input is valid
        if taxonomy.startswith("txid"):
            # test if id is integer
            id=taxonomy.replace("txid",'') 
            try:
                int(id) 
            except:
                print("invalid taxonomy ID, try again.\n")
                continue  # input again
            print("searching proteins sequences from ncbi...") # start program
            #partial or not partial search
            par=input("would you like to contain partial proteins or not? yes/no:\n")
            par=par.replace(' ','').lower()
            if par=="no":
                print("no partial searching...\n")
                ecount='esearch -db protein -query "('+protein+')[Protein Name] AND '+taxonomy+'[Organism] NOT PARTIAL" | grep "Count" > count.txt'
                par_query='NOT PARTIAL'
            else:
                print("full searching...\n")
                ecount='esearch -db protein -query "('+protein+')[Protein Name] AND '+taxonomy+'[Organism]" | grep "Count" > count.txt'
                os.system(ecount)
                par_query=''
            # count the number of sequences searched from NCBI
            count=open("count.txt").read()
            count=int(count.replace('<Count>','').replace('</Count>\n','').strip()) 
            if count>1000: # too many sequences, ask user download or not
                print("Too many sequences("+str(count)+")!  Continue downloading? yes/no:\n")
                test=input()
                if test.lower().replace(' ', '')=='no':
                    print("Bye!\n")
                    quit()
            if count==0: # check if the retrieve result is 0
                print("0 result searched. Try again.\n")
                continue  # input again
            else:
                print(str(count)+" sequences downloading...\n")
            #download sequences
            fastafetch='esearch -db protein -query "('+protein+')[Protein Name] AND '+taxonomy+'[Organism]'+par_query+'" | efetch -format fasta > '+taxonomy+'.fa'
            os.system(fastafetch) 
            print("protein sequences is inside file: "+taxonomy+".fa")
            os.system("grep '>' "+taxonomy+".fa > "+taxonomy+"_title") # create a file containing headers of all sequences
            print("protein headers is inside file: "+taxonomy+"_title")
            break
        else:
            print("invalid taxonomy ID, try again.\n") # invalid, continue retyping loop
    return [protein,taxonomy]


# Check how many species contained in the dataset
def species(taxonomy):
    fas=open(taxonomy+"_title").read()
    fas=fas.rstrip().replace(taxonomy+".fa:>",'').split('\n')
    dic={} # create dictionary { accession : species }
    for pro in fas:
        if pro.find("["): # 
            key=pro.split('\t')[0].replace(">",'')
            value=pro.split('[')[-1].replace("]",'')
            dic[key]=value
    # count the number of proteins and species
    protein_numbers=len(dic.values())
    species_numbers=len(set(dic.values()))
    percent=round(protein_numbers/species_numbers,2)
    print("It contains "+str(protein_numbers)+" sequences and "+str(species_numbers)+" kinds of species.")
    print("About "+str(percent)+" proteins per species.\n")
    # write result to file
    with open(taxonomy+'proteins_species.txt','w') as output:
        output.write('protein\tspecies\n')
        for key,value in dic.items():
            output.write('%s\t%s\n' % (key, value))
    print("output proteins and corresponding species file: "+taxonomy+'proteins_species.txt')
    return dic

import os
import pandas as pd
import matplotlib.pyplot as plt
# count sequence length
def seq_len(taxonomy):
    fas=open(taxonomy+'.fa').read().strip()
    full_seq=fas.split('>')[1:] # delete first empty string
    seq=[]
    raw_seq=[]
    pro=[]
    seqlen=[]
    pro_number=[]
    for i in list(range(0,len(full_seq),1)): 
        # count sequence length for each sequence
        devide=full_seq[i].find("\n") # first index of '\n' to divide header and sequence
        pro.append(full_seq[i][0:devide])
        pro_number.append(full_seq[i].split(' ')[0])
        raw_seq.append(full_seq[i][(devide+1):])
        seq.append(full_seq[i][(devide+1):].replace('\n',''))
        seqlen.append(len(seq[i]))
    df={} # data frame contain the protein name, sequence, sequence length
    df=pd.DataFrame({"protein":pro,"sequence":seq,"length":seqlen,'raw_seq':raw_seq,'pro_number':pro_number})
    df=df.sort_values('length',ascending=False) # sort by length
    return df
