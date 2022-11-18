#!/usr/local/bin/python3
import os
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
            # count the number of sequences searched from NCBI
            ecount='esearch -db protein -query "('+protein+')[Protein Name] AND '+taxonomy+'[Organism]" | grep "Count" > count.txt'
            os.system(ecount)
            count=open("count.txt").read()
            count=int(count.replace('<Count>','').replace('</Count>\n','').strip()) 
            if count>1000: # too many sequences, ask user download or not
                print("So many sequences("+str(count)+")! Continue downloading? yes/no:\n")
                test=input()
                if test=='no':
                    print("Bye!\n")
                    quit()
            if count==0: # check if the retrieve result is 0
                print("0 result searched. Try again.\n")
                continue  # input again
            else:
                print(str(count)+" sequences downloading...\n")
            #download sequences
            fastafetch='esearch -db protein -query "('+protein+')[Protein Name] AND '+taxonomy+'[Organism]" | efetch -format fasta > '+taxonomy+'.fa'
            os.system(fastafetch) 
            print("protein sequences inside "+taxonomy+".fa file")
            os.system("grep '>' "+taxonomy+".fa > "+taxonomy+"_title") # create a file containing headers of all sequences
            print("protein headers inside "+taxonomy+"_title file")
            break
        else:
            print("invalid taxonomy ID, try again.\n") # invalid, continue retyping loop
    return [protein,taxonomy]


# Check how many species contained in the dataset
def species(taxonomy):
    fas=open(taxonomy+"_title").read()
    fas=fas.rstrip().replace(taxonomy+".fa:>",'').split('\n')
    i=0
    dic={} # create dictionary { accession : species }
    for pro in fas:
        key=pro.split(' ')[0].replace(">",'')
        value=pro.split('[')[-1].replace("]",'')
        dic[key]=value
        i+=1
    # count the number of proteins and species
    protein_numbers=len(dic.values())
    species_numbers=len(set(dic.values()))
    print("It contains "+str(protein_numbers)+" sequences and "+str(species_numbers)+" kinds of species.")
    # write result to file
    with open(taxonomy+'pro_spe.txt','w') as output:
        output.write('protein\tspecies')
        for key,value in dic.items():
            output.write('%s\t%s\n' % (key, value))
    print("output")
    return dic

# count sequence length and output fasta file for each sequence
def seq_len_out(taxonomy):
    fas=open(taxonomy+'.fa').read().strip()
    full_seq=fas.split('>')[1:] # delete first empty string
    filename=[]
    seq=[]
    pro=[]
    seqlen=[]
    for i in list(range(0,len(full_seq),1)): 
        # count sequence length for each sequence
        devide=seq[i].find("\n")
        pro[i]=full_seq[i][0:devide]
        seq[i]=full_seq[i][(devide+1):].replace('\n','')
        seqlen[i]=len(seq[i])
        # output fasta file for each sequence
        filename.append(full_seq[i].split(' ')[0])  # set sequence name as file name
        full_fa[i]='>'+full_seq[i]   # add separator('>') back
        with open(filename[i]+'.seq.fa',"w") as outfile:
            outfile.write(full_fa[i]) # write in single sequence
    import pandas as pd
    df={} # data frame contain the protein name, sequence, sequence length
    df["protein"]=pro
    df['sequence']=seq
    df['length']=seqlen 
    return df   


# determine, and plot, the level of protein sequence conservation
def conservation(dataset):
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
    display="plotcon -sformat msf txid4890.align.msf -winsize 25 -gsubtitle '"+protein+" sequences in "+taxonomy+"' -graph xterm"
    os.system(display)

