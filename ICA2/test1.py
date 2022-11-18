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

# choose a appropriate protein subset for alignment by considering the sequence length
def seq_len_subset(df,taxonomy):
    print("Summary of sequences length:\n")
    print(df.describe()) # print out calculations of length
    # display histogram of sequence length
    print("histogram of sequences length.\n")
    y=list(df['length'])
    plt.hist(y, 50)
    plt.grid(True)
    plt.xlabel('Sequence length')
    plt.ylabel('Number of proteins')
    plt.title('Histogram of sequences length in'+taxonomy+" group")
    plt.savefig(taxonomy+".seqLenHis.png")
    print("the plot is saved in "+taxonomy+".seqLenHis.png")
    print("close the plot window to continue.\n")
    plt.show()
    # test if the difference of sequences length is too big
    percentlen=(df.max()[2]-df.min()[2])/df.max()[2]
    if percentlen > 0.1:
        test=input("the difference of length of the sequences is too big, would you choose a subset of these proteins? yes/no:\n")
        if test.lower().replace(' ','')=="yes":
            # choose a subset
            test=input("would you like to choose sequences with length between median ± 100 ? yes/no:\n")
            if test.lower().replace(' ','')=="yes":
                print("generate subset of sequences with length in median ± 100.\n")
                mini=int(df.median()[0]-100)
                maxi=int(df.median()[0]+100)
            else:
                # user input a range of length
                while True:
                    print("Please choose a range of length.\n")
                    mini=input('minimum:\n').replace(' ', '')
                    maxi=input("maximum:\n").replace(' ', '')
                    try:
                        int(mini)
                        int(maxi)
                    except:
                        print("invalid number. Please type in integer number. Try again.\n")
                        continue
                    if(maxi>mini):
                        break
                    else:
                        print("minimum number should be less than maximum number. Try again.\n")
                print("generate subset of sequences with length between "+str(mini)+" and "+str(maxi)+".\n")
                sub_df=df[(df['length']>mini) & (df['length']<maxi)]
                label=taxonomy+'_'+mini+'_'+maxi
                with open(label+'_length.fa','w') as out:
                    for i in list(range(0,len(sub_df),1)):
                        out.write('>'+df['protein'][i]+"\n")
                        out.write(df['raw_seq'][i])
        else:
            label=taxonomy+'_whole'
            print("align for all the proteins.\n")
            os.system('cp '+taxonomy+'.fa '+label+'_length.fa')
    else:
        label=taxonomy+'_whole'
        os.system('cp '+taxonomy+'.fa '+label+'_length.fa')
    return label



# determine, and plot, the level of protein sequence conservation
def conservation(dataset,label):
    protein=dataset[0]
    taxonomy=dataset[1]
    # protein sequence alignment
    print("aligning sequences...\n")
    os.system('clustalo -i '+label+'_length.fa -o '+label+'.align.msf --outfmt=msf --threads=100')
    # basic information about sequences in an input multiple sequence alignment
    print("summarizing alignment results...\n")
    os.system('infoalign '+label+'.align.msf -out '+label+'.infoalign')
    # sort by similarity from low to high, indecating conservation from high to low
    os.system('head -n 1 '+label+'.infoalign > '+label+'.infoalign_sort')
    os.system('tail -n +2 '+label+'.infoalign |  sort -k8,8nr >> '+label+'.infoalign_sort')
    # determine the level of conservation and display on the screen
    print("result of alignment sorted by similarity:\n")
    df = pd.read_csv(label+'.infoalign_sort', sep="\t")
    print(df)
    # plot the the level of conservation
    plot="plotcon -sformat msf "+label+".align.msf -winsize 20 -gsubtitle '"+protein+" sequences in "+taxonomy+"' -graph png -goutfile "+taxonomy
    os.system(plot)
    print("plot image file: "+label+".png\n")
    # display the plot
    try:
        os.system('display '+'taxonomy+".png')
    except: # sometimes "display" just doesn't work
        try:
            plot="plotcon -sformat msf "+label+".align.msf -winsize 20 -gsubtitle '"+protein+" sequences in "+taxonomy+"' -graph xterm -goutfile "+taxonomy
            os.system(plot)
        except:
            print("your terminal can't open image. The image is in "+label+".png. You can open later yourself.\n")


