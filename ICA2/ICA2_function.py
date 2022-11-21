#!/usr/local/bin/python3

#############################################################################
#Script Name	: ICA2_function.py                                                                                          
#Description	: Contain all functions. Run with "ICA2_script.py" file.                                                                            
#Date           : 2022_11_19                                                                                       
#Author       	: B222908                                          
#help document  : B222908-2022.ICA2.manuals.pdf  
#version        : 2.1                                     
############################################################################

# import model
import os
import pandas as pd
import re
import matplotlib.pyplot as plt
import shutil

# user input protein name and taxinomic group
def pro_tx_seq():
    while True: # a loop till correct input
        protein=input("Please type in protein name:\n").strip()
        taxonomy=input("Please type in a taxonomic name or a taxonomy ID. (taxonomy ID for example: txid4890)\n").lower().strip()
        # test if taxonomy input is valid
        if re.search('txid', taxonomy): # if input is taxonomy ID, then
            # test if id is integer
            taxonomy=taxonomy.replace(' ','') # delete spaces
            id=taxonomy.replace("txid",'') # save id number
            try:
                int(id) 
            except:
                print("invalid taxonomy ID, try again.\n")
                continue  # skip to the next loop, input again
            print('query for "'+protein+'" protein and "'+taxonomy+'" group.\n')
        else: # if input is only integer
            try:
                int(taxonomy.replace(' ','')) 
                taxonomy='txid'+taxonomy.replace(' ','') # if integer, add 'txid' in the front
                print('query for "'+protein+'" protein and "'+taxonomy+'" group.\n')
            except: # if not, then it's taxonomic name
                print('generate query for "'+protein+'" protein and "'+taxonomy+'" taxonomic group.')
        #partial or not partial search
        par=input("would you like to contain partial proteins or not? yes/no:\n")
        par=par.replace(' ','').lower()
        # count the number of sequences searched from NCBI
        # check if count.txt file exsist, if not, create file
        try: 
            open("count.txt").read()
        except:
            with open("count.txt",'w') as c:
                c.write('')
        if re.search('n',par): # search not partial proteins
            print("searching no partial protein sequences from ncbi...")
            ecount='esearch -db protein -query "'+protein+'[Protein] AND '+taxonomy+'[organism] NOT PARTIAL" | grep "Count" > count.txt'
            os.system(ecount)
            par_query='NOT PARTIAL'  # for esearch query
            par="no_par"  # for filename label
        else: # search all proteins
            print("searching all protein sequences from ncbi...")
            ecount='esearch -db protein -query "'+protein+'[protein] AND '+taxonomy+'[organism]" | grep "Count" > count.txt'
            os.system(ecount)
            par_query='' # for esearch query
            par='all' # for filename label
        count=open("count.txt").read()
        count=int(count.replace('<Count>','').replace('</Count>\n','').strip()) 
        # check if the number of sequences too high and ask user if they want to download
        if count>1000: 
            print("Too many sequences("+str(count)+")!  Continue downloading? yes/no:\n")
            test=input()
            if re.search('n',test.lower()):
                test2=input("Do you want to retype a group of proteins? yes/no:\n")
                if re.search('y',test.lower()):
                    continue
                else:
                    print("Bye!\n")
                    quit() # quit the script
        # check if the retrieve result is 0
        if count==0:
            print("0 result searched. Try again.\n")
            continue  # skip to the next loop, input again
        else:
            print(str(count)+" sequences downloading...")
            # create label for naming the file
            tax=taxonomy.replace(' ', '_')
            pro=protein.replace(' ', '_')
            protaxpar_label=pro+'_'+tax+'_'+par
            #download sequences
            fastafetch='esearch -db protein -query "'+protein+'[protein] AND '+taxonomy+'[organism] '+par_query+'" | efetch -format fasta > '+protaxpar_label+'.fa'
            os.system(fastafetch) 
            print("protein sequences is inside file: "+protaxpar_label+".fa")
            break
    return [protein,taxonomy,protaxpar_label] # return a list




# Make a dataframe containing contents from fasta file and calculate the sequence length
def create_df(protaxpar_label):
    fas=open(protaxpar_label+".fa").read().strip()
    # get all headers 
    fas_list=fas.split('\n')
    headers=[]
    for line in fas_list:
        if '>' in line:
            headers.append(line)
    # get headers content: accession, protein_name, genus, species
    acc=[]
    pro_name=[]
    genus=[]
    species=[]
    full_species=[]
    for pro in headers:
        # for sequences with regular header:  >accession protein_name [species]
        if pro.find("["): 
            acc.append(re.search(r'>(\S*)\s(.*)\[(.*)\]',pro).group(1).strip()) # accession
            pro_name.append(re.search(r'>(\S*)\s(.*)\[(.*)\]',pro).group(2).strip()) # protein_name
            spec=re.search(r'>(\S*)\s(.*)\[(.*)\]',pro).group(3).strip() # [species]
            genus.append(spec.split(' ')[0]) # genus name
            full_species.append(spec) # species name
            sp_list=spec.split(' ')[1:]
            species.append(' '.join(sp_list))
        # for sequences with strange header, all contents equals to header
        else: 
            acc.append(pro.split(' ')[0].replace(">",''))
            full=re.search(r'(>\S*)\s(.*$)',pro).group(2).strip()
            proname.append(full)
            species.append(full)
            full_species.append(full)
            genus.append(full)
            species.append(full)
    # get sequence and sequence length
    full_seq=fas.split('>')[1:] # delete first empty string
    seq=[]
    raw_seq=[]
    seqlen=[]
    for i in list(range(0,len(full_seq),1)): 
        # count sequence length for each sequence
        devide=full_seq[i].find("\n") # first index of '\n' to divide header and sequence
        seq.append(full_seq[i][(devide+1):].replace('\n',''))
        raw_seq.append(full_seq[i][(devide+1):])
        seqlen.append(len(seq[i]))
    df={} # data frame
    df=pd.DataFrame({'protein_accession':acc,"protein_name":pro_name,
    "genus":genus,"species":species,"full_species_name":full_species,
    "length":seqlen,"header":headers,"sequence":seq,'raw_seq':raw_seq})
    df=df.sort_values('genus',ascending=True) # sort by genus
    df.to_csv(protaxpar_label+"_summary.csv",sep='\t',index=False)
    return df



# Check how many species contained in the dataset
def count_species(protaxpar_label):
    # load dataframe
    df=pd.read_csv(protaxpar_label+"_summary.csv", sep="\t")
    # count the number of proteins, genus and species
    protein_numbers=len(list(df['protein_accession']))
    genus_numbers=len(set(list(df['genus'])))
    species_numbers=len(set(list(df['full_species_name'])))
    sp_percent=round(protein_numbers/species_numbers,2)
    ge_percent=round(protein_numbers/genus_numbers,2)
    print("\nResult of retrieved data:")
    print("1. The result contains "+str(protein_numbers)+" sequences, "+str(genus_numbers)+" kinds of genus and "+str(species_numbers)+" kinds of species.")
    print("2. About "+str(sp_percent)+" proteins per species, "+str(ge_percent)+" proteins per genus.")
    # write species summary to csv file
    sub_df=df[['protein_accession',"protein_name",'genus','species']]
    sub_df.to_csv(protaxpar_label+"_summary_species.csv",sep='\t',index=False)
    print('3. Result inside "'+protaxpar_label+'_summary_species.csv" file.\n')
    print("Here are the first 10 rows of the result sorted by genus name.")
    print(sub_df.head(10)) # print first 10 rows



# choose a appropriate protein subset for alignment by considering the sequence length
def seq_len_subset(dataset):
    # load dataframe
    df=pd.read_csv(protaxpar_label+"_summary.csv", sep="\t")
    print("\nSummary of sequences length:")
    print(df.describe()) # print out calculations of length
    # display histogram of sequence length
    print("\nHistogram of sequences length.")
    print("(If it fails to generate plot, please restart a terminal window and run the program again.)\n")
    y=list(df['length'])
    n, bins, patches = plt.hist(y, 50, color='c')
    plt.grid(True)
    plt.xlabel('Sequence length')
    plt.ylabel('Number of proteins')
    plt.title('Histogram of sequences length')
    plt.savefig(dataset[2]+".seqLenHis.png")
    print("The plot is saved in "+dataset[2]+".seqLenHis.png\n")
    try:
        print("Showing histogram. Close the plot window to continue.(If it failed, try open it yourself: "+dataset[2]+".seqLenHis.png)\n")
        plt.show() # sometimes it failed, restart terminal will fix it
    except:
        print("Failed to open the image, try open the file yourself: "+dataset[2]+".seqLenHis.png\n")
    # test if the difference of sequences length is too big: (max-min)/max > 0.2
    percentlen=(df.max()['length']-df.min()['length'])/df.max()['length']
    if percentlen > 0.2:
        test1=input("the difference of length of the sequences is quite big, would you choose a subset based on length range? yes/no:\n")
        if re.search('y','test1.lower()'):
            # choose a subset
            test2=input("would you like to choose sequences with length within median ± sd ? yes/no:\n")
            if re.search('y','test2.lower()'):
                mini=int(df.median()[0]-df.std()[0])
                maxi=int(df.median()[0]+df.std()[0])
            else:
                # user input a range of length, test the validity
                while True:
                    print("Please choose a range of length.\n")
                    mini=input('minimum:\n').replace(' ', '')
                    maxi=input("maximum:\n").replace(' ', '')
                    try:
                        mini=int(mini)
                        maxi=int(maxi)
                    except:
                        print("invalid number. Please type in an integer number. Try again.\n")
                        continue
                    if(maxi>mini):
                        break
                    else:
                        print("minimum number should be less than maximum number. Try again.\n")
                        continue
            print("generate subset of sequences with length between "+str(mini)+" and "+str(maxi)+".\n")
            sub_df=df[(df['length']>mini) & (df['length']<maxi)]
            sub_label=dataset[2]+'_'+str(mini)+'_'+str(maxi)
            sub_df.to_csv(sub_label+"_subgroup.csv",sep='\t',index=False)
            with open(sub_label+'_len.fa','w') as out:
                for i in list(range(0,len(sub_df),1)):
                    out.write(df['header'][i]+"\n")
                    out.write(df['raw_seq'][i])
        else:
            sub_label=dataset[2]+'_whole'
            print("align for all the proteins.\n")
            shutil.copyfile(dataset[2]+'.fa',sub_label+'_len.fa')
    else:
        sub_label=dataset[2]+'_whole'
        print("align for all the proteins.\n")
        shutil.copyfile(dataset[2]+'.fa',sub_label+'_len.fa')
    return sub_label


dataset=['pyruvate dehydrogenase','ascomycete fungi',"pyruvate_dehydrogenase_ascomycete_fungi_all"] 
sub_label='pyruvate_dehydrogenase_ascomycete_fungi_no_par_whole'
protaxpar_label
# determine, and plot, the level of protein sequence conservation
def conservation(dataset,sub_label):
    protein=dataset[0]
    taxonomy=dataset[1]
    # sequence alignment
    print("Aligning sequences...\n")
    os.system('clustalo -i '+sub_label+'_len.fa -o '+sub_label+'.align.msf --outfmt=msf --threads=100')
    # basic information for the alignment
    print("Summarizing alignment results...\n")
    os.system('infoalign '+sub_label+'.align.msf -out '+sub_label+'.infoalign')
    # dataframe sorted by similarity from high to low
    #aligninfo=pd.read_csv(sub_label+'.infoalign', sep="\t",na_values=[''])
    # fix unmatched header
    #aligninfo=aligninfo.sort_values('Similar',ascending=False)
    #aligninfo.to_csv(sub_label+".infoalign_sort",sep='\t')
    #print('Here is the first 10 rows of aligment summary sorted by similarity.')
    #aligninfo.head(20)
    non=os.system("head -1 "+sub_label+'.infoalign')
    non=os.system("sort -t$'\t' -k7,7nr "+sub_label+'.infoalign | head -10')
    print("Full result inside \""+sub_label+".infoalign\" file.")
    # generate conservation sequence
    os.system('cons '+sub_label+'.align.msf '+sub_label+'_align.cons')
    aligncons=open(sub_label+'_align.cons').read()
    print("Consensus sequence inside \""+sub_label+"_align.cons\" file.\nHere is consensus sequence:")
    print(aligncons)
    # plot the the level of conservation
    plot="plotcon -sformat msf "+sub_label+".align.msf -winsize 20 -gsubtitle '"+protein+" sequences in "+taxonomy+"' -graph png -goutfile "+sub_label
    os.system(plot)
    print("conservation image file: "+sub_label+".png\n")
    # display the plot
    non=os.system('display '+sub_label+'.png') # sometimes "display" just doesn't work
    print("If your terminal can't open image, try open the image yourself. Sometimes restart terminal help\n")


# output fasta file for each sequence
# only
def seq_out(sub_label):
    fas=open(sub_label+'_len.fa').read().strip()
    full_seq=fas.split('>')[1:] # delete first empty string
    filename=[]
    full_fa=[]
    for i in list(range(0,len(full_seq),1)): 
        # output fasta file for each sequence
        filename.append(full_seq[i].split(' ')[0])  # set sequence name as file name
        full_fa.append('>'+full_seq[i])   # add separator('>') back
        with open(filename[i]+'_'+sub_label+'_len_single_seq.fa',"w") as outfile:
            outfile.write(full_fa[i]) # write in single sequence

# scan each sequence with motifs from PROSITE database
def find_motifs(sub_label):
    # load dataframe
    df=pd.read_csv(sub_label+"_subgroup.csv", sep="\t")
    # scan motif for each sequence, parallel
    motif='for fn in *'+sub_label+'_len_single_seq.fa \n do patmatmotifs -full -sequence ${fn} -outfile ${fn//_seq.fa/}.patmatmotifs & \n done \n wait'
    os.system(motif)
    # get all the motif names from all files
    count_motifs="grep 'Motif =' *"+sub_label+"_len_single.patmatmotifs | cut -d ' ' -f 3 | sort > "+sub_label+"_len.motifs.count"
    os.system(count_motifs)
    # count the number of each motif in total
    motifs_n=open(sub_label+"_len.motifs.count").read()
    motifs_n=motifs_count.strip().split('\n')
    motifs_unique=set(motifs_n) 
    moti_count=[]
    moti_name=[]
    for uni in motifs_unique:
        moti_name.append(uni)
        moti_count.append(motifs_n.count(uni))
    # dataframe for motif counts. sort by the amount 
    motif_df=pd.DataFrame({"motif_name":moti_name,"amount":moti_count})
    motif_df=motif_df.sort_values('amount',ascending=False)
    print("Result of patmatmotifs is in "+sub_label+"patmatmotifs_summary.csv file.\n")
    print(motif_df)
    # summary of motifs result and output
    motif_df.to_csv(sub_label+'patmatmotifs_summary.csv',sep='\t',index=False)


