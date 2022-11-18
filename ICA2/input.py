def pro_tx_seq():
    import os
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
            ecount='esearch -db Protein -query "(('+protein+')[Protein Name] AND '+taxonomy+'[Organism]" | grep "Count" > count.txt'
            os.system(ecount)
            count=open("count.txt").read()
            count=int(count.replace('<Count>','').replace('</Count>\n','').strip())
            # too many sequences, ask user 
            if count>1000:
                print("so many sequences! ("+str(count)+"). Continue downloading? yes/no:\n")
                test=input()
                if test=='no':
                    print("Bye!\n")
                    quit()
            #download sequences
            fastafetch='esearch -db Protein -query "(('+protein+')[Protein Name] AND '+taxonomy+'[Organism]" | efetch -format fasta > '+taxonomy+'.fa'
            os.system(fastafetch) 
            print("protein sequences inside "+taxonomy+".fa file.")
            # create a file containing headers of all sequences
            os.system("grep '>'"+taxonomy+".fa >> "+taxonomy+"_title")
            break
        else:
            print("invalid taxonomy ID, try again.\n") # invalid, continue retyping loop
            continue
    return [protein,taxonomy]
