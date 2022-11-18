#!/usr/local/bin/python3
import os
import ICA2_function

# run pro_tx_seq() function to let user input protein name and taxinomic group ID
# check if the input invalid or retreved 0 result, let user type again (in loop)
dataset=ICA2_function.pro_tx_seq()

# run species() function to check how many species contained in the dataset
dic_species=ICA2_function.species(dataset[1]) # return dictionary { accession : species }

# ask user to continue or not
condition=input("Do you want to continue analysing? yes/no:\n")
if condition.lower().replace(' ','')=="no":
    print("Bye!\n")
    quit()


# count length of each sequence
df_seq=ICA2_function.seq_len(dataset[1])
# generate subset based on length of protein sequences
label=seq_len_subset(df_seq,dataset[1])

# determine, and plot, the level of protein sequence conservation
ICA2_function.conservation(dataset,label)
condition=input("Do you want to continue finding motifs?\n")
if condition.lower().replace(' ','')=="no":
    print("Bye!\n")
    quit()

# output fasta file for each sequence for motifs scanning
seq_out(dataset[1])

# choose
find_motifs(taxonomy,label,df)


#3. scan with motifs from PROSITE database , names of motifs
ICA2_function.seq_out(dataset[1])
ICA2_function.find_motifs(dataset[1])






cut -d ' ' -f 3 
os.system("cut -d ' ' -f 3 *count | sort > "+taxonomy+".motifs.txt")

elink -target protein | efilter -organism mammal -source refseq | efetch -format fasta

esearch -db taxonomy -query "mammals" | efetch -format xml  > mammals.xml
esearch -db protein -query "(pyruvate dehydrogenase)[Protein Name] AND txid4890[Organism]"| grep "Count"


esearch -db protein -query "(pyruvate dehydrogenase)[Protein Name] AND txid4890[Organism]" | efetch -format acc | wc -l
esearch -db protein -query "(ABC transporters)[Protein Name] AND mammals[Organism] " | efetch -format acc | wc -l


esearch -db protein -query "((pyruvate dehydrogenase) AND txid4890[Organism])" | efetch -format fasta  > proteins.fa
clustalo -i txid4890.fa -o txid4890.align.msf --outfmt=msf --threads=100
infoalign txid4890.align.msf -out txid4890.infoalign
sort -k10 txid4890.infoalign > txid4890.infoalign.sort
plotcon -sformat msf txid4890.align.msf -winsize 25 -gsubtitle 'pyruvate dehydrogenase sequences in txid4890' -graph svg -goutfile txid4890


plotcon -sformat msf txid40674.align.msf -winsize 25 -gsubtitle 'pyruvate dehydrogenase sequences in txid40674' -graph png -goutfile txid40674


patmatmotifs -full -sequence test.fa1 -outfile txid4890.patmatmotifs

# Create an ambiguous consensus sequence
consambig -sequence txid4890.align.msf -outseq txid4890.aligned.cons
dataset=["ABC transporters",'txid40674']

'esearch -db protein -query "('+protein+')[Protein Name] AND '+taxonomy+'[Organism] NOT PARTIAL'

def find_motifs(taxonomy,label,df):
    motif='for fn in *.seq.fa \n do patmatmotifs -full -sequence ${fn} -outfile ${fn//.seq.fa/}.patmatmotifs & \n done \n wait'
    os.system(motif)
    # count the number of each motif


    with multiprocessing.Pool() as pool:
        pool.map(moti, proteins)


    # count the number of each motif
    proteins=list(df['pro_number'])
    for protein in proteins:
        protein=protein.replace('|','\|')
        os.system('ls '+protein+'*')
        query_motif='patmatmotifs -full -sequence '+protein+'.seq.fa -outfile '+protein+'.patmatmotifs'
        os.system(query_motif)
