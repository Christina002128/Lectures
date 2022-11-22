#!/usr/local/bin/python3
import subprocess
from Bio import Entrez,SeqIO
# set email address
Entrez.email = "s2442072@ed.ac.uk"
Entrez.api_key=subprocess.check_output("echo $NCBI_API_KEY", shell=True).rstrip().decode()

result1=Entrez.read(Entrez.esearch(db='protein',term='COX1[protein] AND mammals[Organism] NOT PARTIAL',retmax=30))
len(result1['IdList']) #13
sum=0
for acc in result1['IdList']:
    gb_f=Entrez.efetch(db='protein',id=acc,rettype='gb')
    record=SeqIO.read(gb_f,'genbank')
    print("COX1",record.id,record.description)
    prolen=len(str(record.seq))
    sum+=prolen
    print("length: "+str(prolen))

ave=round(sum/len(result1['IdList']),2)
print('average length is: '+str(ave))


def entrezSearch(gene,taxon,partial=False,rmax=30):
    if not partial:
        res=Entrez.read(Entrez.esearch(db='nucleotide',term=gene+"[gene] AND "+taxon+"[Organism] NOT PARTIAL",retmax=rmax))
        print('Not partial results searching.')
    else:
        res=Entrez.read(Entrez.esearch(db='nucleotide',term=gene+"[gene] AND "+taxon+"[Organism]",retmax=rmax))
        print('Partial and not partial results searching.')
    print(res['Count'],"results in total, only get first",len(res['IdList']),'results.')
    sum=0
    for acc in res['IdList']:
        gb_f=Entrez.efetch(db='nucleotide',id=acc,rettype='gb')
        record=SeqIO.read(gb_f,'genbank')
        print(gene,record.id,record.description)
        genelen=len(str(record.seq))
        sum+=genelen
        print("length: "+str(genelen))
    ave=round(sum/len(res['IdList']),2)
    print('average length is: '+str(ave))


entrezSearch('cox1','mammals',partial=False,rmax=20)
