#!/usr/local/bin/python3
def profile():
  profile={}
  profile['name']=input("What's your name?\n")
  profile['age']=input("How old are you?\n")
  profile['color']=input("What is your favourite colour?\n")
  profile['likepy']=input("Do you like Python?\n")
  profile['flat']=input("The world is flat: True or False?\n")
  if profile['flat']=="True":
    print("worng! world is round!\n")
  print(profile)


def transla(dna):
  gencode = {
  'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
  'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
  'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
  'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
  'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
  'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
  'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
  'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
  'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
  'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
  'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
  'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
  'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
  'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
  'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
  'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
  pro1=[]
  for i in list(range(0,len(dna),3)):
    if i+3 > len(dna):break  
    codon=dna[i:i+3]
    if 'N' in codon:
      pro1.append('N')  # if N in codon, then amino acid = N
    else:
      pro1.append(gencode.get(codon))
  print(dna)
  print("start at base 1 :",''.join(pro1))
  pro2=[]
  for i in list(range(1,len(dna)-1,3)):
    if i+3 > len(dna):break
    codon=dna[i:i+3]
    if 'N' in codon:
      pro2.append('N')
    else:
      pro2.append(gencode.get(codon))
  print("start at base 2 :",''.join(pro2))
  pro3=[]
  for i in list(range(2,len(dna)-2,3)):
    if i+3 > len(dna):break
    codon=dna[i:i+3]
    if 'N' in codon:
      pro3.append('N')
    else:
      pro3.append(gencode.get(codon))
  print("start at base 3 :",''.join(pro3))
  pro=[''.join(pro1),''.join(pro2),''.join(pro3)]
  return pro


profile()

try:
  transla("ATCGCTAGCTAGCCG")
except: # if fail, print sth, and keep running rest of script
  print("error in transla(\"ATCGCTAGCTAGCCG\")")
try:
  transla("ATCGCTAGCNAGCCG")
except: # if fail, print sth, and keep running rest of script
  print("error in transla(\"ATCGCTAGCNAGCCG\")")
try:
  transla("ATNNNNGCTAGCNAGCCG")
except: # if fail, print sth, and keep running rest of script
  print("error in transla(\"ATNNNNGCTAGCNAGCCG\")")

