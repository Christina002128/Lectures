#!/usr/bin/python3
def protein_content(protein,amino):
  result=protein.upper().count(amino.upper())/len(protein)
  return result*100

assert round(protein_content("MSRSLLLRFLLFLLLLPPLP", "M")) == round(5)
assert round(protein_content("MSRSLLLRFLLFLLLLPPLP", "r")) == round(10)
assert round(protein_content("MSRSLLLRFLLFllllPPLP", "L")) == round(50)
assert round(protein_content("MSRSLLLRFLLFLLLLPPLP", "Y")) == round(0)


def pro_content_list(protein, amino = ['A', 'I', 'L', 'M', 'F', 'W', 'Y', 'V']):
  if set(amino)!=amino:
    print("not unique amino acid! try again")
    return
  count = 0
  for i in amino:
    count += protein.upper().count(i.upper())/len(protein)
  return count*100

pro_content_list("MSRSLLLRFLLFLLLLPPLP", ['F', 'S', 'L'])
pro_content_list("MSRSLLLRFLLFLLLLPPLP")

assert round(pro_content_list("MSRSLLLRFLLFLLLLPPLP", ["M"])) == 5
assert round(pro_content_list("MSRSLLLRFLLFLLLLPPLP", ['F', 'S', 'L'])) == 70
assert round(pro_content_list("MSRSLLLRFLLFLLLLPPLP")) == 65

def dna_undef(dna,thre=0):
  ATGC=dna.upper().count("A")+dna.upper().count("T")+dna.upper().count("G")+dna.upper().count("C")
  undef=(len(dna)-ATGC)/len(dna)
  if undef >= thre:
    return True
  else :
    return False
## better is just return logic:  return undef >= threshold
assert dna_undef("atgctgcaccgtffr", 0.1)
assert dna_undef("atgctgcaccgtffr")
assert dna_undef("atgctgcaccgtact")


def kmer_count(dna,k=2,n=2):
  if k > len(dna) :
   return "Sorry, your kmer length is longer than your DNA (" + str(len(dna)) +" bases)." 
  if k < 2 or k > 50 :
     return "Sorry, inappropriate kmer length, only 2 to 50 accepted here."
  kmer=[]
  returnkmer=[]
  for i in list(range(0,len(dna)-int(k)+1,1)):
    kmer.append(dna[i:i+k])
  for m in set(kmer):
    if sdna.count(m) > n:
      print(m,sdna.count(m))
      returnkmer.append(m)
  return returnkmer

kmer_count("ATGCATCATGATGCTGATCG",2,2)


# user call
def kmer_count():
  dna=input("the sequence of interest:\n")
  k=input("the kmer length for analysis:\n")
  n=input("the threshold frequency of kmers found (i.e. the \"more than this number\" value):\n")
  k=int(k)
  n=int(n)
  kmer=[]
  returnkmer=[]
  for i in list(range(0,len(dna)-int(k)+1,1)):
    kmer.append(dna[i:i+k])
  for m in set(kmer):
    if sdna.count(m) > n:
      print(m,sdna.count(m))
      returnkmer.append(m)
  return returnkmer
