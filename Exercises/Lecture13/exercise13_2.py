#2
sdna="ATGCATCATG"
k=2 # kmer size
n=2 # more than this number found
kmer=[]
for i in list(range(0,len(sdna)-k+1,1)):
  kmer.append(sdna[i:i+k])
for k in set(kmer):
  if sdna.count(k) > n:
    print(k)

print(kmer)


#3
dnas=['ATTGTACGG', 'AATGAACCG', 'AATGAACCC', 'AATGGGAAT']

for i  in list(range(0,3,1)):
  count=0
  for k in list(range(0,8,1)):
    if dnas[i][k] == dnas[i+1][k]:
      count+=1
  print("seq",i+1,'and seq',i+2,', the percentage of identical positions is',100*count//9,'%') 



