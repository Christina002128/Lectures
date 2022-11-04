  
data=open("data.csv").read().rstrip().split("\n")

for x in data:
  if x.split(',')[0] == "Drosophila melanogaster":
    print(x.split(',')[0],x.split(',')[2])

for x in data:
  if len(x.split(',')[1]) > 90 and len(x.split(',')[1]) < 110:
    print(x.split(',')[2],'length of seq:',len(x.split(',')[1]))

for x in data:
  at= (x.split(',')[1].count('a')+ x.split(',')[1].count('t'))/len(x.split(',')[1])
  if  at < 0.5 and int(x.split(",")[3])>200:
    print(x.split(',')[2]+":AT content< 0.5:",at,", expression level> 200:",x.split(',')[3])
#Print out the gene names begins with "k" or "h" except those belonging to Drosophila melanogaster
for x in data:    
  if x.split(',')[0] =="Drosophila simulans":
    if x.split(',')[2].startswith("k") or x.split(',')[2].startswith("h"):
      print(x.split(',')[2])

for x in data:
  at= (x.split(',')[1].count('a')+ x.split(',')[1].count('t'))/len(x.split(',')[1])
  if at>0.65:
    print(x.split(',')[2],'AT content is:high')
  if at>0.45 and at<0.65:
    print(x.split(',')[2],'AT content is:medium')
  if at<0.45:
    print(x.split(',')[2],'AT content is:low')

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
print("kmer size",kmer,"number",n,"count",k)

#3
dnas=['ATTGTACGG', 'AATGAACCG', 'AATGAACCC', 'AATGGGAAT']

for i  in list(range(0,3,1)):
  count=0
  for k in list(range(0,8,1)):
    if dnas[i][k] == dnas[i+1][k]:
      count+=1
  print("seq",i+1,'and seq',i+2,', the percentage of identical positions is',100*count//9,'%') 



