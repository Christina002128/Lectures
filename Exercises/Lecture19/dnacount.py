def dna_count(file="ecoli.txt"):
  ecoli = open(file).read().replace('\n', '').upper()[0:100000]
  window = input("window size? \n")
  content = input("AT content or GC content?\n")
  base_ran=input("what portion (base range) of the genome you want analysed?\n")
  eco=ecoli[0:int(base_ran)]
  at = []
  counter =0
  for start in range(len(eco) - window):
    counter+=1
    print(counter)
    win=ecoli[start:start+window]
    at.append((win.count('A')+win.count('T')) / len(win))
  if content=='AT':
    at=at
  if content=='GC':
    for i in range(0,len(at),1):
      at[i]=1-at[i]
  plt.figure(figsize=(20,10))
  plt.plot(at, label=content+" rep")
  plt.ylabel('Overrepresentation')
  plt.xlabel('Position on genome')
  plt.suptitle("Base composition in the E coli genome",fontsize=20) # super title
  plt.title("Window size of "+str(window),fontsize=14)
  plt.legend()
  plt.savefig(file+base_ran+content+str(window)+".png",transparent=True)
  plt.show()

dna_count(file="ecoli.txt")
20
AT
5000
