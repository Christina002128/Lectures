#!/usr/bin/python3
#1
dna=open('input.txt').read()
clean=dna.replace(dna[0:14],'')
outfile=open('clean_dna.txt','w')
outfile.write(clean)
outfile.close()
print(open('clean_dna.txt').read().rstrip()) # remove blank line at the end
print('the number of lines is',len(open('clean_dna.txt').read().rstrip()))

import subprocess
subprocess.call('cat clean_dna.txt', shell=True)

#2
dna=open('genomic_dna2.txt').read()
expos=open('exons.txt').read().rstrip().replace(',','\n').split(',')
exon=[]
for i in [1,3,5,7]:
    exon1=dna[int(expos[i-1])-1:int(expos[i])]
    exon.append(exon1)

with open("genomic_exon.txt",'w') as ex:
    ex.write(''.join(exon))

#3
import os, shutil
shutil.copy("../Lecture11/AJ223353_noheader.fasta2",".")
sorted(os.listdir())
dna=open("AJ223353_noheader.fasta2").read().rstrip()
k=0
with open("segment.fasta","w") as file:
  print("number\tseg_sequence\tGC(%)\n")
  file.write("number\tseg_sequence\tGC(%)\n")
  for i in list(range(0,len(dna)-3,3)):
    k=k+1
    seg=dna[i:int(i)+30]
    GC=100*(seg.count('G')+seg.count('C'))/len(seg)
    file.write(str(k)+'\t'+seg+'\t'+str(GC)+'\n')

open("segment.fasta").read()

#4
shutil.copy("/localdisk/data/BPSM/Lecture12/*.dna",".")
os.system("wc -l *dna")
os.system("mkdir dna_files")
os.system("mv *.dna ./dna_files")
for file_name in os.listdir('dna_files/'):
  if file_name.endswith('.dna'):
    f'Reading sequences from {file_name}'
    dna_file = open('dna_files/' + file_name)
# This loop then looks at each line and gets the length
# Note the indentation...
    for line in dna_file : 
      dna = line.rstrip('\n') 
      length = len(dna) 
      print(f'\tFound a DNA sequence with length {str(length)}' )

for file_name in  sorted(os.listdir('dna_files')) : 
# Check if the file name ends with .dna
  if file_name.endswith('.dna') : 
    print('Reading sequences from ' + file_name) 
# Open the file and process each line
    dna_file = open('dna_files/' + file_name) 
# Calculate the sequence length 
    for line in dna_file :
      dna = line.rstrip('\n') 
      length = len(dna) 
      print('\tsequence length is ' + str(length)) 
# Go through each size bin and check if the sequence belongs in it
      for bin_lower in list(range(100,1000,100)) : 
        bin_upper = bin_lower + 99 
        if length >= bin_lower and length <= bin_upper : 
          print('\t\tbin is ' + str(bin_lower) + ' to ' + str(bin_upper)) 

for bin_lower in list(range(100,1000,100)) : 
    bin_upper = bin_lower + 99 
    bin_directory_name = str(bin_lower) + '_' + str(bin_upper) 
    os.mkdir(bin_directory_name) 

