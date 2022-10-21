#!/usr/bin/python3
# use the python3 in bash,can't use alias or (which python)
# or just use: python3 script_name, without the upper thing
DNA="ATCGATCGATTACGTATAGTATTTGCTATCATACATATATATCGATGCGTTCAT"
print("number of A:",DNA.count('A')) # need use print in script
print("number of T:",DNA.count('T'))
print("AT%:",(DNA.count('A')+DNA.count('T'))/len(DNA))
com_DNA=DNA.replace('A','t')
com_DNA=com_DNA.replace('T','a')
com_DNA=com_DNA.replace('G','c')
com_DNA=com_DNA.replace('C','g')
print("complement of the DNA:",com_DNA.upper())
DNA="ACTGATCGATTACGTATAGTAGAATTCTATCATACATATATATCGATGCGTTCAT"
print("number of GAATTC:",DNA.count('GAATTC'))
print("start point of GAATTC:",DNA.find('GAATTC'))
cut_point=DNA.find('GAATTC')
print("total length",len(DNA))
print("the length of first fragment",len(DNA[0:cut_point+1])) #exclusive at end
print("the length of second fragment",len(DNA[cut_point+1:]))
#4
DNA="ATCGATCGATCGATCGACTGACTAGTCATAGCTATGCATGTAGCTACTCGATCGATCGATCGATCGATCGATCGATCGATCGATCATGCTATCATCGATCGATATCGATGCATCGACTACTAT"
exon1=DNA[0:63]
exon2=DNA[90:]
print("total length:",len(DNA))
print("first exon",exon1,len(exon1))
print("second exon",exon2,len(exon2))
print("join two exons",exon1+exon2)
join=exon1+exon2
print("length of coding seq:",len(join))
print("percentage of the DNA sequence is coding:",int((len(join)/len(DNA))*100),"%") # round -- int
intron=DNA[63:90]
intron=intron.lower()
print("intron:",intron,len(intron))
print("concatenate:",exon1+intron+exon2)
