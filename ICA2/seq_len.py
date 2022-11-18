def seq_len(taxonomy):
    fas=open(taxonomy+'.fa').read().strip()
    full_seq=fas.split('>')[1:] # delete first empty string
    seq=[]
    pro=[]
    seqlen=[]
    for i in list(range(0,len(full_seq),1)): 
        # count sequence length for each sequence
        devide=full_seq[i].find("\n")
        pro[i]=full_seq[i][0:devide]
        seq[i]=full_seq[i][(devide+1):].replace('\n','')
        seqlen[i]=len(seq[i])  
    df={} # data frame contain the protein name, sequence, sequence length
    df=pd.DataFrame(pro,columns=["protein"])
    df=pd.DataFrame(seq,columns=["sequence"])
    df=pd.DataFrame(seqlen,columns=["length"])
    return df
