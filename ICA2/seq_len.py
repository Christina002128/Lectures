import matplotlib.pyplot as plt
# choose a appropriate protein subset for alignment by considering the sequence length
def seq_len_subset(df,dataset):
    print("\nSummary of sequences length:")
    print(df.describe()) # print out calculations of length
    # display histogram of sequence length
    print("\nHistogram of sequences length.")
    print("(If it fails to generate plot, please restart a terminal window and run the program again.)\n")
    y=list(df['length'])
    n, bins, patches = plt.hist(y, 50, color='b')
    plt.grid(True)
    plt.xlabel('Sequence length')
    plt.ylabel('Number of proteins')
    plt.title('Histogram of sequences length')
    plt.savefig(dataset[2]+".seqLenHis.png")
    print("The plot is saved in "+dataset[2]+".seqLenHis.png\n")
    try:
        print("Showing histogram. Close the plot window to continue.(If it failed, try open it yourself: "+dataset[2]+".seqLenHis.png)\n")
        plt.show() # sometimes it failed, restart terminal will fix it
    except:
        print("Failed to open the image, try open the file yourself: "+dataset[2]+".seqLenHis.png\n")
    # test if the difference of sequences length is too big: (max-min)/max > 0.2
    percentlen=(df.max()['length']-df.min()['length'])/df.max()['length']
    if percentlen > 0.2:
        test1=input("the difference of length of the sequences is quite big, would you choose a subset based on length range? yes/no:\n")
        if test1.lower().replace(' ','')=="yes":
            # choose a subset
            test2=input("would you like to choose sequences with length within median ± sd ? yes/no:\n")
            if test2.lower().replace(' ','')=="yes":
                mini=int(df.median()[0]-df.std()[0])
                maxi=int(df.median()[0]+df.std()[0])
            else:
                # user input a range of length, test the validity
                while True:
                    print("Please choose a range of length.\n")
                    mini=input('minimum:\n').replace(' ', '')
                    maxi=input("maximum:\n").replace(' ', '')
                    try:
                        mini=int(mini)
                        maxi=int(maxi)
                    except:
                        print("invalid number. Please type in an integer number. Try again.\n")
                        continue
                    if(maxi>mini):
                        break
                    else:
                        print("minimum number should be less than maximum number. Try again.\n")
                        continue
            print("generate subset of sequences with length between "+str(mini)+" and "+str(maxi)+".\n")
            sub_df=df[(df['length']>mini) & (df['length']<maxi)]
            sub_label=dataset[2]+'_'+str(mini)+'_'+str(maxi)
            with open(sub_label+'_len.fa','w') as out:
                for i in list(range(0,len(sub_df),1)):
                    out.write(df['header'][i]+"\n")
                    out.write(df['raw_seq'][i])
        else:
            sub_label=dataset[2]+'_whole'
            print("align for all the proteins.\n")
            shutil.copyfile(dataset[2]+'.fa',sub_label+'_len.fa')
    else:
        sub_label=dataset[2]+'_whole'
        print("align for all the proteins.\n")
        shutil.copyfile(dataset[2]+'.fa',sub_label+'_len.fa')
    return sub_label
