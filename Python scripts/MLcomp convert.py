import numpy as np
import pandas as pd
#import sklearn
#from sklearn.preprocessing import scale
#from sklearn.lda import LDA
#from sklearn.decomposition import PCA
#from sklearn import RandomForestClassifier
#from cudatree import RandomForestClassifier
#from hybridforest import RandomFncgorestClassifier
#from PyWiseRF import WiseRF    #premium, part of wise.io and Anaconda Pro

def df2mlcomp(train, target, fname='mlcomp.txt'):
    s = ''
    rw = 0
    for header, row in train.T.iteritems():
        s += str(target[rw])  #headerd
        cl = 0
        for k,v in row.iteritems():
            cl += 1
            if np.isnan(v):
                raise Exception("Value "+str(v)+" at col "+k+" ("+str(cl)+") and row "+header+" ("+str(rw)+") is NaN!")
            s += " "+str(cl)+":"+str(v)  #k
        #print s
        s += '\n'  #\r
        rw += 1
    text_file = open(fname, 'w')
    text_file.write(s)
    text_file.close()

if __name__ == '__main__':

    #%cd 'D:\Workspace\Breath Analysis\'
    a = []
    
    file = 'concentrationArray.txt'
    df = pd.read_csv(file, sep=',', header=None).astype(float)
    sample,molecule=df.shape
    df.ix[:,molecule]='-1'
    df.ix[sample/2:,molecule]='1'
    #df = df.ix[:, 0:molecule-1]   #normalize column count for renaming
    df = df.dropna().T.dropna().T
    a.append(df)
    
    df = pd.concat(a)

    train = df.ix[:,:molecule-1]#[:,:-1]
    target = df.ix[:,molecule].values#df.ix[:,-1:].values.T[0]
    mlcomp = df2mlcomp(train, target, 'binary.txt')

