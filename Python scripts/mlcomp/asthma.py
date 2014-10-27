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

    #freqs = set(np.loadtxt('Wavenumber.txt').astype(str))
    a = []
    elements=37
    
    for grp in range(0, 1+1):  #0: Healthy, 1: Asthma, 2: CF
        for smpl in range(1, 2+1):
            file = 'data of 26-8-2014\data low p only\limG%dS%d.txt' % (grp, smpl)
            df = pd.read_csv(file, delim_whitespace=True, header=0).astype(float).T
            df = df.ix[:, 0:elements-1]   #normalize column count for renaming
            #df.columns = freqs
            df = df.dropna().T.dropna().T
            df['result'] = grp if (grp > 0) else -1 
            a.append(df)
    
    df = pd.concat(a)
    #df.to_csv('combined.csv')

    #processed = sklearn.preprocessing.scale(df)
    #df = pd.DataFrame(processed, index=df.index, columns=df.columns)
    #df.to_csv('scaled.csv')

    train = df.ix[:,:-1]
    target = df.ix[:,-1:].values.T[0]
    mlcomp = df2mlcomp(train, target, 'binary.txt')

'''    
    #clf = sklearn.ensemble.RandomForestClassifier(n_estimators=5, bootstrap=False, verbose=True, n_jobs = 4)   #50
    #clf = RandomForestClassifier(n_estimators=5, bootstrap=False, n_jobs = 1, n_gpus = 0)  #, verbose=True, cpu_classifier = WiseRF
    #clf = LDA() #sklearn.lda.
    clf = PCA()

    clf.fit(train.values, target)	#, bfs_threshold = 4196
    #predicted = clf.predict(target)
    predicted = clf.components_
    np.savetxt('predicted.csv', predicted)
'''
