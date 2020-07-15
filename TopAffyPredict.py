#!/usr/bin/python3
import random
import numpy as np
import collections
import scipy.stats
global kmer5
global kmer8
global emTable
global psudeoEm
global X
global psTua
global psEm

def walkScore(seq,mapK): 
        old = 'S'
        trans = 0
        total = np.zeros(3)
        for i in range(len(seq)-1):
                if table[diOrder[seq[i:i+2]],i] == 0:
                        tau = psTua
                else:
                        tau = table[diOrder[seq[i:i+2]],i]
                if sum(tableEm[AAL[seq[i]],AAL[seq[i+1]],mapK]) == 0:
                        em = psEm[mapK]
                else:
                        em = tableEm[AAL[seq[i]],AAL[seq[i+1]],mapK]
                total += tau * em
        else:
                if table[diOrder[seq[i+1] + 'END'],i] == 0:
                        tau = psTua
                else:
                        tau = table[diOrder[seq[i+1] + 'END'],i]
                if sum(tableEm[AAL[seq[i]],AAL[seq[i+1]],mapK]) == 0:
                        em = psEm[mapK]
                else:
                        em = tableEm[AAL[seq[i]],AAL[seq[i+1]],mapK]                
                total += tau * em

        return sum(total)

def makeHMM(fileName):
    file = open(fileName, 'r')
    seqs = np.array([list(row.strip().split(',')[2]) for row in file ])

    
    for i in X:
        count = 0 
        for j in range(len(i)-1):
            table[diOrder[(i[j:j+2])],count]+=1
            count+=1
        else:
            table[diOrder[i[j+1] + 'END'],count]+=1
    inv_diOrder = {v: k for k, v in diOrder.items()}
    maxV = np.array(table).max()
    for i in range(maxL):            
           table[:,i] = (table[:,i]/maxV)*(sum(table[:,i])/len(X)) #Table Max
           table[table[:,i] != 0,i] = (.5-(table[table[:,i] != 0,i])*(1-table[table[:,i] != 0,i]))*2#Gini          
    
    psTua = np.mean(table[table != 0])
def read8mers(fileName):
	file = open(fileName, 'r')
	for i in file:
		token = i.strip().split(',')
		kmer8.append([int(j) for j in token[1:]])

if __name__ == "__main__":
        kmer8 = []
        read8mers('Full6merComb.csv')
        SeqFile = 'seqOnlyTBOX.csv'
        name = 'MouseTest'
        hmmSeq = 'mouseSeqs.csv'
        kmerLen = 4**6
        print(name)

        file = open(SeqFile, 'r')
        fin = [row.strip() for row in file ]
        trans = {'-':0 ,'A':1, 'R':2,'N':3, 'D':4,  'C':5,  'E':6,  'Q':7,  'G':8,  'H':9,  'I':10, 'L':11, 'K':12, 'M':13, 'F':14, 'P':15, 'S':16, 'T':17, 'W':18, 'Y':19, 'V':20}
        X = []
        xName = []
        
        for i in fin:
            tokens = i.strip().split(',')
            X.append(tokens[2])
            xName.append(tokens[0])
        X = np.array(X)
        avgLen = sum([len(i) for  i in X])/(len(X)*1.0)
        maxL = max([len(i) for  i in X])
        diOrder = {}

        nucL =  ['A','R','N','D','C','E','Q','G','H', 'I','L','K','M', 'F', 'P','S','T','W', 'Y', 'V','END']
        table = np.zeros((len(nucL)**2,maxL))
        count = 0
        AAL =  {'STA':0 ,'A':1, 'R':2,'N':3, 'D':4,  'C':5,  'E':6,  'Q':7,  'G':8,  'H':9,  'I':10, 'L':11, 'K':12, 'M':13, 'F':14, 'P':15, 'S':16, 'T':17, 'W':18, 'Y':19, 'V':20,'END':21}
        tableEm = np.zeros((len(AAL),len(AAL),kmerLen))
        tableTau = np.zeros((len(AAL),len(AAL)))
        psTua = 0
        psEm = np.zeros(kmerLen)
        for i in nucL:
                for j in nucL:
                    diOrder[i+j] = count
                    count = count + 1

        makeHMM(hmmSeq)
        
        count = 1
        predTable = list(range(len(X)))   
        inf = open('emission' + name  +'.txt')
                
        count = count + 1
        tableEm = np.zeros((len(AAL),len(AAL),kmerLen))
        posi = 0
        posj = 0
        nzCount = 0
        for line in inf:
                tableEm[posi,posj] = (line.strip().split(','))[1:]
                if np.mean(tableEm[posi,posj]) != 0:
                        nzCount =  nzCount+1
                        psEm+=tableEm[posi,posj]
                posj += 1
                if posj == len(AAL):
                        posj = 0
                        posi +=1
        psEm = psEm/nzCount
        for j in range(len(X)):
                print("Predicting 8-mer binding profile for sequence: " + str(j+1))
             
                testingKmers = []
                for k in range(32896):
                        pred = walkScore(X[j],kmer8[k])
                        testingKmers.append(pred)
                predTable[j] = testingKmers
     
        np.savetxt('PREDICTIONS' + name +'.csv', np.transpose(np.array(predTable)), delimiter=',')
                
