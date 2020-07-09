#!/usr/bin/python3
import random
import numpy as np
import collections
import scipy.stats
global kmer5
global kmer8
global profiles
global X
global table
global AAL
global tableEm 
global avgLen
import time


def predict(trainSet,outFile):
        conv = False
        kmerSize = .1
        eta=0.0001
        mo = 0.99
        e = .05
        count = 0
        score = 0
        out2 = open('emission'+ outFile +'.txt', 'w')
        out3 = open('error'+ outFile +'.txt', 'w')
        sumError = 0
        
        #print(score)
        go = 10
        plat = False
        while not conv:
                sample = random.sample(trainSet, int(len(trainSet)))
                for i in sample:
                        plat = False
                        if  not plat:
                                testKmers = random.sample(range(len(profiles[0])), int((len(profiles[0]))*kmerSize))
                        else:
                                
                                testKmers = wKmerList[i]
                        plat = False
                        tempKmers = profiles[i]
                        for j in testKmers:
                                pred = walkScore(X[i],kmer8[j])
                                error = tempKmers[j] - pred
                                walkLearn(X[i],(np.array(np.zeros(len(kmer8[j]))+(error*eta))),kmer8[j])
                                out3.write(str(error) + '\n')
                                count = count + 1
                sumError = 0
                      
                sampleTest = random.sample(trainSet, int(len(trainSet)*1))
                
                if go == 10:
                        
                        sumDiff = 0
                        avegDiff = np.zeros(len(tempKmers))
                        wKmerList = range(len(X))
                        for i in sampleTest:
                                tempKmers = profiles[i]
                                testingKmers = []
                                error = 0.0
                                diff = 0.0
                                diffList = []
                                for j in range(len(tempKmers)):
                                        pred = walkScore(X[i],kmer8[j])
                                        testingKmers.append(pred)
                                        diff += abs(tempKmers[j]-pred)
                                        diffList.append(abs(tempKmers[j]-pred))
                                sumDiff += diff/len(tempKmers)
                                sumError += scipy.stats.spearmanr(tempKmers,testingKmers)[0]
                        AvegScore = sumError/len(sampleTest)
                        print("Done 10 epochs")
                        print(score,AvegScore,sumDiff/len(sampleTest))
                        avegDiff = scipy.stats.stats.rankdata(avegDiff)
                        plat = True
                        if   abs(abs(score) - abs(AvegScore)) < e:
                                
                                
                                eta = eta*mo
                                convC = convC + 1
                                if convC == 1:
                                        conv = True
                                        break
                        else:
                                convC = 0
                        go = 1
                else:
                        go = go + 1
                score = AvegScore

        inv_AAL = {v: k for k, v in AAL.items()}
        
        print("SGD done")
        for i in range(len(inv_AAL)):
                for j in range(len(inv_AAL)):
                        out2.write(str(inv_AAL[i]) + str(inv_AAL[j]) + ',' + ','.join(map(str,tableEm[i,j,:])) + '\n')


def makeHMM(outFile,fileName):
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
    out = open(outFile + '.csv', 'w' )
    for i in range(len(table[:,0])):
        out.write(str(inv_diOrder[i]) + ',' + ','.join(map(str,table[i,:])) + '\n')
           
    print('Construction of DBD Sequence graph done')

def read5mers(fileName):
	file = open(fileName, 'r')
	for i in file:
		token = i.strip().split(',')
		kmer5[token[0]] = token[1]

def read8mers(fileName):
	file = open(fileName, 'r')
	for i in file:
		token = i.strip().split(',')
		kmer8.append([int(j) for j in token[1:]])
	

def makeEmTable():
        for i in range(len(AAL)):
                for j in range(len(AAL)):
                        tableEm[i,j] = np.zeros(kmerLen)


def walkScore(seq,mapK): 
        old = 'S'
        trans = 0
        total = np.zeros(3)
        for i in range(len(seq)-1):
                total += table[diOrder[seq[i:i+2]],i] * tableEm[AAL[seq[i]],AAL[seq[i+1]],mapK]
        else:
                
                total += table[diOrder[seq[i+1] + 'END'],i] * tableEm[AAL[seq[i]],AAL[seq[i+1]],mapK]

        return sum(total)

def walkLearn(seq,error,mapK):
        tableEm[AAL['STA'],AAL[seq[0]],mapK] += error
        for i in range(len(seq)-1):
 
                tableEm[AAL[seq[i]],AAL[seq[i+1]],mapK] += error
        else:
                tableEm[AAL[seq[i]],AAL['END'],mapK] += error
                    


        
if __name__ == "__main__":
        SeqFile = 'seqNoTBOX.csv' #Sequence file
        HMMFile = 'mouseSeqs.csv' #HMM file
        profileFile = 'profileNoTBOX.csv' #Profile
        TrainSplit = 10
        kmerLen = 4**6
        kmer5 = {}
        kmer8 = []
        read8mers('Full6merComb.csv') #8 mer comb list
        outFile = 'MouseTest' #name of the output file
        print(outFile)
        table = []        
        profiles = []
        Nprofiles = []
        file = open(profileFile, 'r')
        first = True
        for row in file:
                tokens = row.strip().split(',')
                for i in range(len(tokens)):
                        if first:
                                profiles.append([])
                                Nprofiles.append(tokens[i])
                        else:
                                profiles[i].append(float(tokens[i]))
                first = False
        profiles = np.array(profiles)

        file = open(SeqFile, 'r')
        fin = [row.strip() for row in file ]
        X = []
        xName = []
        for i in fin:
                tokens = i.strip().split(',')
                X.append(tokens[2])
                xName.append(tokens[0])
        X = np.array(X)
        maxL = max([len(i) for  i in X])
    
        diOrder = {}
        nucL =  ['A','R','N','D','C','E','Q','G','H', 'I','L','K','M', 'F', 'P','S','T','W', 'Y', 'V','END']
        table = np.zeros((len(nucL)**2,maxL))
        count = 0
        for i in nucL:
                for j in nucL:
                    diOrder[i+j] = count
                    count = count + 1

 
        avgLen = sum([len(i) for  i in X])/(len(X)*1.0)
        print("Average sequence length " + str(avgLen))
    
        AAL =  {'STA':0 ,'A':1, 'R':2,'N':3, 'D':4,  'C':5,  'E':6,  'Q':7,  'G':8,  'H':9,  'I':10, 'L':11, 'K':12, 'M':13, 'F':14, 'P':15, 'S':16, 'T':17, 'W':18, 'Y':19, 'V':20,'END':21}
        tableEm = np.zeros((len(AAL),len(AAL),kmerLen))

        sample = list(range(len(xName)))
        random.shuffle(sample)
        HMM = makeHMM(outFile,HMMFile)
        makeEmTable()    
 
        predict(sample,outFile)
       
