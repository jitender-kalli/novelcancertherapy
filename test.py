import numpy as np
import pandas as pd
import os
import sys


def uniquepeptidecheck():
    df = pd.read_csv("lungpeptides.csv", sep=',')
    df2 = pd.read_csv("normalpeptides.csv", sep='\t',names=["peptide", "core", "icore", "score"]);
    #usecols=["peptide"];
    #df.columns=df["peptide"];
    lungpeptides=df["peptide"].tolist()
    normalpeptides=df2["peptide"].tolist()
    uniquepeptides = [];
    #cancerpeptides=['ALQPFPAPV', 'SQVEKLVRV', 'YLXXXXXLL'];
    ctr=0;
    for i in lungpeptides :
    	if i in normalpeptides:
    		a=10;
    	else:
    		ctr+=1;
    		score=0;
    		for j in normalpeptides:
    			for x in range(len(j)):
    				if i[x]== j[x]:
    					score=score+1;
    					print(score);
    					break;
    				break;	
    			break;	
    print (ctr);
	
uniquepeptidecheck()
