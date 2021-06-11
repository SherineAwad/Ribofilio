#! /usr/bin/env python
import sys
import argparse
import screed
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os.path
from scipy.stats import t

def getstats(infile1, infile2):
    
    sample1 = [] 
    sample2 = [] 
    count = 0
    for line in open(infile1):
        if count == 0: 
            count+=1 
            continue
        slope1, _, _, _, SE1, _, _, _, n1 = line.split('\t')
    count = 0
    for line in open(infile2): 
        if count == 0:
             count+=1 
             continue
        slope2, _, _, _, SE2, _, _, _, n2 = line.split('\t') 
    df = ( float(n1) + float(n2) )-4
    tscore = (float(slope1) - float(slope2)) / np.sqrt(np.square(float(SE1)) +np.square(float(SE2)) ) 
    pvalue =  (t.sf(abs(tscore),df= df))
    pvalue = np.round(pvalue, decimals=4)
    print("tscore and pvalues are", tscore, pvalue)
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile1', default=False) 
    parser.add_argument('infile2', default=False)
    args = parser.parse_args()
    getstats(args.infile1, args.infile2) 

if __name__ == '__main__':
    main()

   
