#! /usr/bin/env python
import sys
import argparse
import screed
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg') 
from pylab import *
import os.path
from scipy import stats

def getstats(infile):
    slope = [] 
    intercept = [] 
    RMSE = [] 
    r2 = []
    row = 0
    sample =infile.split(".all")[0]
    print(sample)
    for line in open(infile):
            x = line.rstrip().split(' ')
            s1 = round(float(x[0]), 4)  
            if row ==0: 
                realslope = float(x[0])
                print('Real Slope is ', realslope) 
            else:
                slope.append( float(s1) )
            row+=1
    count = 0 
    for i in slope: 
        if i <= realslope: 
            count += 1
    print('Count is', count, 'N of slopes is ', len(slope), 'and PV is ', count/len(slope))

    plt.hist(slope, 20,
         density=False,
         histtype='bar',
         facecolor='teal',
         alpha=0.5)
    plt.title('Histogram of Slopes of Shuffled Footprints of Sample '+str(sample), fontsize=10)
    xtext = 'Real slope = '+ str(round(realslope,4)) +' Count ='+str(count)+' No. of slopes = ' +str(len(slope)) +' PValue = '+str(count/len(slope)) 
    plt.xlabel(str(xtext), fontsize=8)
    
    plt.savefig(str(sample)+".slopes.png", format='png')
    plt.clf()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', default=False) 
    args = parser.parse_args()
    getstats(args.infile) 
if __name__ == '__main__':
    main()

   
