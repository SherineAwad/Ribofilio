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

def is_sane(infile):
    flag =True 
    count = 1
    for line in open(infile):
            x = line.rstrip().split('\t')
            if len(x) != 6:  
                print('Bed file has to be 6 columns: check issues at line', count)
                flag = False
            count+=1 
    if flag == True: 
        print('Bed file seems ok' )

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile',dest='infile', default=False) 
    args = parser.parse_args()
    is_sane(args.infile) 
if __name__ == '__main__':
    main()

   
