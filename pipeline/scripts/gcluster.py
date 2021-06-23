#! /usr/bin/env python
import sys
import argparse
import screed
import math
import numpy as np
import matplotlib.pyplot as plt

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('transcript')
    parser.add_argument('outfile')
    parser.add_argument('-x','--max',dest='max',type=int,default=10000000)
    parser.add_argument('-n','--min',dest='min', type=int,default=0)
    args = parser.parse_args()
    genome = {} 
    gLength = [] 
    
    fout1 =open('yeastGenes.txt', 'w+')

    c = 0.0000001
    for record in screed.open(args.transcript):
        name = record.name.split(' ')[0]
        sequence = record.sequence 
        genome[name] = sequence
        print(name, file=fout1)
    fout1.close()
   
    fout2 =open(args.outfile, 'w+') 
    for i in genome:
        gLength.append(len(genome[i]) )
        if len(genome[i]) <= int(args.max) and len(genome[i]) > int(args.min)   : 
            print(i, file=fout2)
    print('Done reading Transcripts and clustering genes')
    fout2.close()
    plt.hist(gLength, 100,
         density=False,
         histtype='bar',
         facecolor='Navy',
         alpha=0.5)
    plt.xlabel('Genes Length',fontsize=14)
    plt.ylabel('Frequency',fontsize=14)
    plt.title('Histogram of Saccharomyces cerevisiae genome length',fontsize=12)
    figure = plt.gcf()
    figure.set_size_inches(15, 15)
    plt.savefig(args.transcript+".png", format='png', dpi =100) 

if __name__ == '__main__':
    main()

   
