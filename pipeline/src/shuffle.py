#! /usr/bin/env python
import sys
import argparse
import screed
import numpy as np

#------------------------------------------------------------------------
#profile function: estimates the drop rate of ribosomes after binning
#------------------------------------------------------------------------
def shuffle(transcripts,sample,outfile,N):
 
    gLength = {} 

    for record in screed.open(transcripts):
        gname = record.name.split(' ')[0]
        if "mRNA" in gname:
            gname=gname.split('_')[0]
        gLength[gname] = len(record.sequence)

    files = [] 
    for i in range(0,N): 
        fname = outfile +"."+str(i)+".bed" 
        files.append(fname)
    for f in files: 
       fp = open(f, 'w+')
       for line in open(sample):
            l = line.rstrip().split('\t')
            gname = l[0].strip()
            if "_mRNA" in gname : 
                gname = str(l[0]).split('_')[0] 
            read, pos1, pos2 = l[0] , int(l[1]), int(l[2])
            read_size = pos2 - pos1
            randy = np.random.randint(0, gLength[gname]-read_size, size=1) 
            rpos1 = randy[0]
            rpos2 = rpos1 + read_size 
            print (l[0],'\t',str(rpos1),'\t',str(rpos2),'\t',l[3],'\t',l[4],'\t',l[5],file=fp )
       fp.close() 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t','--transcripts',dest='transcripts', default=False)
    parser.add_argument('-b','--bed', dest='bed', default=False)
    parser.add_argument('-o', '--outfile', dest='outfile', default=False) 
    parser.add_argument('-n', dest='N', type= int, default=False) 
    args = parser.parse_args()

    shuffle(args.transcripts, args.bed,args.outfile, args.N) 
if __name__ == '__main__':
    main()

   
