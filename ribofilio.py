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

#------------------------------------------------------------------------
#profile function: estimates the drop rate of ribosomes after binning
#------------------------------------------------------------------------
def profile(sample, BINSIZE, coverage, gLength,posmax, posmin):
    
    gbins = [] 
    c = 0.000001
    gCovered = [0] * 100000000
    pos = [0] * 100000000
    j = 0
    for line in open(sample):
            x = line.rstrip().split('\t')
            gname = x[0].strip()
            if "_mRNA" in gname : 
                gname = str(x[0]).split('_')[0] 
            coverage[str(gname)].append(int(x[2])) 
            pos[int(x[2])] =0
            j +=1
            if (j%100000 ==0):
                 print (".....................",j)
    print('Done filling the counts dictionary')

    N = int (posmax / BINSIZE) + 1
    gbins = [0] * N
    npos = [0] * posmax
    print("Posmax is ", posmax)

    for gene in coverage:
        for i in coverage[gene]:
            if gLength[gene] <= posmax and gLength[gene] > posmin:
                pos[int(i)] +=1
    print("Filling pos is done")

    for gene in coverage:
        if gLength[gene] <= posmax and gLength[gene] > posmin:
            for i in range(0, gLength[gene]+1):
                gCovered[i] +=1
    print("Filling gCoverered is done")

    for i in range(0, posmax):
        npos[int(i)] = pos[int(i)] / (gCovered[int(i)] + c)
    print("Normalizing pos is done")



    index = 0
    a = 0
    b = BINSIZE
    gcounts = 0
    while a< posmax:
        for i in range (a,b+1):
            if i >  (len(npos) - 1):
                break
            posum = float(npos[i])
            gbins[index] +=posum
        gbins[index] = (float(c +  (gbins[index]/BINSIZE) ) )
        index +=1
        a = b+1
        b = b +BINSIZE
    return gbins 
 
#--------------------------------------------------------------------------
#run_sunset function
#This function finds runs the ribosome profile on a subset of genes
#--------------------------------------------------------------------------
def run_subset(transcripts,sample, subset_file, binsize): 
    
    subset = []
    coverage = {}
    gLength = {}
    gmax = -1
    posmax = 0
    posmin = 0
    gmin = 0
    c = 0.000001
    gCovered = [0] * 100000000
    pos = [0] * 100000000

    with open(subset_file) as f:
        gener = f.read().splitlines()
    for i in gener:
        if "_mRNA" in i:
            i = i.split('_')[0]
        subset.append(str(i).strip() )
    print("Starting reading input files")
    for record in screed.open(transcripts):
        gname = record.name.split(' ')[0]
        if "_mRNA" in gname: 
            gname =gname.split('_')[0] 
        if gname in subset:
            coverage[gname] = []
            gLength[gname] = len(record.sequence)
            if gLength[gname] > gmax:
                    gmax = gLength[gname] 
    print('Done reading Transcript')

    posmax = gmax
    posmin = gmin

    print("gmin, gmax, posmin, and posmax are:", gmin,gmax,posmin,posmax)

    j = 0
    for line in open(sample):
            x = line.rstrip().split('\t')
            gname =str(x[0])
            gname = x[0].strip()
            if "_mRNA" in gname :
                gname = str(x[0]).split('_')[0]
            if gname in subset:
                coverage[gname].append(x[2])
                pos[int(x[2])] =0
                j +=1
                if (j%100000 ==0):
                    print (".....................",j)
    print('Done filling the counts dictionary')

    N = int (posmax / binsize) + 1
    print ("Number of bins is ", N)

    gbins = [0] * N
    npos = [0] * posmax
    print("Filling gCovered is done")
    
    for gene in coverage:
        for i in coverage[gene]:
            if gLength[gene] <= posmax and gLength[gene] > posmin:
                pos[int(i)] +=1
    print("Filling pos is done")

    for gene in coverage:
        if gLength[gene] <= posmax and gLength[gene] > posmin:
            for i in range(0, gLength[gene]+1):
                gCovered[i] +=1
    print("Filling gCoverered is done")
    
    for i in range(0, posmax):
        npos[int(i)] = pos[int(i)] / (gCovered[int(i)] + c)
    print("Normalizing pos is done")


    index = 0
    a = 0
    b = binsize 
    gcounts = 0
    while a< posmax:
        for i in range (a,b+1):
            if i >  (len(npos) - 1):
                break
            posum = float(npos[i])
            gbins[index] +=posum
        gbins[index] = (float(c +  (gbins[index]/binsize) ) )
        index +=1
        a = b+1
        b = b +binsize
    print("Binning is done")
    print("N and index, gbins[0],and gbins[index-1] are:", N,index,gbins[0],gbins[index-1])
    return(gbins,posmax,gLength) 
 
#---------------------------------------------------------------------------
#-verbose function, prints summary of running process
#---------------------------------------------------------------------------
def verbose(outName,coverage, gbins, posmax,posmin,N,BINSIZE,pgLength):
    gbinsfile = outName+".gbins"
    fgbins = open(gbinsfile, 'w+')
    for index in range(0, N):
        print(gbins[index], file=fgbins)
    fgbins.close()
    
    coveragefile = outName+ ".coverage"
    fcoverage = open(coveragefile, 'w+')
    for gene in coverage:
        fcoverage.write(gene +",")
        for i in coverage[gene]:
            fcoverage.write(str(i)+",")
        fcoverage.write("\n")
    fcoverage.close()


    summaryfile = outName + '.summary.txt'
    fsummary = open(summaryfile, 'w+')
    print('Summary of run on sample', sample, 'with binsize is', BINSIZE, file=fsummary)
    print('In this run we considered genes with length between', posmin, 'and ', posmax, file=fsummary)
    print('The maximum gene length in this run is ', posmax, file=fsummary)
    print('Number of bins is ', N)
    print('The value of the first bin is ', gbins[0],file=fsummary )
    print('The value of the last bin is' , gbins[index-1],file=fsummary )
    print('Number of genes in this run is', len(pgLength), file=fsummary)
    fsummary.close()
        
#-------------------------------------------------------------------------------------------
#plots function, plot several figures for the ribosome profiling 
#-------------------------------------------------------------------------------------------
def plots(outName,N,coverage,gbins,binsize, posmax, posmin,gLength,ymin,ymax,ylogmin,ylogmax):
    print('Start plotting')  
    BINS = []
    LOGBINS = []
    pgLength = []
    genesfile = outName +".genes"
    for gene in coverage: 
        if gene in gLength:
            if gLength[gene] <= posmax and gLength[gene] > posmin:
                pgLength.append(int(gLength[gene]))
    
    print('N is in main', N) 
    
    #Take log of bins indeces for LOG/LOG plotting
    for i in range (0,N) :
        BINS.append(i)
    LOGBINS.append(0)
    i = 0 
    for i in range (1,N):
        LOGBINS.append(np.log(i) )
    print('Prining BINS, gbins and LOGBINS:')
    print('lengths of bins, logbins, and gbins is:', len(BINS), len(LOGBINS), len(gbins))
    #Genes Length Histogram 
    #-------------------------
    plt.hist(pgLength, 100,density=False,histtype='bar',facecolor='Navy',alpha=0.5)
    plt.xlabel('Genes Length')
    plt.ylabel('Frequency')
    plt.title('Histogram of '+outName+ ' genome length distribution')
    plt.savefig(outName+".Length.Histogram.png", format='png')
    plt.clf()

    #Plot Coverage without Log
    #----------------------------
    label = "Bins (Binsize = " + str(binsize) + ")"
    plt.plot(np.array(BINS),np.array(gbins),color='Navy')
    plt.ylim(ymin,ymax)
    plt.ylabel('Sums of coverage of genes at each position in each bin')
    plt.xlabel(label)
    plt.title(outName+' Coverage per Bin (Sum without log) ')
    plt.savefig(outName+".NoLog.png", format='png')
    plt.clf()

    #Plot Coverage Log/Log 
    #--------------------------
    plt.plot(np.array(LOGBINS),np.array(np.log(gbins)),color='Navy')
    plt.ylabel('Log of sums of coverage of genes at each position in each bin')
    plt.ylim(ylogmin,ylogmax)
    plt.xlabel(label)
    plt.title(outName+' Coverage per Bin (Log/Log) ')
    plt.savefig(outName+".LogLog.png", format='png')
    plt.clf()
    #Plot coverage Log /Linear 
    #-----------------------------
    plt.plot(np.array(BINS),np.log(gbins),color='Navy')
    plt.ylim(ylogmin, ylogmax)
    plt.ylabel('Log of sums of coverage of genes at each position in each bin')
    plt.xlabel(label)
    plt.title(outName+' Coverage per Bin (Log/Linear) ')
    plt.savefig(outName+".LogLinear.png", format='png')
    plt.clf()
    
    #Linear Regression
    #------------------------
    x = np.array(BINS)
    y = np.array(gbins) #x is just the bin
    m,b = np.polyfit(x, y, 1)
    plt.plot(x, y, 'o', x, m*x+b, '--k',color='Navy')
    plt.ylim(ymin, ymax)
    plt.title(outName+' :x is bin index based Regression')
    plt.savefig(outName+".regression.png", format='png')
    plt.clf()
    #Polynomial Regression
    #------------------------
    x = np.array(BINS)
    y = np.array(gbins)
    pol_coeff = np.polyfit(x,y,3)
    yfit = np.poly1d(pol_coeff)
    xnew = np.linspace(x[0],x[-1],100)
    plt.plot(xnew,yfit(xnew), x, y, 'o',color='Navy')
    plt.ylim(ymin, ymax)
    plt.title(outName+' :x is bin index based Polynomial Regression')
    plt.savefig(outName+".Pregression.png", format='png')
    plt.clf()
    verbose(outName,coverage, gbins, posmax,posmin,N,binsize,pgLength)

#---------------------------------------------------------------------------
#Main function: gets input paramters and calls corresponding functions 
#---------------------------------------------------------------------------

def main(): 
    print("Initializing and reading arguments") 
    parser = argparse.ArgumentParser()
    parser.add_argument('-t','--transcripts',dest='transcripts', default=False) 
    parser.add_argument('-f','--footprint', dest='footprint',default=False)
    parser.add_argument('-r','--rna', dest='rna',default="NULL")
    parser.add_argument('-s','--subset',dest='subset', default ="NULL")
    parser.add_argument('-b','--binsize',dest='binsize',type=int, default=50)
    parser.add_argument('--ymin' ,dest='ymin',type=int, default=0) 
    parser.add_argument('--ymax',dest='ymax', type=int,default=30) 
    parser.add_argument('--ylogmin', type=int, default=-5)
    parser.add_argument('--ylogmax', type=int, default=5) 
   
    args = parser.parse_args()

    
    sample = str(args.footprint).split('.')[0] 
    if args.rna != "NULL": 
        sample+=  '_' + str(args.rna).split('.')[0]

    subset = str(args.subset)

    binsize = int(args.binsize)
  
    outName = sample+'.'+str(args.binsize)  
    if subset != "NULL": 
        subset_name = subset.split('.')[0]
        outName += '.'+subset_name 

    print('outName is', outName) 
    coverage = {}
    gLength = {}  
    gene_snps = {} 
    gmax = -1
    posmax = 0
    posmin = 0
    gmin = 0 
    
    #Check required files exist 
    if not os.path.exists(args.transcripts): 
        sys.exit('transcripts file is required and does not exist, exiting !')
    if not os.path.exists(args.footprint):
        sys.exit('footprint file is required and does not exist, exiting !')
    if args.rna !="Null":
        if not os.path.exists(args.rna):
            sys.exit('rna  file does not exist: either run with footprint only option or enter a valid rna file, exiting !')
    
    #Read all the genes in the reference 
    print("Starting reading input files") 
    for record in screed.open(args.transcripts):
        gname = record.name.split(' ')[0]
        if "mRNA" in gname: 
            gname=gname.split('_')[0] 
        coverage[gname] = [] 
        gene_snps[gname] = 0
        gLength[gname] = len(record.sequence)
        if gLength[gname] > gmax: 
            gmax = gLength[gname] 
    gmid = float (gmin + gmax) / 2
    print('Done reading Transcript') 
   
    posmax = gmax 
    posmin = gmin 
   
    print("gmin, gmax, posmin, and posmax are:", gmin, gmax,posmin,posmax) 
    print('We are calling profile for' , args.footprint)
    if args.rna != "NULL": 
        if subset == "NULL":
            gbins1 = profile(args.footprint,binsize, coverage, gLength,posmax, posmin)
            gbins2 = profile(args.rna,binsize, coverage, gLength,posmax, posmin) 
        else: 
            gbins1, posmax, gLength = run_subset(args.transcripts,args.footprint, subset, binsize)
            gbins2, posmax, gLenght = run_subset(args.transcripts,args.rna, subset, binsize)
        print ('Length of gbins1 is ', len(gbins1) ) 
        print ('Lenght of gbins2 is ', len(gbins2) ) 
    else: 
        if subset =="NULL":
            gbins1 = profile(args.footprint,binsize, coverage, gLength,posmax, posmin)
        else: 
            gbins1, posmax, gLength = run_subset(args.transcripts,args.footprint, subset, binsize)  


    N = int (posmax / binsize) + 1


    gbins = gbins1
    if args.rna !="NULL": 
        for i in range(0, len(gbins1) ):
            gbins[i] = float(gbins1[i] / gbins2[i])
        print('Done calculating gbins') 
    #-----------------------------------------
    #Plotting and summarizing  
    #-----------------------------------------
    plots(outName,N,coverage,gbins,binsize, posmax, posmin,gLength,int(args.ymin),int(args.ymax),int(args.ylogmin),int(args.ylogmax) )
    print('Done Plotting and summarazing')
    print("All is done") 
if __name__ == '__main__':
    main()
