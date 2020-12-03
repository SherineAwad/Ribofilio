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
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score

#------------------------------------------------------------------------
#profile function: estimates the drop rate of ribosomes after binning
#------------------------------------------------------------------------
def profile(sample, BINSIZE, coverage, gLength,posmax, posmin,cutoff):
    
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
            if (j%1000000 ==0):
                 print (".....................",j)
    print('Done filling the counts dictionary')

    N = int (posmax / BINSIZE) + 1
    gbins = [0] * N
    npos = [0] * posmax
    bCovered = [0] * N
    print("Posmax is ", posmax)
 
    #Counts genes in each position 
    for gene in coverage:
        for i in coverage[gene]:
            if gLength[gene] <= posmax and gLength[gene] > posmin:
                pos[int(i)] +=1
    print("Filling pos is done")

    #Count how many genes cover each position, so if a position is coverved by x and x <cutoff, then this position is discarded as noise 
    for gene in coverage:
        if gLength[gene] <= posmax and gLength[gene] > posmin:
            for i in range(0, gLength[gene]+1):
                gCovered[i] +=1

    print("Filling gCoverered is done")
    i = 0
    last_pos = 0
    #Normalize bin position  with the number of gene covering that position
    #Discard bins with number of genes less than cutoff and record the position in last_pos 
    while(i < posmax):
        if (gCovered[i] >= cutoff) :
            i +=1
            continue  
        else:     
            npos[int(i)] = pos[int(i)] / (gCovered[int(i)] + c)
            last_pos = i
            i +=1
    print('last_pos is ', last_pos)     
    print("Normalizing pos is done")
    #Here, we group reads into bins and normalize by BINSIZE, we stop at last_pos
    index = 0
    a = 0
    b = BINSIZE
    gcounts = 0
    while a< posmax:
        for i in range (a,b+1):
            if i >  (len(npos) - 1):
                break
            if i > last_pos: 
                break
            posum = float(npos[i])
            gbins[index] +=posum
            bCovered[index] =gCovered[i]
        gbins[index] = (float(c +  (gbins[index]/BINSIZE) ) )
        index +=1
        a = b+1
        b = b +BINSIZE
    return gbins, last_pos, bCovered 

#--------------------------------------------------------------------------
#----------run_subset function
#--------------------------------------------------------------------------

#This function finds runs the ribosome profile on a subset of genes
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
                if (j%1000000 ==0):
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

#-------------------------------------------------------------------------------------------
#plots function, plot several figures for the ribosome profiling 
#-------------------------------------------------------------------------------------------

def plots(plot_flag,outName,N,coverage,all_bins,binsize, posmax, posmin,last_pos,gLength,ymin,ymax,ylogmin,ylogmax,cutoff,end,bCovered,regfile):
    print('Start plotting')  
    BINS = []
    LOGBINS = []
    pgLength = []
    gbins = []
    outName = outName 
    genesfile = outName +".genes" 
    for gene in coverage: 
        if gene in gLength:
            if gLength[gene] <= posmax and gLength[gene] > posmin:
                pgLength.append(int(gLength[gene]))
    print('N is: ', N) 
    if end !=-1: 
        N = N - end
    print('N now is ', N)
    #Take log of bins indeces for LOG/LOG plotting
    for i in range (0,N) :
        BINS.append(i)
        gbins.append(all_bins[i])
    LOGBINS.append(0)
    i = 0 
    for i in range (1,N):
        LOGBINS.append(np.log(i) )
    print('Prining BINS, gbins and LOGBINS:')
    print('lengths of bins, logbins, and gbins is:', len(BINS), len(LOGBINS), len(gbins))
    label = "Bins  (Binsize = " + str(binsize) + " Cutoff = "+ str(cutoff) +" Max gene Length ="+str(last_pos)+" End= "+str(end)+")"

    if plot_flag ==1: 
        #Genes Length Histogram 
        #-------------------------
        plt.hist(pgLength, 100,density=False,histtype='bar',facecolor='Navy',alpha=0.5)
        plt.xlabel('Genes Length')
        plt.ylabel('Frequency')
        plt.title(outName+' genes length distribution', fontsize=12)
        plt.savefig(outName+".Length.Histogram.png", format='png')
        plt.clf()

        #Gene bins Histogram 
        #-----------------------
        plt.hist(bCovered, 100,density=False,histtype='bar',facecolor='Navy',alpha=0.5)
        plt.xlabel('No. of Genes')
        plt.ylabel('Bin') 
        plt.title(outName+'Genes covered in each bin distribution', fontsize =12) 
        plt.savefig(outName+".Bcovered.Histogram.png", format='png')
        plt.clf()
    
        #Plot Coverage without Log
        #----------------------------
        plt.figure(figsize=(15,10))
        plt.plot(np.array(BINS),np.array(gbins),color='Navy')
        plt.ylim(ymin,ymax)
        plt.ylabel('Sums of coverage of genes at each position in each bin')
        plt.xlabel(label)
        plt.title(outName+' Coverage per Bin (Sum without log) ', fontsize=12)
        plt.savefig(outName+"_"+str(binsize)+"_"+str(cutoff)+"_"+str(end)+".NoLog.png", format='png')
        plt.clf()
    
        #Plot Coverage Log/Log 
        #--------------------------
        plt.figure(figsize=(15,10))
        plt.plot(np.array(LOGBINS),np.array(np.log(gbins)),color='Navy')
        plt.ylabel('Log of sums of coverage of genes at each position in each bin')
        plt.ylim(ylogmin,ylogmax)
        plt.xlabel(label)
        plt.title(outName+' Coverage per Bin (Log/Log) ', fontsize=12)
        plt.savefig(outName+"_"+str(binsize)+"_"+str(cutoff)+"_"+str(end)+".LogLog.png", format='png')
        plt.clf()
        #Plot coverage Log /Linear 
        #-----------------------------
        plt.figure(figsize=(15,10))
        plt.plot(np.array(BINS),np.log(gbins),color='Navy')
        plt.ylim(ylogmin, ylogmax)
        plt.ylabel('Log of sums of coverage of genes at each position in each bin')
        plt.xlabel(label)
        plt.title(outName+' Coverage per Bin (Log/Linear) ', fontsize=12)
        plt.savefig(outName+"_"+str(binsize)+"_"+str(cutoff)+"_"+str(end)+".LogLinear.png", format='png')
        plt.clf()
        #print('Printing bins in plotting', gbins) 
    #Linear Regression 
    #------------------
    x = np.array(BINS).reshape(-1, 1)
    y = np.array(np.log(gbins)).reshape(-1, 1)
    regression_model = LinearRegression()
    W = binsize
    weight = BINS
    for i in range(0, len(weight)):
        #weight[i] = 1/bCovered[i] #5002
        weight[i] = bCovered[i] /(i+1) #5003
        #weight[i] = len(weight) - i
        #weight[i] = abs(BINS[i] - W) 
        #weight[i] = BINS[i] * ( len(BINS) - BINS[i])
        #weight[i] = gbins[i]/(BINS[i] + c) 
        #weight[i] = (BINS[i]+c)/gbins[i]
    #print('Weight is ')
    #print(weight)

    # Fit the data(train the model)
    regression_model.fit(x, y,weight)
    # Predict
    y_predicted = regression_model.predict(x)

    # model evaluation

    rmse = mean_squared_error(y, y_predicted,weight)
    r2 = r2_score(y, y_predicted,weight)

    # printing values
    print('Slope:' ,regression_model.coef_)
    print('Intercept:', regression_model.intercept_)
    print('Root mean squared error: ', rmse)
    print('R2 score: ', r2)
    plt.figure(figsize=(15,10))
    plt.ylim(ylogmin,ylogmax) 
    plt.scatter(x, y, s=10)
    xtext = ' Slope: '+str(regression_model.coef_)+' Intercept: '+str(regression_model.intercept_)+' RMSE: '+str(rmse)+' R2 score :'+ str(r2) 
    plt.xlabel(str(label)+'\n'+str(xtext), fontsize=10)
    plt.ylabel('y')
    # predicted values
    plt.plot(x, y_predicted, color='r')
    plt.title(outName+' Linear Regression', fontsize=12)
    plt.savefig(outName+"_"+str(binsize)+"_"+str(cutoff)+"_"+str(end)+".WLR.png", format='png')
    regfp=open(str(regfile),'a+')
    print(str(regression_model.coef_[0][0]),str(regression_model.intercept_[0]),str(rmse),str(r2), file=regfp)
    plt.clf()

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
    parser.add_argument('-c','--cutoff',dest='cutoff', type=int,default=1000000000)
    parser.add_argument('--plots', dest='plots',type=int, default=1)
    parser.add_argument('-n', '--name', dest='regfile_name', default='regfile.log') 
    parser.add_argument('-e', '--end', dest='end', type=int,default=0)
    parser.add_argument('--ymin' ,dest='ymin',type=int, default=0) 
    parser.add_argument('--ymax',dest='ymax', type=int,default=30) 
    parser.add_argument('--ylogmin', type=int, default=-15)
    parser.add_argument('--ylogmax', type=int, default=5) 
    args = parser.parse_args()

    
    sample = str(args.footprint).split('.')[0] 
    if args.rna != "NULL": 
        sample+=  '_' + str(args.rna).split('.')[0]

    subset = str(args.subset)

    binsize = int(args.binsize)
        
    cutoff = int(args.cutoff)

    outName = sample   
    if subset != "NULL": 
        subset_name = subset.split('.')[0]
        outName += '.'+subset_name 
    plot_flag = args.plots
    print('plot_flag is', plot_flag)
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
    if args.rna !="NULL":
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
            gbins1, last_pos, bCovered = profile(args.footprint,binsize, coverage, gLength,posmax, posmin,cutoff)
            gbins2,_,_ = profile(args.rna,binsize, coverage, gLength,posmax, posmin,cutoff) 
        else:
            gbins1, posmax, gLength = run_subset(args.transcripts,args.footprint, subset, binsize)
            gbins2, posmax, gLength = run_subset(args.transcripts,args.rna, subset, binsize)
            last_pos = posmax
        print ('Length of gbins1 is ', len(gbins1) ) 
        print ('Length of gbins2 is ', len(gbins2) ) 
    else:
        if subset =="NULL":
            gbins1,last_pos, bCovered = profile(args.footprint,binsize, coverage, gLength,posmax, posmin,cutoff)
        else: 
            gbins1, posmax, gLength = run_subset(args.transcripts,args.footprint, subset, binsize)  
            last_pos = posmax
    
    N = int (last_pos/binsize) + 1 
    BINS = []
    gbins = gbins1
    if args.rna !="NULL": 
        for i in range(0, len(gbins1) ):
            gbins[i] = float(gbins1[i] / gbins2[i])
    #-----------------------------------------
    #Plotting and summarizing  
    #-----------------------------------------
    plots(plot_flag,outName,N,coverage,gbins,binsize, posmax, posmin,last_pos,gLength,int(args.ymin),int(args.ymax),int(args.ylogmin),int(args.ylogmax),cutoff,int(args.end),bCovered,args.regfile_name )
    print('Done Plotting and summarizing')
    print("All is done") 
if __name__ == '__main__':
    main()
