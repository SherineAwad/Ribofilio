#! /usr/bin/env python
import argparse
import screed
import math
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score

# -----------------------------------------------
# Read and parse input files
# -----------------------------------------------
def get_subset_genes(transcripts, subset_file):
    subset = []
    max_gene_length = -100 
    genes_length = {} 
    for line in open(subset_file):
        gene_name = line.rstrip()
        if "mRNA" in gene_name:
            gene_name=gene_name.split('_')[0]
        subset.append(gene_name) 
    for record in screed.open(transcripts):
        gene_name = record.name.split(' ')[0]
        if "mRNA" in gene_name:
            gene_name=gene_name.split('_')[0]
        if gene_name in subset:
            genes_length[str(gene_name)] = int(len(record.sequence))
            if genes_length[gene_name] > max_gene_length:
                max_gene_length = genes_length[gene_name]
    print ('Maximum gene length is',max_gene_length)
    return max_gene_length, genes_length

def get_genes(transcripts): 
    max_gene_length = -100 
    genes_length = {} 
    itr = 0 
    for record in screed.open(transcripts):
        gene_name = record.name.split(' ')[0]
        itr += 1 
        if (itr % 1000) ==0: 
            print(".....................")
        if "mRNA" in gene_name:
            gname=gname.split('_')[0]
        genes_length[str(gene_name)] = int(len(record.sequence))
        if genes_length[gene_name] > max_gene_length:
            max_gene_length = genes_length[gene_name]
    print ('max_gene_length is',max_gene_length)
    return max_gene_length, genes_length 


def get_reads(infile, genes_length):
    coverage = {} 
    itr = 0 
    for gene in genes_length: 
        coverage[gene] = [] 
    for  line in open(infile):
        record = line.rstrip().split("\t")
        itr += 1
        if (itr % 5000000) ==0:
            print(".....................")
        gene_name = str(record[0])
        if "mRNA" in gene_name:
            gene_name = gene_name.split("_")[0]
        if gene_name in genes_length: 
            coverage[str(gene_name)].append(int(record[2]) )
    return coverage 

def get_gene_coverage_at_bin(max_gene_length, binsize, genes_length):
    num_bins = int(max_gene_length / binsize) + 1
    gene_coverage_at_bin = [0] * num_bins
    # Fill gene_coverage_at_bin
    for gene in genes_length:
        bin_fit = math.ceil(genes_length[gene] / binsize)
        for i in range(0, bin_fit):
            gene_coverage_at_bin[i] += 1
    return gene_coverage_at_bin

def get_gene_coverage_at_pos(max_gene_length, coverage, genes_length) :
    gene_coverage_at_pos = [0] * (max_gene_length + 1)
    for gene in coverage:
            for i in range(1, genes_length[gene] + 1):
                gene_coverage_at_pos[i] += 1
    return gene_coverage_at_pos

def fill_positions (coverage, max_gene_length):
    positions = [0] * (max_gene_length + 1 ) 
    for gene in coverage:
        for i in coverage[gene]:
                positions[int(i)] += 1
    return positions

# ------------------------------------------------------------------------
# binning function: estimates the drop rate of ribosomes after binning
# ------------------------------------------------------------------------
def binning(
    binsize, positions, gene_coverage_at_pos, max_gene_length):
    gene_bins = []
    c = 0.000001
    
    num_bins = int(max_gene_length / binsize) + 1
    gene_bins = [0] * num_bins
    normalized_positions = [0] * (max_gene_length + 1)

    i = 1
    # Normalize bin position  with the number of gene covering that position
    while i <= max_gene_length:
        normalized_positions[int(i)] = positions[int(i)] / (gene_coverage_at_pos[int(i)] + c)
        i += 1
    
    # Here, we group reads into bins and normalize by binsize, we stop at last_pos
    index = 0
    a = 1
    b = binsize
    while a <= max_gene_length:
        for i in range(a, b + 1):
            if i > (len(normalized_positions) - 1):
                break
            position_sum = float(normalized_positions[i])
            gene_bins[index] += position_sum
        gene_bins[index] = c +float( gene_bins[index] / binsize)
        index += 1
        a = b + 1
        b = b + binsize
    return gene_bins

# -------------------------------------------------------------------------------------------
# Regression function, plot several figures for the ribosome profiling
# -------------------------------------------------------------------------------------------
def regression(
    output,
    num_bins,
    coverage,
    all_bins,
    binsize,
    max_gene_length,
    genes_length,
    ymin,
    ymax,
    ylogmin,
    ylogmax,
    gene_coverage_at_bin,plot 
):

    bins = []
    log_gene_bins = []
    gene_bins = []
    output = output
    # Take log of bins indeces for LOG/LOG plotting
    for i in range(0, num_bins):
        bins.append(i)
        gene_bins.append(all_bins[i])
    log_gene_bins.append(0)
    for i in range(1, num_bins):
        log_gene_bins.append(np.log(i))
    label = "Bin Number" 

    # Weighted Linear Regression
    # ------------------
    x = np.array(bins).reshape(-1, 1)
    y_in_log = np.array(np.log(gene_bins)).reshape(-1, 1)
    y_in_linear = np.array(gene_bins).reshape(-1, 1)
    regression_model = LinearRegression()
    weight = [0] * num_bins  
    norm_weight = [0] * num_bins
    for i in range(0, len(weight)):
        weight[i] = gene_coverage_at_bin[i]
        norm_weight[i] = cbrt(weight[i])
    # Fit the data(train the model)
    regression_model.fit(x, y_in_log, sample_weight = norm_weight)
    # Predict
    y_in_log_predicted = regression_model.predict(x)

    # model evaluation
    rmse_in_log = mean_squared_error(y_in_log, y_in_log_predicted, sample_weight =weight)
    r2_in_log = r2_score(y_in_log, y_in_log_predicted, sample_weight =weight)
    
    regfp = open(str(output)+".regression.log", "a+")
    sumYminusYbar = 0
    
    for i in range (0, num_bins):
        sumYminusYbar  = (np.square(y_in_log[i] - y_in_log_predicted[i]))
    SEE_log = sqrt(sumYminusYbar / (num_bins-2) )
    print(str(regression_model.coef_[0][0]),str(regression_model.intercept_[0]),str(rmse_in_log),str(r2_in_log), str(SEE_log),file=regfp)
    # printing values
    print("----------------------------------------")
    print("Dropoff Rate:", regression_model.coef_)
    print("Standard Error is ", SEE_log)
    print("Root mean squared error: ", rmse_in_log)
    print("R2 score: ", r2_in_log)
    print("----------------------------------------")
    if plot ==1 : 

        plt.figure(figsize=(10, 10))
        plt.ylim(ylogmin, ylogmax)
        plt.scatter(x, y_in_log, s=norm_weight)
        xtext = (" Slope: "+str(regression_model.coef_)+ " RMSE: "+ str(rmse_in_log)+ " R2 score: "+ str(r2_in_log) +" SEE: " +str(SEE_log) )
        plt.xlabel(str(label) + "\n" + str(xtext), fontsize=12)
        plt.ylabel("Bin Value", fontsize=12)
        # predicted values
        plt.plot(x, y_in_log_predicted, color="r")
        plt.title(output + " Weighted Linear Regression", fontsize=12)
        plt.savefig(output + ".Log.WLR.png", format="png")
        plt.clf()

    #WLR in linear mode 
    regression_model.fit(x, y_in_linear, sample_weight = weight)
    # Predict
    y_in_linear_predicted = regression_model.predict(x)

    # model evaluation
    rmse_in_linear = mean_squared_error(y_in_linear, y_in_linear_predicted, sample_weight =weight)
    r2_in_linear = r2_score(y_in_linear, y_in_linear_predicted, sample_weight =weight)
    sumYminusYbar = 0
    for i in range (0, num_bins):
        sumYminusYbar  = (np.square(y_in_linear[i] - y_in_linear_predicted[i]))
    SEE_linear = sqrt(sumYminusYbar / (num_bins-2) )
    print('Regression in linear mode') 
    print("----------------------------------------")
    print("Dropoff Rate:", regression_model.coef_)
    print("Standard Error is ", SEE_linear)
    print("Root mean squared error: ", rmse_in_linear)
    print("R2 score: ", r2_in_linear)
    print("----------------------------------------")
    if plot ==1 :

        plt.figure(figsize=(10, 10))
        plt.ylim(ylogmin, ylogmax)
        plt.scatter(x, y_in_linear, s=norm_weight)
        xtext = (" Slope: "+str(regression_model.coef_)+ " RMSE: "+ str(rmse_in_linear)+ " R2 score: "+ str(r2_in_linear)+ " SEE: "+str(SEE_linear) )
        plt.xlabel(str(label) + "\n" + str(xtext), fontsize=12)
        plt.ylabel("Bin Value", fontsize=12)
        # predicted values
        plt.plot(x, y_in_linear_predicted, color="r")
        plt.title(output + " Weighted Linear Regression", fontsize=12)
        plt.savefig(output + ".Linear.WLR.png", format="png")
        plt.clf()
    
# ---------------------------------------------------------------------------
# Main function: gets input paramters and calls corresponding functions
# ---------------------------------------------------------------------------


def main():
    print("Initializing and reading arguments")
    parser = argparse.ArgumentParser()
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-t','--transcripts', required=True, help="Transcripts file")
    optional.add_argument('--optional_arg')
    required.add_argument('-f','--footprint', default=False, help="Ribosomes bed file")
    optional.add_argument("-r", "--rnaseq", dest="rnaseq", default="NULL",help="mRNA bed file, if mRNA bed is not provided, ribosomes dropoff won't be normalized with mRNA")
    optional.add_argument("-s", "--subset", dest="subset", default="NULL",help="subset of genes to run the analysis on")
    optional.add_argument("-b", "--binsize", dest="binsize", type=int, default=50, help="Bin size default is 50")
    optional.add_argument("-o", "--out", dest="output", default="", help="Output file name")
    optional.add_argument("-p", "--plot", dest="plot",type=int, default=1, help="Plotting mode is on by default, use --plot 0 to turn off plots")
    optional.add_argument("--ymin", dest="ymin", type=int, default=-3,help="ymin for the y axis min position in linear plot")
    optional.add_argument("--ymax", dest="ymax", type=int, default=2, help="ymax for the y axis max position in linear plot")
    optional.add_argument("--ylogmin", type=int, default=-3, help="ylogmin for y axis min position in log plot")
    optional.add_argument("--ylogmax", type=int, default=2, help="ylogmax for y axis max position in log plot")
    args = parser.parse_args()
    plot = args.plot 
    if plot ==1:
        print("Plot mode is on, regression plots will be printed")
    else: 
        print("Plot mode is off, no plots will be printed") 
    '''
     # Check required files exist
    if not os.path.exists(args.transcripts):
        sys.exit("transcripts file is required and does not exist, exiting !")
    if not os.path.exists(args.footprint):
        sys.exit("footprint file is required and does not exist, exiting !")
    if args.rnaseq != "NULL":
        if not os.path.exists(args.rnaseq):
            sys.exit(
                "rna  file does not exist: either run with footprint only option or enter a valid rna file, exiting !"
            )
   ''' 
    sample = str(args.footprint).split(".")[0]
    if args.rnaseq != "NULL":
        sample += "_" + str(args.rnaseq).split(".")[0]
    
    plot = args.plot

    output = args.output
    if output == "":
        output = sample
        if args.subset != "NULL":
            subset_name = args.subset.split(".")[0]
            output += "." + subset_name

    binsize = int(args.binsize)
    
    print("Reading transcripts")
    if args.subset == "NULL":
        max_gene_length, genes_length = get_genes(args.transcripts)
    else:
        max_gene_length, genes_length = get_subset_genes(args.transcripts, args.subset)
    print("Done reading transcripts")  
    gene_coverage_at_bin = get_gene_coverage_at_bin(max_gene_length, binsize, genes_length)
    print("Reading ribosome footprint", args.footprint)
    fp_coverage = get_reads(args.footprint, genes_length) 
    print("Reading mRNA",args.rnaseq) 
    mRNA_coverage  = get_reads(args.rnaseq, genes_length) 
    print("Filling coverage matrix") 
    fp_gene_coverage_at_pos = get_gene_coverage_at_pos(max_gene_length, fp_coverage, genes_length) 
    mRNA_gene_coverage_at_pos = get_gene_coverage_at_pos(max_gene_length, mRNA_coverage, genes_length)
   
    print("Filling positions matrix")

    fp_positions = fill_positions(fp_coverage, max_gene_length) 
    mRNA_positions = fill_positions(mRNA_coverage, max_gene_length) 

    print("Started binning process") 
    ribosomes_gene_bins = binning(
        binsize,
        fp_positions,
        fp_gene_coverage_at_pos,
        max_gene_length
    )
    
    mRNA_gene_bins = binning(
        binsize,
        mRNA_positions,
        mRNA_gene_coverage_at_pos,
        max_gene_length 
    )
    

    print("No. of bins:", len(ribosomes_gene_bins))

    num_bins = int(max_gene_length / binsize) + 1
    gene_bins = ribosomes_gene_bins
    if args.rnaseq != "NULL":
        for i in range(0, num_bins):
            gene_bins[i] = float(ribosomes_gene_bins[i] / mRNA_gene_bins[i])
    print("Done binning, starting regression and plotting")
    # -----------------------------------------
    # Plotting and summarizing
    # -----------------------------------------
    regression(
        output,
        num_bins,
        fp_coverage,
        gene_bins,
        binsize,
        max_gene_length,
        genes_length,
        int(args.ymin),
        int(args.ymax),
        int(args.ylogmin),
        int(args.ylogmax),
        gene_coverage_at_bin,plot 
    )
    
    print("All is done")
    

if __name__ == "__main__":
    main()
