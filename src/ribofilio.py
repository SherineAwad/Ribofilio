#! /usr/bin/env python
import sys
import argparse
import screed
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use("Agg")
from pylab import *
import os.path
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score

# -----------------------------------------------
# Read and parse input files
# -----------------------------------------------
def get_subset_genes(transcripts_file, subset_file):
    subset = []
    max_gene_length = -100 
    genes_length = {} 
    for line in open(subset_file):
        gene_name = line.rstrip()
        if "mRNA" in gene_name:
            gene_name=gene_name.split('_')[0]
        subset.append(gene_name) 
    for record in screed.open(transcripts_file):
        gene_name = record.name.split(' ')[0]
        if "mRNA" in gene_name:
            gene_name=gene_name.split('_')[0]
        if gene_name in subset:
            genes_length[str(gene_name)] = int(len(record.sequence))
            if genes_length[gene_name] > max_gene_length:
                max_gene_length = genes_length[gene_name]
    print ('max_gene_length is',max_gene_length)
    return max_gene_length, genes_length




def get_genes(transcripts_file): 
    max_gene_length = -100 
    genes_length = {} 
    for record in screed.open(transcripts_file):
        gene_name = record.name.split(' ')[0]
        if "mRNA" in gene_name:
            gname=gname.split('_')[0]
        genes_length[str(gene_name)] = int(len(record.sequence))
        if genes_length[gene_name] > max_gene_length:
            max_gene_length = genes_length[gene_name]
    print ('max_gene_length is',max_gene_length)
    return max_gene_length, genes_length 


def get_reads(infile, genes_length):
    coverage = {} 
    for gene in genes_length: 
        coverage[gene] = [] 
    for  line in open(infile):
        record = line.rstrip().split("\t")
        gene_name = str(record[0])
        if "mRNA" in gene_name:
            gene_name = gene_name.split("_")[0]
        if gene_name in genes_length: 
            coverage[str(gene_name)].append(int(record[2]) )
    return coverage 

# ------------------------------------------------------------------------
# ribosomes_profile function: estimates the drop rate of ribosomes after binning
# ------------------------------------------------------------------------
def ribosomes_profile(
     bin_size, coverage, genes_length, max_gene_length, min_gene_length
):
    positions = [0] * 100000000  
    gene_bins = []
    c = 0.000001
    gene_coverage_at_pos = [0] * 100000000
    j = 0
    
    num_bins = int(max_gene_length / bin_size) + 1
    gene_bins = [0] * num_bins
    normalized_positions = [0] * max_gene_length
    gene_coverage_at_bin = [0] * num_bins

    # Fill gene_coverage_at_bin
    for gene in genes_length:
        bin_fit = math.ceil(genes_length[gene] / bin_size)
        for i in range(0, bin_fit):
            gene_coverage_at_bin[i] += 1

    # Counts genes in each position
    for gene in coverage:
        for i in coverage[gene]:
            if (
                genes_length[gene] <= max_gene_length
                and genes_length[gene] > min_gene_length
            ):
                positions[int(i)] += 1
    print("Filling pos is done")

    for gene in coverage:
        if (
            genes_length[gene] <= max_gene_length
            and genes_length[gene] > min_gene_length
        ):
            for i in range(0, genes_length[gene] + 1):
                gene_coverage_at_pos[i] += 1
    print("Filling gCoverered is done")
    i = 0
    last_pos = 0
    # Normalize bin position  with the number of gene covering that position
    while i < max_gene_length:
        normalized_positions[int(i)] = positions[int(i)] / (gene_coverage_at_pos[int(i)] + c)
        last_pos = i
        i += 1
    print("last_pos is ", last_pos)
    print("Normalizing pos is done")
    
    # Here, we group reads into bins and normalize by bin_size, we stop at last_pos
    index = 0
    a = 0
    b = bin_size
    while a < max_gene_length:
        for i in range(a, b + 1):
            if i > (len(normalized_positions) - 1):
                break
            if i > last_pos:
                break
            position_sum = float(normalized_positions[i])
            gene_bins[index] += position_sum
        gene_bins[index] = float(c + (gene_bins[index] / bin_size))
        index += 1
        a = b + 1
        b = b + bin_size
    return gene_bins, last_pos, gene_coverage_at_bin
# -------------------------------------------------------------------------------------------
# plots function, plot several figures for the ribosome profiling
# -------------------------------------------------------------------------------------------
def plots(
    plot_flag,
    output,
    num_bins,
    coverage,
    all_bins,
    bin_size,
    max_gene_length,
    min_gene_length,
    last_pos,
    genes_length,
    y_min,
    y_max,
    ylogmin,
    ylogmax,
    gene_coverage_at_bin, 
    regfile,
):

    print("Start plotting")
    bins = []
    log_gene_bins = []
    gene_bins = []
    output = output
    genesfile = output + ".genes"
    print("num_bins is: ", num_bins)
    # Take log of bins indeces for LOG/LOG plotting
    for i in range(0, num_bins):
        bins.append(i)
        gene_bins.append(all_bins[i])
    log_gene_bins.append(0)
    i = 0
    for i in range(1, num_bins):
        log_gene_bins.append(np.log(i))
    print("Prining bins, gene_bins and logene_bins:")
    print(
        "lengths of bins, logene_bins, and gene_bins is:",
        len(bins),
        len(log_gene_bins),
        len(gene_bins),
    )
    label = "Bins  (Binsize = " + str(bin_size) + " Max gene Length =" + str(last_pos)

    if plot_flag == 1:
        # Plot Coverage without Log
        # ----------------------------
        plt.figure(figsize=(15, 10))
        plt.plot(np.array(bins), np.array(gene_bins), color="Navy")
        plt.ylim(y_min, y_max)
        plt.ylabel("Sums of coverage of genes at each position in each bin")
        plt.xlabel(label)
        plt.title(output + " Coverage per Bin (Sum without log) ", fontsize=12)
        plt.savefig(output + "_" + str(bin_size) + ".NoLog.png", format="png")
        plt.clf()

        # Plot Coverage Log/Log
        # --------------------------
        plt.figure(figsize=(15, 10))
        plt.plot(np.array(log_gene_bins), np.array(np.log(gene_bins)), color="Navy")
        plt.ylabel("Log of sums of coverage of genes at each position in each bin")
        plt.ylim(ylogmin, ylogmax)
        plt.xlabel(label)
        plt.title(output + " Coverage per Bin (Log/Log) ", fontsize=12)
        plt.savefig(output + "_" + str(bin_size) + ".LogLog.png", format="png")
        clf()
        # Plot coverage Log /Linear
        # -----------------------------
        plt.figure(figsize=(15, 10))
        plt.plot(np.array(bins), np.log(gene_bins), color="Navy")
        plt.ylim(ylogmin, ylogmax)
        plt.ylabel("Log of sums of coverage of genes at each position in each bin")
        plt.xlabel(label)
        plt.title(output + " Coverage per Bin (Log/Linear) ", fontsize=12)
        plt.savefig(output + "_" + str(bin_size) + ".LogLinear.png", format="png")
        plt.clf()
    # Linear Regression
    # ------------------
    x = np.array(bins).reshape(-1, 1)
    y = np.array(np.log(gene_bins)).reshape(-1, 1)
    regression_model = LinearRegression()
    W = bin_size
    weight = bins
    for i in range(0, len(weight)):
        weight[i] = gene_coverage_at_bin[i]
    # Fit the data(train the model)
    regression_model.fit(x, y, weight)
    # Predict
    y_predicted = regression_model.predict(x)

    # model evaluation

    rmse = mean_squared_error(y, y_predicted, weight)
    r2 = r2_score(y, y_predicted, weight)

    # printing values
    print("Slope:", regression_model.coef_)
    print("Intercept:", regression_model.intercept_)
    print("Root mean squared error: ", rmse)
    print("R2 score: ", r2)
    plt.figure(figsize=(15, 10))
    plt.ylim(ylogmin, ylogmax)
    plt.scatter(x, y, s=10)
    xtext = (
        " Slope: "
        + str(regression_model.coef_)
        + " Intercept: "
        + str(regression_model.intercept_)
        + " RMSE: "
        + str(rmse)
        + " R2 score :"
        + str(r2)
    )
    plt.xlabel(str(label) + "\n" + str(xtext), fontsize=10)
    plt.ylabel("y")
    # predicted values
    plt.plot(x, y_predicted, color="r")
    plt.title(output + " Linear Regression", fontsize=12)
    plt.savefig(output + "_" + str(bin_size) + ".WLR.png", format="png")
    regfp = open(str(regfile), "a+")
    print(
        str(regression_model.coef_[0][0]),
        str(regression_model.intercept_[0]),
        str(rmse),
        str(r2),
        file=regfp,
    )
    plt.clf()

# ---------------------------------------------------------------------------
# Main function: gets input paramters and calls corresponding functions
# ---------------------------------------------------------------------------


def main():
    print("Initializing and reading arguments")
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-t", "--transcripts_file", dest="transcripts_file", default=False
    )
    parser.add_argument("-f", "--footprint", dest="footprint", default=False)
    parser.add_argument("-r", "--rnaseq", dest="rnaseq", default="NULL")
    parser.add_argument("-s", "--subset_file", dest="subset_file", default="NULL")
    parser.add_argument("-b", "--bin_size", dest="bin_size", type=int, default=50)
    parser.add_argument("--plots", dest="plots", type=int, default=1)
    parser.add_argument("-o", "--out", dest="output", default="")
    parser.add_argument("-n", "--name", dest="reg_file_name", default="regfile.log")
    parser.add_argument("--y_min", dest="y_min", type=int, default=0)
    parser.add_argument("--y_max", dest="y_max", type=int, default=30)
    parser.add_argument("--ylogmin", type=int, default=-15)
    parser.add_argument("--ylogmax", type=int, default=5)
    args = parser.parse_args()
    
    '''
     # Check required files exist
    if not os.path.exists(args.transcripts_file):
        sys.exit("transcripts file is required and does not exist, exiting !")
    if not os.path.exists(args.footprint):
        sys.exit("footprint file is required and does not exist, exiting !")
    if args.rnaseq != "NULL":
        if not os.path.exists(args.rnaseq):
            sys.exit(
                "rna  file does not exist: either run with footprint only option or enter a valid rna file, exiting !"
            )
   ''' 
    plot_flag = args.plots
    print("plot_flag is", plot_flag)
    
    sample = str(args.footprint).split(".")[0]
    if args.rnaseq != "NULL":
        sample += "_" + str(args.rnaseq).split(".")[0]
    
    output = args.output
    if output == "":
        output = sample
        if args.subset_file != "NULL":
            subset_name = args.subset_file.split(".")[0]
            output += "." + subset_name

    bin_size = int(args.bin_size)

    if args.subset_file == "NULL":
        max_gene_length, genes_length = get_genes(args.transcripts_file)
    else:
        max_gene_length, genes_length = get_subset_genes(args.transcripts_file, args.subset_file)
   
    
    print(args.footprint)
    print(args.rnaseq) 


    fp_coverage = get_reads(args.footprint, genes_length) 
    mRNA_coverage  = get_reads(args.rnaseq, genes_length) 


    min_gene_length = 0
    
    print("We are calling ribosomes_profile for", args.footprint)
    ribosomes_gene_bins, last_pos, gene_coverage_at_bin = ribosomes_profile(
        bin_size,
        fp_coverage,
        genes_length,
        max_gene_length,
        min_gene_length
    )
    
    mRNA_gene_bins, _, _ = ribosomes_profile(
        bin_size,
        mRNA_coverage, 
        genes_length, 
        max_gene_length, 
        min_gene_length
    )


    print("Length of ribosomes_gene_bins is ", len(ribosomes_gene_bins))
    print("Length of mRNA_gene_bins is ", len(mRNA_gene_bins))
    

    num_bins = int(last_pos / bin_size) + 1
    gene_bins = ribosomes_gene_bins
    if args.rnaseq != "NULL":
        for i in range(0, num_bins):
            gene_bins[i] = float(ribosomes_gene_bins[i] / mRNA_gene_bins[i])
    print('gene_bins is', gene_bins)
    # -----------------------------------------
    # Plotting and summarizing
    # -----------------------------------------
    plots(
        plot_flag,
        output,
        num_bins,
        fp_coverage,
        gene_bins,
        bin_size,
        max_gene_length,
        min_gene_length,
        last_pos,
        genes_length,
        int(args.y_min),
        int(args.y_max),
        int(args.ylogmin),
        int(args.ylogmax),
        gene_coverage_at_bin, 
        args.reg_file_name,
    )
    
    print("Done Plotting and summarizing")
    print("All is done")
    

if __name__ == "__main__":
    main()
