#! /usr/bin/env python
import argparse
import os
import sys
import math
import screed
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import t

# Read subsets of genes --subset option


def get_subset_genes(transcripts, subset_file):
    print("Processing ", subset_file)
    subset = []
    found = 0
    max_gene_length = -100
    genes_length = {}
    for line in open(subset_file):
        gene_name = line.rstrip()
        if "mRNA" in gene_name:
            gene_name = gene_name.split('_')[0]
        subset.append(gene_name)
    for record in screed.open(transcripts):
        gene_name = record.name.split(' ')[0]
        if gene_name in subset:
            found += 1
            genes_length[str(gene_name)] = int(len(record.sequence))
            if genes_length[gene_name] > max_gene_length:
                max_gene_length = genes_length[gene_name]
    try:
        assert found > 0
    except AssertionError:
        print("All genes in", subset_file,
              "doesn't exist in transcripts file.")
        sys.exit()
    return max_gene_length, genes_length

# Read Transcripts


def get_genes(transcripts):
    max_gene_length = -100
    genes_length = {}
    itr = 0
    for record in screed.open(transcripts):
        gene_name = record.name.split(' ')[0]
        itr += 1
        if (itr % 1000) == 0:
            print(".....................")
        if "mRNA" in gene_name:
            gene_name = gene_name.split('_')[0]
        genes_length[str(gene_name)] = int(len(record.sequence))
        if genes_length[gene_name] > max_gene_length:
            max_gene_length = genes_length[gene_name]
    return max_gene_length, genes_length

# Read bed files to get reads and positions of alignment


def get_reads(infile, genes_length):
    print("Processing ", infile)
    coverage = {}
    itr = 0
    found = 0
    for gene in genes_length:
        coverage[gene] = []
    for line in open(infile):
        record = line.rstrip().split("\t")
        itr += 1
        if (itr % 5000000) == 0:
            print(".....................")
        gene_name = str(record[0])
        if "mRNA" in gene_name:
            gene_name = gene_name.split("_")[0]
        if gene_name in genes_length:
            coverage[str(gene_name)].append(int(record[2]))
            found += 1
    try:
        assert found > 0
    except AssertionError:
        print("All genes in", infile, "doesn't exist in transcripts file.")
        sys.exit()
    return coverage

# Fill bins with genes coverage in each bin


def get_gene_coverage_at_bin(max_gene_length, binsize, genes_length):
    num_bins = int(max_gene_length / binsize) + 1
    gene_coverage_at_bin = [0] * num_bins
    # Fill gene_coverage_at_bin
    for gene in genes_length:
        bin_fit = math.ceil(genes_length[gene] / binsize)
        for i in range(0, bin_fit):
            gene_coverage_at_bin[i] += 1
    return gene_coverage_at_bin

# How many positions a gene can cover based on its length


def get_gene_coverage_at_pos(max_gene_length, coverage, genes_length):
    gene_coverage_at_pos = [0] * (max_gene_length + 1)
    for gene in coverage:
        for i in range(1, genes_length[gene] + 1):
            gene_coverage_at_pos[i] += 1
    return gene_coverage_at_pos

# Fill position matrix with genese' reads covering each position


def fill_positions(coverage, max_gene_length):
    positions = [0] * (max_gene_length + 1)
    for gene in coverage:
        for i in coverage[gene]:
            positions[int(i)] += 1
    return positions


# ------------------------------------------------------------------------
# binning function: estimates the drop rate of ribosomes after binning
# ------------------------------------------------------------------------


def binning(binsize, positions, gene_coverage_at_pos, max_gene_length):
    gene_bins = []
    const_c = 0.000001
    num_bins = int(max_gene_length / binsize) + 1
    gene_bins = [0] * num_bins
    normalized_positions = [0] * (max_gene_length + 1)
    i = 1
    # Normalize bin position  with the number of gene covering that position
    while i <= max_gene_length:
        normalized_positions[int(i)] = (
             positions[int(i)] / (gene_coverage_at_pos[int(i)] + const_c))
        i += 1
    # Group reads into bins and normalize by binsize,
    # We stop at last_pos
    index = 0
    window_a = 1
    window_b = binsize
    while window_a <= max_gene_length:
        for i in range(window_a, window_b + 1):
            if i > (len(normalized_positions) - 1):
                break
            position_sum = float(normalized_positions[i])
            gene_bins[index] += position_sum
        gene_bins[index] = const_c + float(gene_bins[index] / binsize)
        index += 1
        window_a = window_b + 1
        window_b = window_b + binsize
    return gene_bins


# -------------------------------------------------------------------------------------------
# Regression function, plot several figures for the ribosome profiling
# -------------------------------------------------------------------------------------------


def regression(output, num_bins, gene_bins,
               binsize, ylogmin, ylogmax,
               gene_coverage_at_bin, plot):
    bins = []
    # Take log of bins indeces for LOG/LOG plotting
    for i in range(0, num_bins):
        bins.append(i)
    codons_bin = 0.0   # codons_bin is number of codons per bin
    codons_bin = float(binsize / 3)
    # Weighted Linear Regression
    x_value = np.array(bins).reshape(-1, 1)
    y_value = np.array(np.log(gene_bins)).reshape(-1, 1)
    regression_model = LinearRegression()
    weight = [0] * num_bins
    norm_weight = [0] * num_bins
    for i in range(0, len(weight)):
        weight[i] = gene_coverage_at_bin[i]
        norm_weight[i] = np.cbrt(weight[i])
    #  Fit the data(train the model)
    regression_model.fit(x_value, y_value, sample_weight=weight)
    #  Predict
    y_predicted = regression_model.predict(x_value)
    fp = open(str(output)+".bins.csv", "a+")
    print("X", '\t', "Y", file=fp)
    for i in range(0, len(y_predicted)):
        print(x_value[i][0], '\t', y_predicted[i][0], file=fp)
    #  Model evaluation
    rmse = mean_squared_error(y_value, y_predicted, sample_weight=weight)
    rsquare = r2_score(y_value, y_predicted, sample_weight=weight)
    # Calculate dropoff_codon is dropoff rate per codon
    dropoff_codon = 1 - pow((1 - regression_model.coef_), (1/codons_bin))
    dropoff_rate = regression_model.coef_[0][0]
    # Calculate Standard error, margin error, and a confidence interval of 95%
    ebsilon = 0
    df = 2
    alpha = 1 - (95 / 100)
    critical_p = 1 - (alpha/2)
    critical_p = t.ppf(critical_p, (num_bins-df))
    mean = np.mean(np.array(bins))
    for i in bins:
        ebsilon += np.square(i-mean)
    xmean = np.sqrt(ebsilon)
    sumy_ypredicted = 0
    for i in range(0, num_bins):
        sumy_ypredicted = (sumy_ypredicted
                           + np.square(y_value[i] - y_predicted[i]))
    stand_error = np.sqrt((sumy_ypredicted / (num_bins - df))) / xmean
    margin_error = (critical_p * stand_error)
    # Calculate tscore and pvalue of
    # how different the slope is from a slope of zero
    tscore = regression_model.coef_[0][0]/stand_error
    pvalue = t.sf(abs(tscore), df=(num_bins-df))
    # Do some rounding and print to both file and screen
    stand_error = np.round(stand_error, decimals=4)
    margin_error = np.round(margin_error, decimals=4)
    rsquare = np.round(rsquare, decimals=4)
    rmse = np.round(rmse, decimals=4)
    dropoff_rate = np.round(regression_model.coef_[0][0], decimals=4)
    dropoff_codon = np.round(dropoff_codon, decimals=4)
    pvalue = np.round(pvalue, decimals=4)
    regfp = open(str(output)+".regression.log", "a+")
    print("Dropoff\tDropoff per codon \tRMSE\tRsquare\tSE" +
          "\tMargin Error\ttscore\tpvalue\tNo.of Bins", file=regfp)
    print(str(dropoff_rate) + "\t" + str(dropoff_codon[0][0]) +
                              "\t" + str(rmse) +
                              "\t" + str(rsquare) +
                              "\t" + str(stand_error[0]) +
                              "\t" + str(margin_error[0]) +
                              "\t" + str(tscore[0]) +
                              "\t" + str(pvalue[0]) +
                              "\t" + str(num_bins), file=regfp)
    print("----------------------------------------")
    print("Weighted Linear Regression")
    print("----------------------------------------")
    print("Dropoff Rate :", dropoff_rate)
    print("Dropoff per codon:", dropoff_codon)
    print("Standard Error is (SE):", stand_error)
    print("Confidence Interval is: [", dropoff_rate, "+/-", margin_error, "]")
    print("Root Mean Squared Error (RMSE): ", rmse)
    print("R2 Score: ", rsquare)
    print("T-test score: ", tscore)
    print("Pvalue:", pvalue)
    regfp.close()
    if plot == 1:
        plot_regression(x_value, y_value, y_predicted, norm_weight,
                        str(dropoff_rate), str(dropoff_codon[0][0]),
                        str(rmse), str(rsquare),
                        str(stand_error[0]), output,
                        ylogmin, ylogmax)
    return (dropoff_rate, dropoff_codon, stand_error,
            margin_error, rmse, rsquare, tscore, pvalue)
# ------------------------------------------
# Plot Regression
# -----------------------------------------


def plot_regression(x_value, y_value, y_predicted,
                    norm_weight, dropoff_rate, dropoff_codon, rmse, rsquare,
                    stand_error, output, ymin, ymax):
    label = "Bin Number"
    fig = plt.figure(figsize=(10, 10))
    plt.ylim(ymin, ymax)
    plt.scatter(x_value, y_value, s=norm_weight)
    xtext = (" Dropoff: " + str(dropoff_rate) +
             " Dropoff per codon: " + str(dropoff_codon) +
             "\n RMSE: " + str(rmse) +
             " RSquare: " + str(rsquare) + " SE: " + str(stand_error))
    plt.xlabel(str(label) + "\n" + str(xtext), fontsize=10)
    plt.ylabel("Bin Value", fontsize=10)
    plt.plot(x_value, y_predicted, color="r")
    plt.title(output + " Weighted Linear Regression", fontsize=10)
    plt.savefig(output + ".Log.WLR.png", format="png")
    plt.clf()
    return fig

# --------------------------------------------
# Call mRNA
# --------------------------------------------


def call_mRNA(rnaseq, genes_length, max_gene_length, binsize):
    rna_gene_bins = []
    if rnaseq != "NULL":
        print("Reading and binning mRNA ...")
        rna_coverage = get_reads(rnaseq, genes_length)
        rna_gene_coverage_at_pos = (get_gene_coverage_at_pos(
                                    max_gene_length, rna_coverage,
                                    genes_length))
        rna_positions = fill_positions(rna_coverage, max_gene_length)
        rna_gene_bins = binning(binsize, rna_positions,
                                rna_gene_coverage_at_pos, max_gene_length)
    return rna_gene_bins

# -------------------------------------------
# Call footprint
# -------------------------------------------


def call_footprints(footprint, genes_length, max_gene_length, binsize):
    fp_coverage = get_reads(footprint, genes_length)
    fp_gene_coverage_at_pos = (get_gene_coverage_at_pos
                               (max_gene_length, fp_coverage, genes_length))
    fp_positions = fill_positions(fp_coverage, max_gene_length)
    ribosomes_gene_bins = binning(binsize, fp_positions,
                                  fp_gene_coverage_at_pos,
                                  max_gene_length)
    return ribosomes_gene_bins


# -----------------------------------------------
# Normalize Footprint using mRNA
# -----------------------------------------------


def normalize(ribosomes_gene_bins, rna_gene_bins, num_bins):
    gene_bins = ribosomes_gene_bins
    for i in range(0, num_bins):
        gene_bins[i] = float(ribosomes_gene_bins[i] / rna_gene_bins[i])
    return gene_bins

# ---------------------------------------------
# Parse arguments
# ---------------------------------------------


def get_arguments():
    print("Initializing and reading arguments")
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--transcripts',
                        required=True, help="Transcripts file")
    parser.add_argument('-f', '--footprint',
                        required=True, help="Ribosomes bed file")
    parser.add_argument("-r", "--rnaseq", dest="rnaseq", default="NULL",
                        help="mRNA bed file, if mRNA bed is not provided," +
                        "ribosomes dropoff won't be normalized with mRNA")
    parser.add_argument("-s", "--subset", dest="subset", default="NULL",
                        help="subset of genes to run the analysis on")
    parser.add_argument("-b", "--binsize", dest="binsize",
                        type=int, default=50,
                        help="Bin size default is 50")
    parser.add_argument("-o", "--out", dest="output", default="",
                        help="Output file name")
    parser.add_argument("-p", "--plot", dest="plot", type=int, default=1,
                        help="Plotting mode is on by default," +
                        "use --plot 0 to turn off plots")
    parser.add_argument("--ylogmin", type=int, default=-3,
                        help="ylogmin for y axis min position in log plot")
    parser.add_argument("--ylogmax", type=int, default=2,
                        help="ylogmax for y axis max position in log plot")
    return parser
# ---------------------------------------------------------------------------
# Main function: gets input paramters and calls corresponding functions
# ---------------------------------------------------------------------------


def main():
    parser = get_arguments()
    args = parser.parse_args()
    # Check required files exist
    try:
        assert os.path.exists(args.transcripts)
        assert os.path.getsize(args.transcripts) > 0
    except AssertionError:
        sys.exit("transcripts file does not exist or empty, exiting !")
    try:
        assert os.path.exists(args.footprint)
        assert os.path.getsize(args.footprint) > 0
    except AssertionError:
        sys.exit("footprint file does not exist or empty, exiting !")
    if args.rnaseq != "NULL":
        try:
            assert os.path.exists(args.rnaseq)
            assert os.path.getsize(args.rnaseq) > 0
        except AssertionError:
            sys.exit("mRNA  file does not exist or empty. " +
                     "Either run with footprint only option" +
                     " or enter a valid rna file, exiting !")
    if args.subset != "NULL":
        try:
            assert os.path.exists(args.subset)
            assert os.path.getsize(args.subset) > 0
        except AssertionError:
            sys.exit("Subset file is empty or doesn't exist")
    sample = str(args.footprint).split(".")[0]
    if args.rnaseq != "NULL":
        sample += "_" + str(args.rnaseq).split(".")[0]
    plot = args.plot
    output = args.output
    # Assign a default output name if none is input
    if output == "":
        output = sample
        if args.subset != "NULL":
            subset_name = args.subset.split(".")[0]
            output += "." + subset_name
    binsize = int(args.binsize)
    if args.subset == "NULL":
        max_gene_length, genes_length = get_genes(args.transcripts)
    else:
        max_gene_length, genes_length = (get_subset_genes
                                         (args.transcripts, args.subset))
    gene_coverage_at_bin = (get_gene_coverage_at_bin(max_gene_length,
                            binsize, genes_length))
    # Call footprint will do the counts and binning
    print("Reading and binning footprints ...")
    ribosomes_gene_bins = call_footprints(args.footprint,
                                          genes_length,
                                          max_gene_length, binsize)
    # Call mRNA if exist will do the counts and binning
    if args.rnaseq != "NULL":
        rna_gene_bins = call_mRNA(args.rnaseq,
                                  genes_length,
                                  max_gene_length, binsize)
    gene_bins = ribosomes_gene_bins
    # Normalize footprint with mRNA
    num_bins = int(max_gene_length / binsize) + 1
    if args.rnaseq != "NULL":
        gene_bins = normalize(
                              ribosomes_gene_bins,
                              rna_gene_bins, num_bins)
    # Plotting and summarizing
    regression(
        output,
        num_bins,
        gene_bins,
        binsize,
        int(args.ylogmin),
        int(args.ylogmax),
        gene_coverage_at_bin, plot
    )
    print("All is done")


if __name__ == "__main__":
    main()
