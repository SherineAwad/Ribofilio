import pytest
import screed
import os
import sys 
import numpy as np 

path = os.getcwd() 
path  = os.path.join(path,"src") 
sys.path.append(path)

print(path)
import ribofilio as rb  


def test_get_genes():
    path = os.getcwd ()
    path =  os.path.join(path, "tests/test-data","transcripts_sample.fa")
    max_gene_length, gene_length = rb.get_genes(path) 
    assert (max_gene_length == 60 ) 
    assert (gene_length ==  {"YBR024W":16, "YDL245C":60,"YBR021W":8}) 

def test_get_subsets(): 
    path = os.getcwd ()
    file1 =  os.path.join(path, "tests/test-data","transcripts_sample.fa")
    file2 = os.path.join(path, "tests/test-data","subset.txt")
    max_gene_length, gene_length = rb.get_subset_genes(file1,file2) 
    
    assert (max_gene_length ==16) 
    assert (gene_length ==  {"YBR024W":16,"YBR021W":8})

def test_get_reads():
    path = os.getcwd()
    path = os.path.join(path, "tests/test-data", "sample.bed") 
    coverage = rb.get_reads(path,  {"YKL152C":50}) 
    assert (coverage == {"YKL152C":[93]})

def test_get_gene_coverage_at_bin():
    path = os.getcwd()
    genes_length = {"YBR024W":52,"YBR021W":45}
    max_gene_length = 52 
    bin_size =50
    assert(rb.get_gene_coverage_at_bin(max_gene_length, bin_size, genes_length) == [2,1])

def test_get_gene_coverage_at_pos():
    max_gene_length = 6
    coverage = {"YBR024W":[1,2,3], "YBR021W":[5,6,7]}
    gene_coverage_at_pos = [0] * (max_gene_length +1)
    genes_length = {"YBR024W":6,"YBR021W":4, "YKL152C":23}
    assert (rb.get_gene_coverage_at_pos(max_gene_length, coverage, genes_length) == [0,2,2,2,2,1,1]) 
   

def test_fill_positions ():
    max_gene_length = 9
    positions = [0] * (max_gene_length + 1 )
    coverage = {"YKL152C":[5,7,9], "YBR021W":[3,6,7], "YBR024W":[1,2,3]}
    assert(rb.fill_positions(coverage, max_gene_length) == [0,1,1,2,0,1,1,2,0,1] )


def test_binning():
     positions =[0,2/3, 1,2/3] 
     gene_coverage_at_pos = [0,2,3,2]
     max_gene_length = 3 
     genes_bin  = rb.binning(2,positions, gene_coverage_at_pos, max_gene_length)
     round_genes_bin  = np.round(genes_bin, 6)
     print(round_genes_bin)
     assert (round_genes_bin  == [0.333334, 0.166668]).all()


def test_regression(): 
    output ="test"
    num_bins = 18
    all_bins = [6612, 6612, 6569, 6481, 6384, 6251, 6130, 5757, 5446, 5226, 5039, 4888, 4716, 4550, 4393, 4233, 4105, 3969]
    binsize = 2
    plot = 0 
    gene_coverage_at_bin = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] 
    ylogmin =-2  
    ylogmax=2
    dropoff_rate, dropoff_codon, stand_error, margin_error, rmse, rsquare, tscore, pvalue = rb.regression(output, num_bins, all_bins,
               binsize, ylogmin, ylogmax,
               gene_coverage_at_bin, plot)
    assert(dropoff_rate == -0.0334)
    assert(dropoff_codon ==[[-0.0506]]) 
    assert(stand_error == [0.0087]) 
    assert(margin_error == [0.0185])
    assert(rmse == 0.033) 
    assert(rsquare == -0.0682)
    assert(np.round(tscore[0],4) == -3.8232)
    assert(pvalue == [0.0007])
