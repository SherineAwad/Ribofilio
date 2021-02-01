import pytest
import screed
import os
import sys 

path = os.getcwd() 
path  = os.path.join(path,"src") 
sys.path.append(path)

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
