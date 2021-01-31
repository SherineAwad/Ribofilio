import pytest
import screed

def test_get_transcipts():
    assert get_transcripts("/rds/project/yhbl2/rds-yhbl2-genehunter/SM/RibosomalProfiling/ribofilio/tests/transctipts_sample.fa") == {"YBR024W":"ATGTTGAATAGTTCAA"}

def test_get_genes():
    max_gene_length, gene_length = get_genes("/rds/project/yhbl2/rds-yhbl2-genehunter/SM/RibosomalProfiling/ribofilio/tests/transctipts_sample.fa")
    assert (max_gene_length == 16 ) 
    assert (gene_length ==  {"YBR024W":16}) 

def test_get_reads():
    coverage = get_reads("/rds/project/yhbl2/rds-yhbl2-genehunter/SM/RibosomalProfiling/ribofilio/tests/sample.bed",  {"YKL152C":50}) 
    assert (coverage == {"YKL152C":[93]}) 


def get_transcripts(transcripts_file):
    transcripts = {}
    for record in screed.open(transcripts_file):
        gene_name = record.name.split(" ")[0]
        if "mRNA" in gene_name:
            gene_name = gene_name.split("_")[0]
        transcripts[gene_name] = record.sequence
    return transcripts





def get_genes(transcripts_file):
    max_gene_length = -100
    genes_length = {}
    for record in screed.open(transcripts_file):
        gname = record.name.split(' ')[0]
        if "mRNA" in gname:
            gname=gname.split('_')[0]
        genes_length[str(gname)] = int(len(record.sequence))
        if genes_length[gname] > max_gene_length:
            max_gene_length = genes_length[gname]
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
        coverage[str(gene_name)].append(int(record[2]) )
    return coverage

