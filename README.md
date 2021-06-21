[![Snakemake](https://img.shields.io/badge/snakemake-≥6.0.2-brightgreen.svg)](https://snakemake.github.io)
[![Documentation Status](https://readthedocs.org/projects/ribofilio/badge/)](http://dammit.readthedocs.io/en/latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


# Ribofilio: A tool to measure ribosomal profiling dropoff rate

Ribofilio is a tool to measure ribosomal profiling drop rate that has been tested in ecoli and yeast so far. 

## Installation 

#### Manual Installation 

      git clone https://github.com/SherineAwad/ribofilio.git
      cd src 
      
## Dependencies

	screed

	numpy
	
	matplotlib
	
	sklearn 
        
	scipy


## How to run Ribofilio
 
To run ribofilio:


	python ribofilio.py --transcripts --footprint footprint.bed --rnaseq rnaseq.bed --binsize binsize --output output 
   
 
or simply:


	python ribofilio.py -t transcripts.bed -f footprint.bed -r rnaseq.bed  -b binsize -o output 


##### Parameters  



   ``--transcripts or -t for transcripts (required)`` 


   ``--footprint or -f for footprint bed file (required)`` 


   ``--rnaseq or -r for rnaseq bed file (if not available, dropoff rate won't be normalized by mRNA)`` 


   ``--subset or -s is a list of genes in file to run ribofilio on this subset only``


   ``--binsize or -b for binsize (default: 50)`` 


   ``--output or -o for output name`` 


   ``--ylogmin is the minimum y axis for log plots (default: -3)``


   ``--ylogmax is the maximum y axis for log plots (default: 2)``


### Example 

Running ribofilio on all gene: 

   
	python ribofilio.py -t yeast.fa -f SRR5945809.bed -r SRR5945808.bed 


Where yeast.fa is the transcripts, SRR5945809.bed is the bed file of footprints of sample, SRR5945808.bed is the mRNA bed file

here binsize used is 50 as no other binsize is passed. 

To run ribofilio on a subset of genes:
 

	python ribofilio.py -t yeast.fa -f SRR5945809.bed -r SRR5945808.bed -s subsetofgenes.txt 


Where subsetofgenes.txt is a list of genes: 

        YDL067C
   
        YGL187C
   
        YGL191W
   
        YHR051W
   
        YIL111W
   
        YLR038C
   
        YLR395C
   
        YMR256C
   
        YNL052W


#### Bed file format 

6 columns bed file format is required, a sample is as follows:

      YGL135W_mRNA    95      125     SRR5090936.1.1  1       +
      YKL009W_mRNA    70      98      SRR5090936.1.5  42      +
      YNL045W_mRNA    1664    1692    SRR5090936.1.11 40      +
      YNR050C_mRNA    32      62      SRR5090936.1.18 42      +
      YLR159C-A_mRNA  26      54      SRR5090936.1.20 0       -
      YHR133C_mRNA    443     473     SRR5090936.1.23 42      +

Refer to [Ensemble Bed format](https://m.ensembl.org/info/website/upload/bed.html) for more details regarding bed file formats.

#### Output 

Dropoff rate, dropoff rate per codon along with standard error, 95% confidence interval, root mean squared error (RMSE), R2 error, a t-test score of comparing the slope to a slope of zero and the corresponding pvalue will be print to both screen and output.regression.log file, and number of bins. If the plot mode is on, a weighted linear regression plot will be saved into output.Log.WLR.png. 

A sample of output.regression.log:

       Dropoff Dropoff per codon       RMSE    Rsquare SE      Margin Error    tscore  pvalue  No.of Bins
       -0.0051 -0.0003 0.0143  0.4907  0.0006  0.0011  -9.089848129268738      0.0     295

and output.Log.WLR.png:

<img src="https://github.com/SherineAwad/ribofilio/blob/master/notebooks/D1.Log.WLR.png" width="15" height="15">



## Documentation 

Looking for a more detailed tutorial: Take a look into [this tutorial](https://ribofilio.readthedocs.io/en/latest/protocol.html) for a complete Ribosomal profiling protocol using Ribofilio.


## Pipeline 

To run our pipeline on either our data in Datasets or your own data, use:

     snakemake --cores --use-conda 

