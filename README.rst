=================================================================
**Ribofilio: A tool to measure ribosomal profiling drop rate**
=================================================================

Ribofilio is a tool to measure ribosomal profiling drop rate that has been tested in ecoli and yeast so far.


Dependecies: 
       
       python 3.7 


       numpy


       screed 


       And for plotting: 


       matplotlib 


       scikit-learn  


To Install Ribofilio::


   git clone https://github.com/SherineAwad/ribofilio.git


To run ribofilio::


    python ribofilio.py --transcript --footprint footprint_sample --mRNA mRNA_sample --binsize (default is 50) --cutoff (default is 0)
    
or simply::


    python ribofilio.py -t transcripts -f -r  -b  -c 

 
Where: 


   ``--transcript or -t for transcripts (required)`` 


   ``--footprint or -f for footprint bed file (required)`` 


   ``--mRNA or -r for mRNA bed file (if not available, droprate won't be normalized by mRNA)`` 


   ``--subset or -s is a list of genes in file to run ribofilio on this subset only``


   ``--binsize or -b for binsize (default: 50)`` 


   ``--cutoff or -c  for cutoff or minimum number of genes required to contribute to a position to be counted (default: 0)``


   ``--ymin is the minimum y axis for linear plots (default: 0)`` 


   ``--ymax is the maximum y axis for linear plots (default: 30)``


   ``--ylogmin is the minimum y axis for log plots (default: -5)``


   ``--ylogmax is the maximum y axis for log plots (default: 5)``


Example 
########

We will create two folders for galore and fastq:: 
   
    python ribofilio.py -t yeast.fa -f SRR5945809.bed -r SRR5945808.bed -b 50 -c 50 

Where yeast.fa is the transcripts, is the bed file of footprints of sample SRR5945809.bed, SRR5945808.bed is the mRNA bed file, binsize is 50 and no cutoff is 50 which means
at least 50 genes should contribute to the reads in a position to be considered in bins. 



To run ribofilio on a subset of genes:: 


    python ribofilio.py -t yeast.fa -f SRR5945809.bed -r SRR5945808.bed -b 50 -c 50 -s subsetofgenes.txt 
