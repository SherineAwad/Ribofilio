
# Ribofilio: A tool to measure ribosomal profiling dropoff rate

Ribofilio is a tool to measure ribosomal profiling drop rate that has been tested in ecoli and yeast so far.


## Dependencies

	screed

	numpy
	
	matplotlib
	
	sklearn 
   
        scipy

## Tutorial

Take a look into [this tutorial](https://ribofilio.readthedocs.io/en/latest/protocol.html) for a complete Ribosomal profiling protocol using Ribofilio.


## How to run Ribofilio
 
To run ribofilio:


	python ribofilio.py --transcripts --footprint footprint.bed --rnaseq rnaseq.bed --binsize binsize --output output 
   
 
or simply:


	python ribofilio.py -t transcripts -f footprint.bed -r rnaseq.bed  -b binsize -o output 

 
Where: 


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


