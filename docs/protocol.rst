========================================================================================
**Ribosomal Profiling Protocol**
========================================================================================

This is a ribosomal profiling protocol to estimate ribosomes' dropoff rate using Ribofilio.

Download Data
------------------
We will use a part of the GEO set  GSE102837::

    curl -O https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR5090934/SRR5090934.1
    curl -O https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR5090936/SRR5090936.1

Lets prepare our samples:: 

    fastq-dump --split-files SRR5090936.1
    fastq-dump --split-files SRR5090934.1
    mv SRR5090934.1_1.fastq  SRR5090934.fastq
    mv SRR5090936.1_1.fastq  SRR5090936.fastq

Trimming
-----------

Now, lets trim our data::
 
    mkdir fastqc 
    mkdir galore 
    trim_galore --gzip --retain_unpaired --trim1 -a "CTGTAGGCACCATCAAT" --fastqc --fastqc_args "--outdir fastqc" -o galore SRR5090936.fastq 
    trim_galore --gzip --retain_unpaired --trim1 -a "CTGTAGGCACCATCAAT" --fastqc --fastqc_args "--outdir fastqc" -o galore SRR5090934.fastq  

Prepare Reference
-------------------

We will need to align our samples, but we need to download and index our reference transcripts first:: 

    curl -O http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Saccharomyces_cerevisiae_Ensembl_R64-1-1.tar.gz
    gunzip Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
    mv  Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa yeast.fa

Then index our reference::

   bowtie2-build  yeast.fa yeast


Aligning 
-----------------

Now, we can align our samples:: 

   bowtie2 -x yeast -U galore/SRR5090936_trimmed.fq.gz  -S SRR5090936.bowtie2.sam
   samtools view -S -b SRR5090936.bowtie2.sam > SRR5090936.bowtie2.bam


And repeat for the second sample:: 

   bowtie2 -x yeast -U galore/SRR5090934_trimmed.fq.gz -S SRR5090934.bowtie2.sam
   samtools view -S -b SRR5090934.bowtie2.sam > SRR5090934.bowtie2.bam

Convert to Bed
-----------------

Now, we need to convert our aligned bam to bed format::

    bedtools bamtobed -i .SRR5090936.bowtie2.bam > SRR5090936.bed 
    bedtools bamtobed -i .SRR5090934.bowtie2.bam > SRR5090934.bed

Run Ribofilio 
-------------------
Now, we are all set to run ribofilio, we are here using the default parameters::

    python ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed  -o SRR5090936_SRR5090934 


Let's change the binsize from 50 (default) to 100:: 

    python ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -b 100 -o SRR5090936_SRR5090934


We can also change the default indeces::

   python ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -b 100 --ylogmin -5 --ylogmax 5 -o SRR5090936_SRR5090934


To run Ribofilio without plots:: 

   python ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed --plot 0 -o SRR5090936_SRR5090934 

We can also run Ribofilio on a subset of genes. Let's say we have this file, GO0016458.txt which has two genes for this GO0016458:: 

    YOR140W
    YBL079W

We can run Ribofilio on this subset only::

    python ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -s GO0016458.txt -b 100 -o SRR5090936_SRR5090934_subset 


Run Ribofilio without mRNA normalization
--------------------------------------------

By default, Ribofilio normalize the dropoff rate of ribosomes' foot print using the corresponding mRNA reads. To estimate ribosomes' dropoff rate without mRNA normalization::

   python ribofilio.py -t yeast.fa -f SRR5090936.bed -b 100 -o SRR5090936_SRR5090934 


