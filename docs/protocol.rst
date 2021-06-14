========================================================================================
**Ribosomal profiling protocol to estimate ribosomes dropoff rate using Ribofilio**
========================================================================================


Download Data
------------------
We will use a part of the GEO set  GSE102837:

    wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR5090934/SRR5090934.1
    wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR5090936/SRR5090936.1

Lets prepare our samples: 

    fastq-dump --split-files SRR5090936.1
    fastq-dump --split-files SRR5090934.1
    mv SRR5090934.1_1.fastq  SRR5090934.fastq
    mv SRR5090936.1_1.fastq  SRR5090936.fastq

Trimming
-----------

Now, lets trim our data: 
 
    mkdir fastqc 
    mkdir galore 
    trim_galore --gzip --retain_unpaired --trim1 -a "CTGTAGGCACCATCAAT" --fastqc --fastqc_args "--outdir fastqc" -o galore SRR5090936.fastq 
    trim_galore --gzip --retain_unpaired --trim1 -a "CTGTAGGCACCATCAAT" --fastqc --fastqc_args "--outdir fastqc" -o galore SRR5090934.fastq  

Prepare Reference
-------------------

We will need to align our samples, but we need to download and index our reference transcripts first: 


    wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Saccharomyces_cerevisiae_Ensembl_R64-1-1.tar.gz
    gunzip Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
    mv  Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa yeast.fa

Then index our reference: 

    
   bowtie2-build  yeast.fa yeast


Aligning 
-----------------

Now, we can align our samples: 

   bowtie2 -x yeast -U galore/SRR5090936_trimmed.fq.gz  -S SRR5090936.bowtie2.sam
   samtools view -S -b SRR5090936.bowtie2.sam > SRR5090936.bowtie2.bam


And repeat for the second sample: 

   bowtie2 -x yeast -U galore/SRR5090934_trimmed.fq.gz -S SRR5090934.bowtie2.sam
   samtools view -S -b SRR5090934.bowtie2.sam > SRR5090934.bowtie2.bam

Convert to Bed
-----------------

Now, we need to convert our aligned bam to bed format: 


    bedtools bamtobed -i .SRR5090936.bowtie2.bam > SRR5090936.bed 
    bedtools bamtobed -i .SRR5090934.bowtie2.bam > SRR5090934.bed

Run Ribofilio 
-------------------
Now, we are all set to run ribofilio: 


    python ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed  -b 50 -o SRR5090936_SRR5090934 

 
