
yeast.fa:
	curl -O ftp://ftp.ensembl.org/pub/release-98/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
	gunzip Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
	mv  Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa yeast.fa

#Dataset GSE91068
D3FP.fastq:
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR5090934/SRR5090934.1
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR5090935/SRR5090935.1
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR5090936/SRR5090936.1
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR5090937/SRR5090937.1
	fastq-dump --split-files SRR5090934.1	
	fastq-dump --split-files SRR5090935.1
	fastq-dump --split-files SRR5090936.1
	fastq-dump --split-files SRR5090937.1
	mv SRR5090934.1_1.fastq  D1mRNA.fastq
	mv SRR5090935.1_1.fastq  D3mRNA.fastq
	mv SRR5090936.1_1.fastq  D1FP.fastq
	mv SRR5090937.1_1.fastq  D3FP.fastq

#Dataset GSE134152
D10FP.fastq:
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-21/SRR9670820/SRR9670820.1
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-19/SRR9670821/SRR9670821.1
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-21/SRR9670822/SRR9670822.1
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-21/SRR9670823/SRR9670823.1
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-20/SRR9670814/SRR9670814.1
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-20/SRR9670815/SRR9670815.1
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-19/SRR9670816/SRR9670816.1 
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-20/SRR9670817/SRR9670817.1 
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-19/SRR9670818/SRR9670818.1
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-19/SRR9670819/SRR9670819.1
	fastq-dump --split-files SRR9670820
	fastq-dump --split-files SRR9670821
	fastq-dump --split-files SRR9670822
	fastq-dump --split-files SRR9670823
	fastq-dump --split-files SRR9670814
	fastq-dump --split-files SRR9670815
	fastq-dump --split-files SRR9670816
	fastq-dump --split-files SRR9670817
	fastq-dump --split-files SRR9670818
	fastq-dump --split-files SRR9670819 
	mv SRR9670820_1.fastq D2mRNA.fastq
	mv SRR9670821_1.fastq D4mRNA.fastq
	mv SRR9670822_1.fastq D2FP.fastq
	mv SRR9670823_1.fastq D4FP.fastq
	mv SRR9670814_1.fastq D9mRNA.fastq 
	mv SRR9670815_1.fastq D10mRNA.fastq 
	cat SRR9670816_1.fastq SRR9670817_1.fastq > D9FP.fastq 
	cat SRR9670818_1.fastq SRR9670819_1.fastq > D10FP.fastq  	

#Dataset GSE13750
D8FP.fastq: 
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR014376/SRR014376.2 	
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR014375/SRR014375.3
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR014374/SRR014374.3
	fastq-dump --split-files SRR014376
	fastq-dump --split-files SRR014375
	fastq-dump --split-files SRR014374
	cat SRR014376_1.fastq SRR014375_1.fastq SRR014374_1.fastq > D5FP.fastq
 
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR014377/SRR014377.2
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR014378/SRR014378.2
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR014379/SRR014379.2
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR014380/SRR014380.2
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR014381/SRR014381.2
	fastq-dump --split-files SRR014377
	fastq-dump --split-files SRR014378
	fastq-dump --split-files SRR014379
	fastq-dump --split-files SRR014380
	fastq-dump --split-files SRR014381
	cat SRR014377_1.fastq SRR014381_1.fastq SRR014380_1.fastq SRR014379_1.fastq SRR014378_1.fastq > D6FP.fastq  
	
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR014369/SRR014369.3
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR014368/SRR014368.3
	fastq-dump --split-files SRR014369
	fastq-dump --split-files SRR014368
	cat SRR014369_1.fastq SRR014368_1.fastq > D7FP.fastq  

	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR014373/SRR014373.2
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR014372/SRR014372.2
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR014371/SRR014371.3
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR014370/SRR014370.2
	fastq-dump --split-files SRR014373
	fastq-dump --split-files SRR014372
	fastq-dump --split-files SRR014371
	fastq-dump --split-files SRR014370
	cat SRR014373_1.fastq SRR014372_1.fastq SRR014371_1.fastq SRR014370_1.fastq > D8FP.fastq 

	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR014385/SRR014385.3
	fastq-dump --split-files SRR014385
	mv SRR014385_1.fastq D5mRNA.fastq 

	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR028774/SRR028774.3
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR014387/SRR014387.3
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR014386/SRR014386.2
	fastq-dump --split-files SRR028774
	fastq-dump --split-files SRR014387
	fastq-dump --split-files SRR014386
	cat SRR028774_1.fastq SRR014387_1.fastq SRR014386_1.fastq > D6mRNA.fastq 
	
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR014382/SRR014382.3
	fastq-dump --split-files SRR014382
	mv SRR014382_1.fastq D7mRNA.fastq 

	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR014384/SRR014384.2
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR014383/SRR014383.2
	fastq-dump --split-files SRR014384
	fastq-dump --split-files SRR014383
	cat SRR014384_1.fastq SRR014383_1.fastq > D8mRNA.fastq 

clean: 
	rm SRR014368_1.fastq
	rm SRR014369_1.fastq
	rm SRR014370_1.fastq
	rm SRR014371_1.fastq
	rm SRR014372_1.fastq
	rm SRR014373_1.fastq
	rm SRR014374_1.fastq
	rm SRR014375_1.fastq
	rm SRR014376_1.fastq
	rm SRR014377_1.fastq
	rm SRR014378_1.fastq
	rm SRR014379_1.fastq
	rm SRR014380_1.fastq
	rm SRR014381_1.fastq
	rm SRR014383_1.fastq
	rm SRR014384_1.fastq
	rm SRR014386_1.fastq
	rm SRR014387_1.fastq
	rm SRR028774_1.fastq
	rm SRR9670816_1.fastq
	rm SRR9670817_1.fastq
	rm SRR9670818_1.fastq
	rm SRR9670819_1.fastq	
