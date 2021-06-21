
yeast.fa:
	curl -O ftp://ftp.ensembl.org/pub/release-98/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
	gunzip Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
	mv  Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa yeast.fa

#Dataset GSE91068
SRR5090937.fastq:
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR5090934/SRR5090934.1
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR5090935/SRR5090935.1
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR5090936/SRR5090936.1
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR5090937/SRR5090937.1
	fastq-dump --split-files SRR5090934.1	
	fastq-dump --split-files SRR5090935.1
	fastq-dump --split-files SRR5090936.1
	fastq-dump --split-files SRR5090937.1
	mv SRR5090934.1_1.fastq  SRR5090934.fastq
	mv SRR5090935.1_1.fastq  SRR5090935.fastq
	mv SRR5090936.1_1.fastq  SRR5090936.fastq
	mv SRR5090937.1_1.fastq  SRR5090937.fastq
	rm SRR5090934.1
	rm SRR5090935.1
	rm SRR5090936.1
	rm SRR5090937.1

#Dataset GSE134152
SRR9670823.fastq:
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-21/SRR9670820/SRR9670820.1
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-19/SRR9670821/SRR9670821.1
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-21/SRR9670822/SRR9670822.1
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-21/SRR9670823/SRR9670823.1
	fastq-dump --split-files SRR9670820
	fastq-dump --split-files SRR9670821
	fastq-dump --split-files SRR9670822
	fastq-dump --split-files SRR9670823
	mv SRR9670820_1.fastq SRR9670820.fastq
	mv SRR9670821_1.fastq SRR9670821.fastq
	mv SRR9670822_1.fastq SRR9670822.fastq
	mv SRR9670823_1.fastq SRR9670823.fastq
clean: 
	rm SRR5090934.1
	rm SRR5090935.1
	rm SRR5090936.1
	rm SRR5090937.1
	rm SRR9670820.1
	rm SRR9670821.1
	rm SRR9670822.1
	rm SRR9670823.1
