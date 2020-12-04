

#Dataset GSE102837
GSE102837/SRR5945808.fastq.gz:
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR594/009/SRR5945809/SRR5945809.fastq.gz
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR594/008/SRR5945808/SRR5945808.fastq.gz
	mv SRR5945809.fastq.gz GSE102837
	mv SRR5945808.fastq.gz GSE102837

#Dataset GSE91068
GSE91068/SRR5090937.fastq:
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR5090934/SRR5090934.1
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR5090935/SRR5090935.1
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR5090936/SRR5090936.1
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR5090937/SRR5090937.1
	fastq-dump --split-files SRR5090934.1
	fastq-dump --split-files SRR5090935.1
	fastq-dump --split-files SRR5090936.1
	fastq-dump --split-files SRR5090937.1
	mv SRR5090934.1_1.fastq  GSE91068/SRR5090934.fastq
	mv SRR5090935.1_1.fastq  GSE91068/SRR5090935.fastq
	mv SRR5090936.1_1.fastq  GSE91068/SRR5090936.fastq
	mv SRR5090937.1_1.fastq  GSE91068/SRR5090937.fastq


#Dataset GSE134152
GSE134152/SRR9670816.fastq:
	wget https://sra-download.ncbi.nlm.nih.gov/traces/sra53/SRR/009444/SRR9670817
	wget https://sra-download.ncbi.nlm.nih.gov/traces/sra59/SRR/009444/SRR9670821
	wget https://sra-download.ncbi.nlm.nih.gov/traces/sra67/SRR/009444/SRR9670822
	wget https://sra-download.ncbi.nlm.nih.gov/traces/sra25/SRR/009444/SRR9670814
	wget https://sra-download.ncbi.nlm.nih.gov/traces/sra59/SRR/009444/SRR9670819
	wget https://sra-download.ncbi.nlm.nih.gov/traces/sra50/SRR/009444/SRR9670820
	wget https://sra-download.ncbi.nlm.nih.gov/traces/sra53/SRR/009444/SRR9670818
	wget https://sra-download.ncbi.nlm.nih.gov/traces/sra2/SRR/009444/SRR9670815
	wget https://sra-download.ncbi.nlm.nih.gov/traces/sra2/SRR/009444/SRR9670823
	wget https://sra-download.ncbi.nlm.nih.gov/traces/sra19/SRR/009444/SRR9670816
	fastq-dump --split-files SRR9670817
	fastq-dump --split-files SRR9670821
	fastq-dump --split-files SRR9670822
	fastq-dump --split-files SRR9670820
	fastq-dump --split-files SRR9670818
	fastq-dump --split-files SRR9670815
	fastq-dump --split-files SRR9670823
	fastq-dump --split-files SRR9670816

#Dataset GSE81966
GSE81966/SRR3591816.fastq:  
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-7/SRR3591816/SRR3591816.1
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-7/SRR3591812/SRR3591812.1
	there is more, checking this for now
	fastq-dump --split-files SRR3591816.1
	fastq-dump --split-files SRR3591812.1
	mv SRR3591816.1_1.fastq GSE81966/SRR3591816.fastq
	mv SRR3591812.1_1.fastq GSE81966/SRR3591812.fastq

clean:
	rm SRR5090934.1  
	rm SRR5090935.1  
	rm SRR5090936.1  
	rm SRR5090937.1	
	rm SRR3591812.1  
	rm SRR3591816.1
