

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
	fasterq-dump --split-files SRR5090934.1
	fasterq-dump --split-files SRR5090935.1
	fasterq-dump --split-files SRR5090936.1
	fasterq-dump --split-files SRR5090937.1
	mv SRR5090934.1.fastq  GSE91068/SRR5090934.fastq
	mv SRR5090935.1.fastq  GSE91068/SRR5090935.fastq
	mv SRR5090936.1.fastq  GSE91068/SRR5090936.fastq
	mv SRR5090937.1.fastq  GSE91068/SRR5090937.fastq


