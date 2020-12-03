BBMAP = /home/sa941/software/bbmap

YEAST = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/annotations/Saccharomyces_cerevisiae/Ensembl/R64-1-1

DATASETS = /home/ec2-user/ribofilio/datasets
SCRIPTS = /home/ec2-user/ribofilio/scripts
CURRENT  =  /home/ec2-user/ribofilio/pipeline

subset1 := $(shell seq 1 500)
subset2 := $(shell seq 501 1000)


#--------------------------------------------------------------------------------------------------------------------------------
GO = GO_0006119 GO_0006406 GO_0006412 GO_0006950 GO_0007049 GO_0009651 GO_0016458 GO_0031047 GO_0031990 GO_0042254
CHRs = chrI chrII chrIII chrIV chrV chrVI chrVII chrVIII chrVIII chrIX chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI
#--------------------------------------------------------------------------
#Get yeasts Transcripts and prepare index and bed files
#--------------------------------------------------------------------------
Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa:
	curl -O ftp://ftp.ensembl.org/pub/release-98/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
	gunzip Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
	mv  Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa yeast.fa

Saccharomyces_cerevisiae_Ensembl_R64-1-1.tar.gz:
	wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Saccharomyces_cerevisiae_Ensembl_R64-1-1.tar.gz

yeast.bed:
	gtf2bed < ${YEAST}/Annotation/Genes/genes.gtf > yeast.bed
yeast:
	bowtie2-build  yeast.fa yeast


#----------------------------------------------------------------------------------------------------------------------
#Get GO, GC, orthogolous, and Chromosomes from Ensembl Biomart required for all GO, orthogolus, and Chr experiments
#----------------------------------------------------------------------------------------------------------------------
yeastchrXVI.txt:
	Rscript getBiomart.R  #This script will download the GO, Chromosomes, GC, and orthogolus data 

#---------------------------------------------------------------------------------------------------
#Prepare Salmon index for running Salmon required fo  TPM experiments 
#-----------------------------------------------------------------------------------------------------
SalmonYeastIndex:
	salmon index --index SalmonYeastIndex --transcripts yeast.fa


#--------------------------------------------------------------------------------------------------------------------------------
SRR9670816.fastq: #get GSE134152  and GSE91068
	bash getSRA.sh 

#---------------------------------------------------------------------------
#Dataset: GSE102837 
#---------------------------------------------------------------------------

#GSM2747049:SRR5945809
#GSM2747048:SRR5945808

galore/SRR5945808.trimmed.fastq.gz:
	mkdir ${CURRENT}/galore
	mkdir ${CURRENT}/fastqc 	
	trim_galore --gzip --retain_unpaired --trim1 -a "CTGTAGGCACCATCAATAGATCGGAAGAGCACACGTCTGAACTCCAGTCA" --fastqc --fastqc_args "--outdir ${CURRENT}/fastqc" -o ${CURRENT}/galore ${DATASETS}/GSE102837/SRR5945809.fastq.gz
	trim_galore --gzip --retain_unpaired --trim1 -a "CTGTAGGCACCATCAATAGATCGGAAGAGCACACGTCTGAACTCCAGTCA" --fastqc --fastqc_args "--outdir ${CURRENT}/fastqc" -o ${CURRENT}/galore ${DATASETS}/GSE102837/SRR5945808.fastq.gz


SRR5945808.bowtie2.bam: 
	bowtie2 -x yeast -U galore/SRR5945809_trimmed.fq.gz  -S SRR5945809.bowtie2.sam
	samtools view -S -b SRR5945809.bowtie2.sam > SRR5945809.bowtie2.bam
	bowtie2 -x yeast -U galore/SRR5945808_trimmed.fq.gz  -S SRR5945808.bowtie2.sam
	samtools view -S -b SRR5945808.bowtie2.sam > SRR5945808.bowtie2.bam


SRR5945808.bed: SRR5945809.bowtie2.bam SRR5945808.bowtie2.bam 
	bedtools bamtobed -i SRR5945809.bowtie2.bam > SRR5945809.bed 
	bedtools bamtobed -i SRR5945808.bowtie2.bam > SRR5945808.bed 


#Prepare simulated shuffled Footprint bed files
SRR5945809_1000.bed: SRR5945809.bed 
	 $(foreach i, $(subset1), python ${SCRIPTS}/shuffle.py -t yeast.fa -b SRR5945809.bed  -o SRR5945809_$(i).bed;) 
	 $(foreach i, $(subset2), python ${SCRIPTS}/shuffle.py -t yeast.fa -b SRR5945809.bed  -o SRR5945809_$(i).bed;) 

#Run ribofilio on main set and on shuffled set 
#Store results in  SRR5945809_regfile.log 
#run stats.py on  SRR5945809_regfile.log to calculate pvalues 
SRR5945809_regfile.log: SRR5945809_1000.bed 
	python ${SCRIPTS}/ribofilioWLR.py -t yeast.fa -f SRR5945809.bed -r SRR5945808.bed -b 50 -n SRR5945809_regfile.log 
	$(foreach i, $(subset1), python ${SCRIPTS}/ribofilioWLR.py -t yeast.fa -f SRR5945809_$(i).bed -r SRR5945808.bed -b 50 --plots 0 -n SRR5945809_regfile.log;)
	$(foreach i, $(subset2), python ${SCRIPTS}/ribofilioWLR.py -t yeast.fa -f SRR5945809_$(i).bed -r SRR5945808.bed -b 50 --plots 0 -n SRR5945809_regfile.log;)

SRR5945809_Bregfile.log: SRR5945809_1000.bed
	python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5945809.bed -r SRR5945808.bed -b 50 -n SRR5945809_Bregfile.log
	$(foreach i, $(subset1), python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5945809_$(i).bed -r SRR5945808.bed -b 50 --plots 0 -n SRR5945809_Bregfile.log;)
	$(foreach i, $(subset2), python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5945809_$(i).bed -r SRR5945808.bed -b 50 --plots 0 -n SRR5945809_Bregfile.log;)	


#Run ribofilio on GO 
SRR5945809_GO0006950_regfile.log: SRR5945809_1000.bed yeastchrXVI.txt 
	python ${SCRIPTS}/ribofilioWLR.py -t yeast.fa -f SRR5945809.bed -r SRR5945808.bed -b 50 -s GO_0006950.txt  -n SRR5945809_GO0006950_regfile.log 
	$(foreach i, $(subset1),python ${SCRIPTS}/ribofilioWLR.py -t yeast.fa -f SRR5945809_$(i).bed -r SRR5945808.bed -b 50 --plots 0 -s GO_0006950.txt -n  SRR5945809_GO0006950_regfile.log;) 
	$(foreach i, $(subset2),python ${SCRIPTS}/ribofilioWLR.py -t yeast.fa -f SRR5945809_$(i).bed -r SRR5945808.bed -b 50 --plots 0 -s GO_0006950.txt -n  SRR5945809_GO0006950_regfile.log;)	
SRR5945809_GO0006119_regfile.log: SRR5945809_1000.bed yeastchrXVI.txt
	python ${SCRIPTS}/ribofilioWLR.py -t yeast.fa -f SRR5945809.bed -r SRR5945808.bed -b 50 -s GO_0006119.txt  -n SRR5945809_GO0006119_regfile.log
	$(foreach i, $(subset1),python ${SCRIPTS}/ribofilioWLR.py -t yeast.fa -f SRR5945809_$(i).bed -r SRR5945808.bed -b 50 --plots 0 -s GO_0006119.txt -n SRR5945809_GO0006119_regfile.log;)
	$(foreach i, $(subset2),python ${SCRIPTS}/ribofilioWLR.py -t yeast.fa -f SRR5945809_$(i).bed -r SRR5945808.bed -b 50 --plots 0 -s GO_0006119.txt -n SRR5945809_GO0006119_regfile.log;)




#------------------------------------------------------------------------------------------
#Case Study:GSE91068
#------------------------------------------------------------------------------------------
s31 = SRR5090934 #GSM2420486
s32 = SRR5090935 #GSM2420487
s33 = SRR5090936 #GSM2420488
s34 = SRR5090937 #GSM2420489

SAMPLES3 = ${s31} ${s32} ${s33} ${s34}

ribofilio/datasets/SRR5090937.fastq:
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR5090934/SRR5090934.1
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR5090935/SRR5090935.1
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR5090936/SRR5090936.1
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR5090937/SRR5090937.1
	fasterq-dump --split-files SRR5090934.1
	fasterq-dump --split-files SRR5090935.1
	fasterq-dump --split-files SRR5090936.1
	fasterq-dump --split-files SRR5090937.1
	mv SRR5090934.1.fastq  SRR5090934.fastq 
	mv SRR5090935.1.fastq  SRR5090935.fastq
	mv SRR5090936.1.fastq  SRR5090936.fastq
	mv SRR5090937.1.fastq  SRR5090937.fastq

#Adpaters reported by dataset author @jiashun@genome.ucsf.edu 
SRR5090937_trimmed.fastq.gz: ${DATASETS}/SRR5090937.fastq
	$(foreach i, $(SAMPLES3), trim_galore --gzip --retain_unpaired --trim1 -a "CTGTAGGCACCATCAAT" --fastqc --fastqc_args "--outdir ${CURRENT}/fastqc" -o ${CURRENT}/galore  ${CURRENT}/$(i).fastq ; )

SRR5090937.bowtie2.bam: SRR5090937_trimmed.fastq.gz
	$(foreach i, $(SAMPLES3), bowtie2 -x yeast -U galore/$(i)_trimmed.fq.gz  -S $(i).bowtie2.sam;)
	$(foreach i, $(SAMPLES3), samtools view -S -b $(i).bowtie2.sam > $(i).bowtie2.bam;)

SRR5090937.bed: SRR5090937.bowtie2.bam
	$(foreach i, $(SAMPLES3), bedtools bamtobed -i $(i).bowtie2.bam > $(i).bed;)

SRR5090936_SRR5090934.yeastGenes15000_50_0.LR.png:
	$(foreach i, $(LEN), python ${SCRIPTS}/ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed -b 50 -s yeastGenes$(i).txt;)
	$(foreach i, $(LEN), python ${SCRIPTS}/ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -b 50 -s yeastGenes$(i).txt;)

SRR5090936_SRR5090934.200.LR.png:  
	$(foreach i, $(CUTOFFS), python ${SCRIPTS}/ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed -b 50 -c ${i};)
	$(foreach i, $(CUTOFFS), python ${SCRIPTS}/ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -b 50 -c ${i};)

SRR5090936_SRR5090934_50_0_200.LR.png:
	$(foreach i, $(Ns), python ${SCRIPTS}/ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed -b 50 -e $(i);)
	$(foreach i, $(Ns), python ${SCRIPTS}/ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -b 50 -e $(i);)

#Get TPM and high and low expressed genes  
SRR5090934_l50.tpm: #SalmonYeastIndex
	salmon quant -i SalmonYeastIndex --libType U  -r ${CURRENT}/galore/SRR5090935_trimmed.fq.gz \
	-o SalmonSRR5090935 -p 3 --writeUnmappedNames --seqBias --gcBias --validateMappings
	@cmd1 = `cut -f1,4 SalmonSRR5090935/quant.sf > SRR5090935.tpm` ; echo $$cmd1
	@cmd2 = `awk '$$2 >50 {print $$0}' SRR5090935.tpm > SRR5090935_g50.tpm` ; echo $$cmd2
	@cmd3 = `awk '$$2 <50 {print $$0}' SRR5090935.tpm > SRR5090935_l50.tpm` ; echo $$cmd3
	salmon quant -i SalmonYeastIndex --libType U  -r ${CURRENT}/galore/SRR5090934_trimmed.fq.gz \
	-o SalmonSRR5090934 -p 3 --writeUnmappedNames --seqBias --gcBias --validateMappings
	@cmd1 = `cut -f1,4 SalmonSRR5090934/quant.sf > SRR5090934.tpm` ; echo $$cmd1
	@cmd2 = `awk '$$2 >50 {print $$0}' SRR5090934.tpm > SRR5090934_g50.tpm` ; echo $$cmd2
	@cmd3 = `awk '$$2 <50 {print $$0}' SRR5090934.tpm > SRR5090934_l50.tpm` ; echo $$cmd3

#Run ribofilio on a subset of genes: high and low expressed genes 
SRR5090936.SRR5090934.l50tpm.LogLinear.png:  
	python ${SCRIPTS}/ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed -s SRR5090935_g50.tpm -b 50
	python ${SCRIPTS}/ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed -s SRR5090935_l50.tpm -b 50
	python ${SCRIPTS}/ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -s SRR5090934_g50.tpm -b 50
	python ${SCRIPTS}/ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -s SRR5090934_l50.tpm -b 50

#Run ribofilio on genes clustered by GO  
SRR5090936.SRR5090934.GO.LR.png: 
	$(foreach i, $(GO),python ${SCRIPTS}/ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed -b 50 -s ${i}.txt ;)
	$(foreach i, $(GO),python ${SCRIPTS}/ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -b 50 -s ${i}.txt ;)


#Run ribofilio on genes clustered by Chromosomes
SRR5090936_SRR5090934.yeastchrXVI_50_0.LR.png:
	$(foreach i, $(CHRs),python ${SCRIPTS}/ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed -b 50 --subset yeast${i}.txt;)
	$(foreach i, $(CHRs),python ${SCRIPTS}/ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -b 50 --subset yeast${i}.txt;)

#Run ribofilio on orthogolus genes to human and mouse 
SRR5090936_SRR5090934.yeastorthuman_50_0.png:  
	python ${SCRIPTS}/ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed  -b 50 --subset yeastorthuman.txt
	python ${SCRIPTS}/ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed  -b 50 --subset yeastortmouse.txt
	python ${SCRIPTS}/ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed  -b 50 --subset yeastorthuman.txt
	python ${SCRIPTS}/ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed  -b 50 --subset yeastortmouse.txt

#Generate vcf file 
SRR5090934.vcf: SRR5090935.bowtie2.bam
	samtools sort -T SRR5090935.tmp.sorted -o SRR5090935.bowtie2.sorted.bam SRR5090935.bowtie2.bam
	samtools index SRR5090935.bowtie2.sorted.bam SRR5090935.bowtie2.sorted.bam.bai
	samtools mpileup -uD -s -S -f yeast.fa SRR5090935.bowtie2.sorted.bam | bcftools  call -m -v -Ob > SRR5090935.bcf
	bcftools view SRR5090935.bcf > SRR5090935.vcf 
	samtools sort -T SRR5090934.tmp.sorted -o SRR5090934.bowtie2.sorted.bam SRR5090934.bowtie2.bam
	samtools index SRR5090934.bowtie2.sorted.bam SRR5090934.bowtie2.sorted.bam.bai
	samtools mpileup -uD -s -S -f yeast.fa SRR5090934.bowtie2.sorted.bam | bcftools  call -m -v -Ob > SRR5090934.bcf
	bcftools view SRR5090934.bcf > SRR5090934.vcf 

#Prepare subset file one with genes with zero snps and one with genes with at least one snp
#outputs genes list with zero snps and genes with at least one snp
SRR5090934_1.snps:
	python countsnps.py  yeast.fa SRR5090935.vcf --cutoff 1
	python countsnps.py  yeast.fa SRR5090934.vcf --cutoff 1

#Run ribofilio on genes with zero snps and one or more snp

SRR5090936_SRR5090934.SRR5090934_1_50_0.LR.png:
	python ${SCRIPTS}/ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed  -b 50 --subset SRR5090935_zerosnps.txt 
	python ${SCRIPTS}/ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed  -b 50 --subset SRR5090935_1snps.txt
	python ${SCRIPTS}/ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed  -b 50 --subset SRR5090934_zerosnps.txt
	python ${SCRIPTS}/ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed  -b 50 --subset SRR5090934_1snps.txt


