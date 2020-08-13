BBMAP = /home/sa941/software/bbmap

YEAST = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/annotations/Saccharomyces_cerevisiae/Ensembl/R64-1-1


WHERE  = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/RibosomalProfiling
#https://www.ebi.ac.uk/ena/data/view/SRR8360209 


#--------------------------------------------------------------------------------------------------------------------------------
GO = GO_0006119 GO_0006406 GO_0006412 GO_0006950 GO_0007049 GO_0009651 GO_0016458 GO_0023052 GO_0031047 GO_0031990 GO_0042254
CHRs = chrI chrII chrIII chrIV chrV chrVI chrVII chrVIII chrVIII chrIX chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI
CUTOFFS = 0 100 200 500 1000 2000 
#--------------------------------------------------------------------------------------------------------------------------------


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

yeast.bowtie:
	bowtie-build yeast.fa yeast.bowtie

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
#--------------------------------------------------------------------------------------------------------------------------------

#--------------------------
#Case Study:GSE134152 
#--------------------------

s21 = SRR9670814 #mRNA SD sample 1  
s22 = SRR9670815 #mRNA GR sample 2 
s23 = SRR9670816 #footprint SD sample 3-1
s24 = SRR9670817 #footprint SD sample 3-2
s25 = SRR9670818 #footprint GR sample 4-1
s26 = SRR9670819 #footprint GR sample 4-2 
s27 = SRR9670820  #mRNA SD sample SD 5 
s28 = SRR9670821  #mRNA sample GR  6
s29 = SRR9670822  #footprint SD  sample 7 
s30 = SRR9670823  #footprint GR sample 8 


SAMPLES2 = ${s21} ${s22} ${s23} ${s24} ${s25} ${s26} ${s27} ${s28} ${s29} ${s30} 
galore/SRR9670823_trimmed.fastq.gz:
	$(foreach i, $(SAMPLES2), trim_galore --gzip --retain_unpaired --trim1 -a "CTGTAGGCACCATCAAT" --fastqc --fastqc_args "--outdir ${WHERE}/fastqc" -o ${WHERE}/galore  ${WHERE}/$(i).fastq ; )
	
SRR9670823.bowtie2.bam: galore/SRR9670823_trimmed.fastq.gz 
	$(foreach i, $(SAMPLES2), bowtie2 -x yeast -U galore/$(i)_trimmed.fq.gz  -S $(i).bowtie2.sam;)
	$(foreach i, $(SAMPLES2), samtools view -S -b $(i).bowtie2.sam > $(i).bowtie2.bam;)

SRR9670823.bed: SRR9670823.bowtie2.bam 
	$(foreach i, $(SAMPLES2), bedtools bamtobed -i $(i).bowtie2.bam > $(i).bed;) 
	
gr4.bowtie2.bam:
	samtools merge sd3.bowtie2.bam SRR9670816.bowtie2.bam SRR9670817.bowtie2.bam 
	samtools merge gr4.bowtie2.bam SRR9670818.bowtie2.bam SRR9670819.bowtie2.bam 

gr4.bed: 
	bedtools bamtobed -i sd3.bowtie2.bam > sd3.bed
	bedtools bamtobed -i gr4.bowtie2.bam > gr4.bed 
	

SRR9670822_SRR9670820_50_500.LR.png: 
	$(foreach i, $(CUTOFFS),python ribofilio.py -t yeast.fa -f SRR9670822.bed -r SRR9670820.bed -b 50 -c ${i};)
	$(foreach i, $(CUTOFFS),python ribofilio.py -t yeast.fa -f SRR9670823.bed -r SRR9670821.bed -b 50 -c ${i};) 

gr3_SRR9670815_50_500.LR.png:
	$(foreach i, $(CUTOFFS),python ribofilio.py -t yeast.fa -f sd3.bed -r SRR9670814.bed -b 50 -c ${i};)
	$(foreach i, $(CUTOFFS),python ribofilio.py -t yeast.fa -f gr4.bed -r SRR9670815.bed -b 50 -c ${i};)


#Run ribofilio on genes clustered by GO
gr4_SRR9670815.GO_0042254_50_0.LR.png:
	$(foreach i, $(GO),python ribofilio.py -t yeast.fa -f SRR9670822.bed -r SRR9670820.bed -b 50 --subset ${i}.txt;)
	$(foreach i, $(GO),python ribofilio.py -t yeast.fa -f SRR9670823.bed -r SRR9670821.bed -b 50 --subset ${i}.txt;)
	$(foreach i, $(GO),python ribofilio.py -t yeast.fa -f sd3.bed -r SRR9670814.bed -b 50 --subset ${i}.txt;)
	$(foreach i, $(GO),python ribofilio.py -t yeast.fa -f gr4.bed -r SRR9670815.bed -b 50 --subset ${i}.txt;)



#Run ribofilio on genes clustered by Chromosomes
gr4_SRR9670815.yeastchrXVI_50_0.LR.png:
	$(foreach i, $(CHRs),python ribofilio.py -t yeast.fa -f SRR9670822.bed -r SRR9670820.bed -b 50 --subset yeast${i}.txt;)
	$(foreach i, $(CHRs),python ribofilio.py -t yeast.fa -f SRR9670823.bed -r SRR9670821.bed -b 50 --subset yeast${i}.txt;)
	$(foreach i, $(CHRs),python ribofilio.py -t yeast.fa -f sd3.bed -r SRR9670814.bed -b 50 --subset yeast${i}.txt;)
	$(foreach i, $(CHRs),python ribofilio.py -t yeast.fa -f gr4.bed -r SRR9670815.bed -b 50 --subset yeast${i}.txt;)


#Run ribofilio on orthogolus genes to human and mouse
gr4_SRR9670815.yeastorthuman_50_0.png: yeastchrXVI.txt
	python ribofilio.py -t yeast.fa -f SRR9670822.bed -r SRR9670820.bed -b 50 --subset yeastorthuman.txt
	python ribofilio.py -t yeast.fa -f SRR9670823.bed -r SRR9670821.bed -b 50 --subset yeastorthuman.txt
	python ribofilio.py -t yeast.fa -f sd3.bed -r SRR9670814.bed -b 50 --subset yeastorthuman.txt
	python ribofilio.py -t yeast.fa -f gr4.bed -r SRR9670815.bed -b 50 --subset yeastorthuman.txt 
	
	python ribofilio.py -t yeast.fa -f SRR9670822.bed -r SRR9670820.bed -b 50 --subset yeastortmouse.txt 
	python ribofilio.py -t yeast.fa -f SRR9670823.bed -r SRR9670821.bed -b 50 --subset yeastortmouse.txt
	python ribofilio.py -t yeast.fa -f sd3.bed -r SRR9670814.bed -b 50 --subset yeastortmouse.txt
	python ribofilio.py -t yeast.fa -f gr4.bed -r SRR9670815.bed -b 50 --subset yeastortmouse.txt

#Get TPM and high and low expressed genes
SRR9670821_l50.tpm:
	salmon quant -i SalmonYeastIndex --libType U -r galore/SRR9670814_trimmed.fq.gz  -o SalmonSRR9670814 \
	-p 3 --writeUnmappedNames --seqBias --gcBias --validateMappings
	@cmd1 = `cut -f1,4 SalmonSRR9670814/quant.sf > SRR9670814.tpm` ; echo $$cmd1
	@cmd2 = `awk '$$2 >50 {print $$0}' SRR9670814.tpm > SRR9670814_g50.tpm` ; echo $$cmd2
	@cmd3 = `awk '$$2 <50 {print $$0}' SRR9670814.tpm > SRR9670814_l50.tpm` ; echo $$cmd3

	salmon quant -i SalmonYeastIndex --libType U -r galore/SRR9670815_trimmed.fq.gz  -o SalmonSRR9670815 \
	-p 3 --writeUnmappedNames --seqBias --gcBias --validateMappings
	@cmd1 = `cut -f1,4 SalmonSRR9670815/quant.sf > SRR9670815.tpm` ; echo $$cmd1
	@cmd2 = `awk '$$2 >50 {print $$0}' SRR9670815.tpm > SRR9670815_g50.tpm` ; echo $$cmd2
	@cmd3 = `awk '$$2 <50 {print $$0}' SRR9670815.tpm > SRR9670815_l50.tpm` ; echo $$cmd3

	salmon quant -i SalmonYeastIndex --libType U -r galore/SRR9670820_trimmed.fq.gz  -o SalmonSRR9670820 \
	-p 3 --writeUnmappedNames --seqBias --gcBias --validateMappings
	@cmd1 = `cut -f1,4 SalmonSRR9670820/quant.sf > SRR9670820.tpm` ; echo $$cmd1
	@cmd2 = `awk '$$2 >50 {print $$0}' SRR9670820.tpm > SRR9670820_g50.tpm` ; echo $$cmd2
	@cmd3 = `awk '$$2 <50 {print $$0}' SRR9670820.tpm > SRR9670820_l50.tpm` ; echo $$cmd3


	salmon quant -i SalmonYeastIndex --libType U -r galore/SRR9670821_trimmed.fq.gz  -o SalmonSRR9670821 \
	-p 3 --writeUnmappedNames --seqBias --gcBias --validateMappings	
	@cmd1 = `cut -f1,4 SalmonSRR9670821/quant.sf > SRR9670821.tpm` ; echo $$cmd1
	@cmd2 = `awk '$$2 >50 {print $$0}' SRR9670821.tpm > SRR9670821_g50.tpm` ; echo $$cmd2	
	@cmd3 = `awk '$$2 <50 {print $$0}' SRR9670821.tpm > SRR9670821_l50.tpm` ; echo $$cmd3

#Run ribofilio on a subset of genes: high and low expressed genes
SRR9670823.SRR9670821.l50tpm.LogLinear.png: 
	python ribofilio.py -t yeast.fa -f sd3.bed -r SRR9670814.bed -s SRR9670814_g50.tpm -b 50
	python ribofilio.py -t yeast.fa -f sd3.bed -r SRR9670814.bed -s SRR9670814_l50.tpm -b 50
	python ribofilio.py -t yeast.fa -f gr4.bed -r SRR9670815.bed -s SRR9670815_g50.tpm -b 50
	python ribofilio.py -t yeast.fa -f gr4.bed -r SRR9670815.bed -s SRR9670815_l50.tpm -b 50
	python ribofilio.py -t yeast.fa -f SRR9670822.bed -r SRR9670820.bed -s SRR9670820_g50.tpm -b 50
	python ribofilio.py -t yeast.fa -f SRR9670822.bed -r SRR9670820.bed -s SRR9670820_l50.tpm -b 50
	python ribofilio.py -t yeast.fa -f SRR9670823.bed -r SRR9670821.bed -s SRR9670821_g50.tpm -b 50
	python ribofilio.py -t yeast.fa -f SRR9670823.bed -r SRR9670821.bed -s SRR9670821_l50.tpm -b 50


#Generate vcf file
SRR9670821.vcf: SRR5945809.bowtie2.bam
	samtools sort -T SRR9670814.tmp.sorted -o SRR9670814.bowtie2.sorted.bam SRR9670814.bowtie2.bam
	samtools index SRR9670814.bowtie2.sorted.bam SRR9670814.bowtie2.sorted.bam.bai
	samtools mpileup -uD -s -S -f yeast.fa SRR9670814.bowtie2.sorted.bam | bcftools  call -v -Ob > SRR9670814.bcf
	bcftools view SRR9670814.bcf > SRR9670814.vcf
	
	samtools sort -T SRR9670815.tmp.sorted -o SRR9670815.bowtie2.sorted.bam SRR9670815.bowtie2.bam
	samtools index SRR9670815.bowtie2.sorted.bam SRR9670815.bowtie2.sorted.bam.bai
	samtools mpileup -uD -s -S -f yeast.fa SRR9670815.bowtie2.sorted.bam | bcftools  call -v -Ob > SRR9670815.bcf
	bcftools view SRR9670815.bcf > SRR9670815.vcf
	
	samtools sort -T SRR9670820.tmp.sorted -o SRR9670820.bowtie2.sorted.bam SRR9670820.bowtie2.bam
	samtools index SRR9670820.bowtie2.sorted.bam SRR9670820.bowtie2.sorted.bam.bai
	samtools mpileup -uD -s -S -f yeast.fa SRR9670820.bowtie2.sorted.bam | bcftools  call  -v -Ob > SRR9670820.bcf
	bcftools view SRR9670820.bcf > SRR9670820.vcf
	
	samtools sort -T SRR9670821.tmp.sorted -o SRR9670821.bowtie2.sorted.bam SRR9670821.bowtie2.bam
	samtools index SRR9670821.bowtie2.sorted.bam SRR9670821.bowtie2.sorted.bam.bai
	samtools mpileup -uD -s -S -f yeast.fa SRR9670821.bowtie2.sorted.bam | bcftools  call  -v -Ob > SRR9670821.bcf
	bcftools view SRR9670821.bcf > SRR9670821.vcf


#Prepare subset file one with genes with zero snps and one with genes with at least one snp
SRR9670821_1.snps: #outputs genes list with zero snps and genes with at least one snp
	python countsnps.py  yeast.fa SRR9670814.vcf --cutoff 1
	python countsnps.py  yeast.fa SRR9670815.vcf --cutoff 1
	python countsnps.py  yeast.fa SRR9670820.vcf --cutoff 1
	python countsnps.py  yeast.fa SRR9670821.vcf --cutoff 1

#Run ribofilio on genes with zero snps and one or more snp 

SRR9670823_SRR9670821.SRR9670821_1_50_0.LR.png: 
	python ribofilio.py -t yeast.fa -f  sd3.bed -r SRR9670814.bed  -b 50 --subset SRR9670814_zerosnps.txt 
	python ribofilio.py -t yeast.fa -f  sd3.bed -r SRR9670814.bed  -b 50 --subset SRR9670814_1snps.txt
	
	python ribofilio.py -t yeast.fa -f  gr4.bed -r SRR9670815.bed  -b 50 --subset SRR9670815_zerosnps.txt 
	python ribofilio.py -t yeast.fa -f  gr4.bed -r SRR9670815.bed  -b 50 --subset SRR9670815_1snps.txt
		
	python ribofilio.py -t yeast.fa -f  SRR9670822.bed -r SRR9670820.bed  -b 50 --subset SRR9670820_zerosnps.txt
	python ribofilio.py -t yeast.fa -f  SRR9670822.bed -r SRR9670820.bed  -b 50 --subset SRR9670820_1snps.txt
	
	python ribofilio.py -t yeast.fa -f  SRR9670823.bed -r SRR9670821.bed  -b 50 --subset SRR9670821_zerosnps.txt 
	python ribofilio.py -t yeast.fa -f  SRR9670823.bed -r SRR9670821.bed  -b 50 --subset SRR9670821_1snps.txt

#--------------------------------------------
#Case Study: GSE102837 
#--------------------------------------------

SRR5945808.fastq.gz:
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR594/009/SRR5945809/SRR5945809.fastq.gz
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR594/008/SRR5945808/SRR5945808.fastq.gz

galore/SRR5945808.trimmed.fastq.gz: 
	trim_galore --gzip --retain_unpaired --trim1 -a "CTGTAGGCACCATCAATAGATCGGAAGAGCACACGTCTGAACTCCAGTCA" --fastqc --fastqc_args "--outdir ${WHERE}/fastqc" -o ${WHERE}/galore SRR5945809.fastq.gz
	trim_galore --gzip --retain_unpaired --trim1 -a "CTGTAGGCACCATCAATAGATCGGAAGAGCACACGTCTGAACTCCAGTCA" --fastqc --fastqc_args "--outdir ${WHERE}/fastqc" -o ${WHERE}/galore SRR5945808.fastq.gz


SRR5945808.bowtie2.bam: 
	bowtie2 -x yeast -U galore/SRR5945809_trimmed.fq.gz  -S SRR5945809.bowtie2.sam
	samtools view -S -b SRR5945809.bowtie2.sam > SRR5945809.bowtie2.bam
	bowtie2 -x yeast -U galore/SRR5945808_trimmed.fq.gz  -S SRR5945808.bowtie2.sam
	samtools view -S -b SRR5945808.bowtie2.sam > SRR5945808.bowtie2.bam


SRR5945808.bed: 
	bedtools bamtobed -i SRR5945809.bowtie2.bam > SRR5945809.bed 
	bedtools bamtobed -i SRR5945808.bowtie2.bam > SRR5945808.bed 

SRR5945809_SRR5945808_50_500.LR.png: SRR5945808.bed
	$(foreach i, $(CUTOFFS),python ribofilio.py -t yeast.fa -f SRR5945809.bed -r SRR5945808.bed -b 50 -c ${i};)


#Get TPM and high and low expressed genes
SRR5945808_l50.tpm: SalmonYeastIndex
	salmon quant -i SalmonYeastIndex --libType U  -r ${WHERE}/galore/SRR5945808_trimmed.fq.gz \
	-o SalmonSRR5945808 -p 3 --writeUnmappedNames --seqBias --gcBias --validateMappings
	@cmd1 = `cut -f1,4 SalmonSRR5945808/quant.sf > SRR5945808.tpm` ; echo $$cmd1
	@cmd2 = `awk '$$2 >50 {print $$0}' SRR5945808.tpm > SRR5945808_g50.tpm` ; echo $$cmd2
	@cmd3 = `awk '$$2 <50 {print $$0}' SRR5945808.tpm > SRR5945808_l50.tpm` ; echo $$cmd3
       

#Run ribofilio on a subset of genes: high and low expressed genes
SRR5945809.SRR5945808.l50tpm.LogLinear.png: SRR5945808_l50.tpm
	python ribofilio.py -t yeast.fa -f SRR5945809.bed -r SRR5945808.bed -s SRR5945808_g50.tpm -b 50
	python ribofilio.py -t yeast.fa -f SRR5945809.bed -r SRR5945808.bed -s SRR5945808_l50.tpm -b 50

#Run ribofilio on genes clustered by GO
SRR5945809_SRR5945808.GO_0042254_50_0.LR.png:yeastchrXVI.txt
	$(foreach i, $(GO),python ribofilio.py -t yeast.fa -f SRR5945809.bed -r SRR5945808.bed -b 50 --subset ${i}.txt;)


#Run ribofilio on genes clustered by Chromosomes 
SRR5945809_SRR5945808.yeastchrXVI_50_0.LR.png:
	$(foreach i, $(CHRs),python ribofilio.py -t yeast.fa -f SRR5945809.bed -r SRR5945808.bed -b 50 --subset yeast${i}.txt;) 


#Run ribofilio on orthogolus genes to human and mouse
SRR5945809_SRR5945808.yeastorthuman_50_0.png: yeastchrXVI.txt
	python ribofilio.py -t yeast.fa -f SRR5945809.bed -r SRR5945808.bed  -b 50 --subset yeastorthuman.txt
	python ribofilio.py -t yeast.fa -f SRR5945809.bed -r SRR5945808.bed  -b 50 --subset yeastortmouse.txt

#Generate vcf file 
SRR5945808.vcf: SRR5945809.bowtie2.bam
	samtools sort -T SRR5945808.tmp.sorted -o SRR5945808.bowtie2.sorted.bam SRR5945808.bowtie2.bam
	samtools index SRR5945808.bowtie2.sorted.bam SRR5945808.bowtie2.sorted.bam.bai
	samtools mpileup -uD -s -S -f yeast.fa SRR5945808.bowtie2.sorted.bam | bcftools  call -m -v -Ob > SRR5945808.bcf
	bcftools view SRR5945808.bcf > SRR5945808.vcf 

#Prepare subset file one with genes with zero snps and one with genes with at least one snp  
SRR5945808_1.snps: #outputs genes list with zero snps and genes with at least one snp  
	python countsnps.py  yeast.fa SRR5945808.vcf --cutoff 1

#Run ribofilio on genes with zero snps and one or more snp 

SRR5945808_SRR5945808.SRR5945808_1_50_0.LR.png: 
	python ribofilio.py -t yeast.fa -f  SRR5945809.bed -r SRR5945808.bed  -b 50 --subset SRR5945808_zerosnps.txt 
	python ribofilio.py -t yeast.fa -f  SRR5945809.bed -r SRR5945808.bed  -b 50 --subset SRR5945808_1snps.txt

#---------------------------------------------
#Case Study:GSE91068
#----------------------------------------------
s31 = SRR5090934
s32 = SRR5090935
s33 = SRR5090936
s34 = SRR5090937

SAMPLES3 = ${s31} ${s32} ${s33} ${s34}

SRR5090937.fastq:
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
SRR5090937_trimmed.fastq.gz: SRR5090937.fastq
	$(foreach i, $(SAMPLES3), trim_galore --gzip --retain_unpaired --trim1 -a "CTGTAGGCACCATCAAT" --fastqc --fastqc_args "--outdir ${WHERE}/fastqc" -o ${WHERE}/galore  ${WHERE}/$(i).fastq ; )

SRR5090937.bowtie2.bam: SRR5090937_trimmed.fastq.gz
	$(foreach i, $(SAMPLES3), bowtie2 -x yeast -U galore/$(i)_trimmed.fq.gz  -S $(i).bowtie2.sam;)
	$(foreach i, $(SAMPLES3), samtools view -S -b $(i).bowtie2.sam > $(i).bowtie2.bam;)

SRR5090937.bed: SRR5090937.bowtie2.bam
	$(foreach i, $(SAMPLES3), bedtools bamtobed -i $(i).bowtie2.bam > $(i).bed;)

SRR5090936_SRR5090934.200.LR.png: #SRR5090937.bed 
	$(foreach i, $(CUTOFFS), python ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed -b 50 -c ${i};)
	$(foreach i, $(CUTOFFS), python ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -b 50 -c ${i};)

#Get TPM and high and low expressed genes  
SRR5090934_l50.tpm: SalmonYeastIndex
	salmon quant -i SalmonYeastIndex --libType U  -r ${WHERE}/galore/SRR5090935_trimmed.fq.gz \
	-o SalmonSRR5090935 -p 3 --writeUnmappedNames --seqBias --gcBias --validateMappings
	@cmd1 = `cut -f1,4 SalmonSRR5090935/quant.sf > SRR5090935.tpm` ; echo $$cmd1
	@cmd2 = `awk '$$2 >50 {print $$0}' SRR5090935.tpm > SRR5090935_g50.tpm` ; echo $$cmd2
	@cmd3 = `awk '$$2 <50 {print $$0}' SRR5090935.tpm > SRR5090935_l50.tpm` ; echo $$cmd3
	salmon quant -i SalmonYeastIndex --libType U  -r ${WHERE}/galore/SRR5090934_trimmed.fq.gz \
	-o SalmonSRR5090934 -p 3 --writeUnmappedNames --seqBias --gcBias --validateMappings
	@cmd1 = `cut -f1,4 SalmonSRR5090934/quant.sf > SRR5090934.tpm` ; echo $$cmd1
	@cmd2 = `awk '$$2 >50 {print $$0}' SRR5090934.tpm > SRR5090934_g50.tpm` ; echo $$cmd2
	@cmd3 = `awk '$$2 <50 {print $$0}' SRR5090934.tpm > SRR5090934_l50.tpm` ; echo $$cmd3

#Run ribofilio on a subset of genes: high and low expressed genes 
SRR5090936.SRR5090934.l50tpm.LogLinear.png: #SRR5090934_l50.tpm 
	python ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed -s SRR5090935_g50.tpm -b 50
	python ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed -s SRR5090935_l50.tpm -b 50
	python ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -s SRR5090934_g50.tpm -b 50
	python ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -s SRR5090934_l50.tpm -b 50

#Run ribofilio on genes clustered by GO  
SRR5090936.SRR5090934.GO.LR.png: yeastchrXVI.txt
	$(foreach i, $(GO),python ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed -b 50 -s ${i}.txt ;)
	$(foreach i, $(GO),python ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -b 50 -s ${i}.txt ;)


#Run ribofilio on genes clustered by Chromosomes
SRR5090936_SRR5090934.yeastchrXVI_50_0.LR.png:
	$(foreach i, $(CHRs),python ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed -b 50 --subset yeast${i}.txt;)
	$(foreach i, $(CHRs),python ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -b 50 --subset yeast${i}.txt;)

#Run ribofilio on orthogolus genes to human and mouse 
SRR5090936_SRR5090934.yeastorthuman_50_0.png: #yeastchrXVI.txt 
	python ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed  -b 50 --subset yeastorthuman.txt
	python ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed  -b 50 --subset yeastortmouse.txt
	python ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed  -b 50 --subset yeastorthuman.txt
	python ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed  -b 50 --subset yeastortmouse.txt

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
SRR5090934_1.snps: #outputs genes list with zero snps and genes with at least one snp
	python countsnps.py  yeast.fa SRR5090935.vcf --cutoff 1
	python countsnps.py  yeast.fa SRR5090934.vcf --cutoff 1

#Run ribofilio on genes with zero snps and one or more snp

SRR5090936_SRR5090934.SRR5090934_1_50_0.LR.png:
	python ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed  -b 50 --subset SRR5090935_zerosnps.txt 
	python ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed  -b 50 --subset SRR5090935_1snps.txt
	python ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed  -b 50 --subset SRR5090934_zerosnps.txt
	python ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed  -b 50 --subset SRR5090934_1snps.txt


#-----------------------------------------
#GSE13750
#-----------------------------------------

s35 = SRR014368
s36 = SRR014369
s37 = SRR014370
s38 = SRR014371
s39 = SRR014372
s40 = SRR014373
s41 = SRR014374
s42 = SRR014375
s43 = SRR014376
s44 = SRR014377
s45 = SRR014378
s46 = SRR014379
s47 = SRR014380
s48 = SRR014381
s49 = SRR014382
s50 = SRR014383
s51 = SRR014384
s52 = SRR014385
s53 = SRR014386
s54 = SRR014387
s55 = SRR028774	

SRR014387.fastq: SRR9670816.fastq
	mv SRR014368.3.fastq SRR014368.fastq
	mv SRR014369.3.fastq SRR014369.fastq
	mv SRR014370.2.fastq SRR014370.fastq
	mv SRR014371.3.fastq SRR014371.fastq
	mv SRR014372.2.fastq SRR014372.fastq
	mv SRR014373.2.fastq SRR014373.fastq
	mv SRR014374.3.fastq SRR014374.fastq
	mv SRR014375.3.fastq SRR014375.fastq
	mv SRR014376.2.fastq SRR014376.fastq
	mv SRR014377.2.fastq SRR014377.fastq
	mv SRR014378.2.fastq SRR014378.fastq
	mv SRR014379.2.fastq SRR014379.fastq
	mv SRR014380.2.fastq SRR014380.fastq
	mv SRR014381.2.fastq SRR014381.fastq
	mv SRR014382.3.fastq SRR014382.fastq
	mv SRR014383.2.fastq SRR014383.fastq
	mv SRR014384.2.fastq SRR014384.fastq
	mv SRR014385.3.fastq SRR014385.fastq
	mv SRR014386.2.fastq SRR014386.fastq
	mv SRR014387.3.fastq SRR014387.fastq 
	mv SRR028774.3.fastq SRR028774.fastq 

SAMPLES4 = ${s35} ${s36} ${s37}	${s38} ${s39} ${s40} ${s41} ${s42} ${s43} ${s44} ${s45} ${s46} ${s47} ${s48} ${s49} ${s50} ${s51} ${s52} ${s53} ${s54} ${s55}

SRR028774_trimmed.fastq.gz:
	$(foreach i, $(SAMPLES4), trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir ${WHERE}/fastqc" -o ${WHERE}/galore  ${WHERE}/$(i).fastq ; )

SRR028774.bowtie2.bam:
	$(foreach i, $(SAMPLES4), bowtie2 -x yeast -U galore/$(i)_trimmed.fq.gz  -S $(i).bowtie2.sam;)
	$(foreach i, $(SAMPLES4), samtools view -S -b $(i).bowtie2.sam > $(i).bowtie2.bam;)

#https://www.ncbi.nlm.nih.gov/sra?term=SRP000637
fp-starved-1.bam:
	samtools merge mRNA-rich-2.bam SRR014386.bowtie2.bam SRR014387.bowtie2.bam SRR028774.bowtie2.bam
	mv SRR014385.bowtie2.bam  mRNA-rich-1.bam 
	samtools merge mRNA-starved-2.bam SRR014383.bowtie2.bam SRR014384.bowtie2.bam 	
	mv SRR014382.bowtie2.bam mRNA-starved-1.bam 
	samtools merge fp-rich-2.bam SRR014377.bowtie2.bam SRR014378.bowtie2.bam SRR014379.bowtie2.bam SRR014380.bowtie2.bam SRR014381.bowtie2.bam
	samtools merge fp-rich-1.bam SRR014374.bowtie2.bam SRR014375.bowtie2.bam SRR014376.bowtie2.bam
	samtools merge fp-starved-2.bam SRR014370.bowtie2.bam SRR014371.bowtie2.bam SRR014372.bowtie2.bam SRR014373.bowtie2.bam
	samtools merge fp-starved-1.bam SRR014368.bowtie2.bam SRR014369.bowtie2.bam 

fp-starved-1.bed:
	bedtools bamtobed -i mRNA-rich-2.bam > mRNA-rich-2.bed 
	bedtools bamtobed -i mRNA-rich-1.bam > mRNA-rich-1.bed 
	bedtools bamtobed -i mRNA-starved-2.bam > mRNA-starved-2.bed 
	bedtools bamtobed -i mRNA-starved-1.bam >  mRNA-starved-1.bed 
	bedtools bamtobed -i fp-rich-2.bam > fp-rich-2.bed 
	bedtools bamtobed -i fp-rich-1.bam > fp-rich-1.bed 
	bedtools bamtobed -i fp-starved-2.bam > fp-starved-2.bed 
	bedtools bamtobed -i fp-starved-1.bam > fp-starved-1.bed 

fp-starved-1_mRNA-starved-1.50.2000.*.LR.png:
	$(foreach i, $(CUTOFFS), python ribofilio.py -t yeast.fa -f fp-rich-1.bed -r mRNA-rich-1.bed -b 50 -c ${i};) 
	$(foreach i, $(CUTOFFS),python ribofilio.py -t yeast.fa -f fp-starved-1.bed -r mRNA-starved-1.bed -b 50 -c ${i};) 


fp-starved-2_mRNA-starved-2.50.2000.*.LR.png:
	$(foreach i, $(CUTOFFS), python ribofilio.py -t yeast.fa -f fp-rich-2.bed -r mRNA-rich-2.bed -b 50 -c ${i};)
	$(foreach i, $(CUTOFFS), python ribofilio.py -t yeast.fa -f fp-starved-2.bed -r mRNA-starved-2.bed -b 50 -c ${i};)


#Get TPM and high and low expressed genes
mRNA-starved-1_l50.tpm: SalmonYeastIndex
	salmon quant -i SalmonYeastIndex --libType U -r galore/SRR014385_trimmed.fq.gz -o SalmonmRNA-rich-1 \
	-p 3 --writeUnmappedNames --seqBias --gcBias --validateMappings
	@cmd1 = `cut -f1,4 SalmonmRNA-rich-1/quant.sf > mRNA-rich-1.tpm` ; echo $$cmd1
	@cmd2 = `awk '$$2 >50 {print $$0}' mRNA-rich-1.tpm > mRNA-rich-1_g50.tpm` ; echo $$cmd2
	@cmd3 = `awk '$$2 <50 {print $$0}' mRNA-rich-1.tpm > mRNA-rich-1_l50.tpm` ; echo $$cmd3
	
	salmon quant -i SalmonYeastIndex --libType U -r galore/SRR014382_trimmed.fq.gz -o SalmonmRNA-starved-1 \
	-p 3 --writeUnmappedNames --seqBias --gcBias --validateMappings
	@cmd1 = `cut -f1,4 SalmonmRNA-starved-1/quant.sf > mRNA-starved-1.tpm` ; echo $$cmd1
	@cmd2 = `awk '$$2 >50 {print $$0}' mRNA-starved-1.tpm > mRNA-starved-1_g50.tpm` ; echo $$cmd2
	@cmd3 = `awk '$$2 <50 {print $$0}' mRNA-starved-1.tpm > mRNA-starved-1_l50.tpm` ; echo $$cmd3

mRNA-starved-2_l50.tpm:
	#cat galore/SRR014386_trimmed.fq.gz  galore/SRR014387_trimmed.fq.gz  galore/SRR028774_trimmed.fq.gz > mRNA-rich-2.fastq.gz 
	salmon quant -i SalmonYeastIndex --libType U -r mRNA-rich-2.fastq.gz -o SalmonmRNA-rich-2 \
	-p 3 --writeUnmappedNames --seqBias --gcBias --validateMappings 
	@cmd1 = `cut -f1,4 SalmonmRNA-rich-2/quant.sf > mRNA-rich-2.tpm` ; echo $$cmd1
	@cmd2 = `awk '$$2 >50 {print $$0}' mRNA-rich-2.tpm > mRNA-rich-2_g50.tpm` ; echo $$cmd2
	@cmd3 = `awk '$$2 <50 {print $$0}' mRNA-rich-2.tpm > mRNA-rich-2_l50.tpm` ; echo $$cmd3
	#cat galore/SRR014383_trimmed.fq.gz galore/SRR014384 _trimmed.fq.gz > mRNA-starved-2.fastq.gz  
	
	salmon quant -i SalmonYeastIndex --libType U -r mRNA-starved-2.fastq.gz -o SalmonmRNA-starved-2 \
	-p 3 --writeUnmappedNames --seqBias --gcBias --validateMappings 
	@cmd1 = `cut -f1,4 SalmonmRNA-starved-2/quant.sf > mRNA-starved-2.tpm` ; echo $$cmd1
	@cmd2 = `awk '$$2 >50 {print $$0}' mRNA-starved-2.tpm > mRNA-starved-2_g50.tpm` ; echo $$cmd2
	@cmd3 = `awk '$$2 <50 {print $$0}' mRNA-starved-2.tpm > mRNA-starved-2_l50.tpm` ; echo $$cmd3

#Run ribofilio on a subset of genes: high and low expressed genes
fp-starved-1.mRNA-starved-1.l50tpm.LogLinear.png: #mRNA-starved-1_l50.tpm
	python ribofilio.py -t yeast.fa -f fp-rich-1.bed -r mRNA-rich-1.bed -s mRNA-rich-1_g50.tpm -b 50
	python ribofilio.py -t yeast.fa -f fp-rich-1.bed -r mRNA-rich-1.bed -s mRNA-rich-1_l50.tpm -b 50
	python ribofilio.py -t yeast.fa -f fp-starved-1.bed -r mRNA-starved-1.bed -s mRNA-starved-1_g50.tpm -b 50
	python ribofilio.py -t yeast.fa -f fp-starved-1.bed -r mRNA-starved-1.bed -s mRNA-starved-1_l50.tpm -b 50




#Run ribofilio on genes clustered by GO
fp-starved-1_mRNA-starved-1.GO_0006950.LR.png: #yeastchrXVI.txt
	$(foreach i, $(GO),python ribofilio.py -t yeast.fa -f fp-rich-1.bed -r mRNA-rich-1.bed -b 50 -s ${i}.txt ;)
	$(foreach i, $(GO),python ribofilio.py -t yeast.fa -f fp-starved-1.bed -r mRNA-starved-1.bed -b 50 -s ${i}.txt ;)


#Run ribofilio on genes clustered by Chromosomes
fp-rich-1_mRNA-rich-1.yeastchrXVI_50_0.LR.png:
	$(foreach i, $(CHRs),python ribofilio.py -t yeast.fa -f fp-rich-1.bed -r mRNA-rich-1.bed -b 50 --subset yeast${i}.txt;)
	$(foreach i, $(CHRs),python ribofilio.py -t yeast.fa -f fp-starved-1.bed -r mRNA-starved-1.bed -b 50 --subset yeast${i}.txt;)


#Run ribofilio on orthogolus genes to human and mouse
fp-rich-1_mRNA-rich-1.yeastorthuman_50_0.png: #yeastchrXVI.txt
	python ribofilio.py -t yeast.fa -f fp-rich-1.bed -r mRNA-rich-1.bed  -b 50 --subset yeastorthuman.txt
	python ribofilio.py -t yeast.fa -f fp-rich-1.bed -r mRNA-rich-1.bed  -b 50 --subset yeastortmouse.txt
	python ribofilio.py -t yeast.fa -f fp-starved-1.bed -r mRNA-starved-1.bed  -b 50 --subset yeastorthuman.txt
	python ribofilio.py -t yeast.fa -f fp-starved-1.bed -r mRNA-starved-1.bed  -b 50 --subset yeastortmouse.txt

#Generate vcf file 
mRNA-starved-1.vcf: mRNA-starved-1.bowtie2.bam
	samtools sort -T mRNA-rich-1.tmp.sorted -o mRNA-rich-1.bowtie2.sorted.bam mRNA-rich-1.bam
	samtools index mRNA-rich-1.bowtie2.sorted.bam mRNA-rich-1.bowtie2.sorted.bam.bai
	samtools mpileup -uD -s -S -f yeast.fa mRNA-rich-1.bowtie2.sorted.bam | bcftools  call -m -v -Ob > mRNA-rich-1.bcf
	bcftools view mRNA-rich-1.bcf > mRNA-rich-1.vcf 
	samtools sort -T mRNA-starved-1.tmp.sorted -o mRNA-starved-1.bowtie2.sorted.bam mRNA-starved-1.bam
	samtools index mRNA-starved-1.bowtie2.sorted.bam mRNA-starved-1.bowtie2.sorted.bam.bai
	samtools mpileup -uD -s -S -f yeast.fa mRNA-starved-1.bowtie2.sorted.bam | bcftools  call -m -v -Ob > mRNA-starved-1.bcf
	bcftools view mRNA-starved-1.bcf > mRNA-starved-1.vcf 

#Prepare subset file one with genes with zero snps and one with genes with at least one snp
mRNA-starved-1.snps: #outputs genes list with zero snps and genes with at least one snp
	python countsnps.py  yeast.fa mRNA-rich-1.vcf --cutoff 1
	python countsnps.py  yeast.fa mRNA-starved-1.vcf --cutoff 1

#Run ribofilio on genes with zero snps and one or more snp

fp-starved-1_mRNA-starved-1.mRNA-starved-1_1_50_0.LR.png:
	python ribofilio.py -t yeast.fa -f fp-rich-1.bed -r mRNA-rich-1.bed  -b 50 --subset mRNA-rich-1_zerosnps.txt 
	python ribofilio.py -t yeast.fa -f fp-rich-1.bed -r mRNA-rich-1.bed  -b 50 --subset mRNA-rich-1_1snps.txt
	python ribofilio.py -t yeast.fa -f fp-starved-1.bed -r mRNA-starved-1.bed  -b 50 --subset mRNA-starved-1_zerosnps.txt
	python ribofilio.py -t yeast.fa -f fp-starved-1.bed -r mRNA-starved-1.bed  -b 50 --subset mRNA-starved-1_1snps.txt

