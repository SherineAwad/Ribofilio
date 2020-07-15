BBMAP = /home/sa941/software/bbmap

YEAST = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/annotations/Saccharomyces_cerevisiae/Ensembl/R64-1-1


WHERE  = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/RibosomalProfiling
#https://www.ebi.ac.uk/ena/data/view/SRR8360209 


#----------------------------------------------------------------------------------
#Download Samples
#----------------------------------------------------------------------------------
SRR8360209.fastq.gz:
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR836/009/SRR8360209/SRR8360209.fastq.gz
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR836/004/SRR8360214/SRR8360214.fastq.gz
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR721/006/SRR7215536/SRR7215536.fastq.gz
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR721/007/SRR7215537/SRR7215537.fastq.gz
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR721/002/SRR7215542/SRR7215542.fastq.gz
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR721/003/SRR7215543/SRR7215543.fastq.gz
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR645/007/SRR6452757/SRR6452757.fastq.gz
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR645/008/SRR6452758/SRR6452758.fastq.gz
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR726/003/SRR7265163/SRR7265163.fastq.gz
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR726/004/SRR7265164/SRR7265164.fastq.gz
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR726/009/SRR7265169/SRR7265169.fastq.gz
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR726/000/SRR7265170/SRR7265170.fastq.gz
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR726/001/SRR7265151/SRR7265151.fastq.gz
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR726/002/SRR7265152/SRR7265152.fastq.gz
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR726/007/SRR7265157/SRR7265157.fastq.gz
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR726/008/SRR7265158/SRR7265158.fastq.gz
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR594/009/SRR5945809/SRR5945809.fastq.gz
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR594/008/SRR5945808/SRR5945808.fastq.gz


s1 = SRR8360209
s2 = SRR8360214
s3 = SRR7215536
s4 = SRR7215537
s5 = SRR7215542
s6 = SRR7215543
s7 = SRR6452757
s8 = SRR6452758
s9 = SRR7265163
s10 = SRR7265164
s11 = SRR7265169
s12 = SRR7265170
s13 = SRR7265151
s14 = SRR7265152
s15 = SRR7265157
s16 = SRR7265158
s17 = SRR5945808 
s18 = SRR5945809
s19 = WT_FP
s20 = WT_RNA
 
SAMPLES = ${s1} ${s2} ${s3} ${s4} ${s5} ${s6} ${s7} ${s8} ${s9} ${s10} ${s11} ${s12} ${s13} ${s14} ${s15} ${s16} ${s17} ${s18} ${s19} ${s20}  
GO = GO_0006119 GO_0006406 GO_0006412 GO_0006950 GO_0007049 GO_0009651 GO_0016458 GO_0023052 GO_0031047 GO_0031990 GO_0042254
CHRs = chrI chrII chrIII chrIV chrV chrVI chrVII chrVIII chrVIII chrIX chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI 

#SAMPLES = "$(shell ls *fastq.gz)"; echo $$SAMPLES ; 

#SAMPLES = $(wildcard *.fastq.gz)

#--------------------------------------------------------------------------
#Get yeats Trabscripts and prepare index and bed files 
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


galore/WT_RNA.trimmed.fastq.gz: 
	$(foreach i, $(SAMPLES), trim_galore --gzip --retain_unpaired --trim1 -a "CTGTAGGCACCATCAATAGATCGGAAGAGCACACGTCTGAACTCCAGTCA" --fastqc --fastqc_args "--outdir ${WHERE}/fastqc" -o ${WHERE}/galore  ${WHERE}/$(i) ; ) 


#-------------------------------------------------------------
#Aligning to transcripts 
#-------------------------------------------------------------
SRR5945809.bowtie2.bam: 
	$(foreach i, $(SAMPLES), bowtie2 -x yeast -U galore/$(i)_trimmed.fq.gz  -S $(i).bowtie2.sam;)
	$(foreach i, $(SAMPLES), samtools view -S -b $(i).bowtie2.sam > $(i).bowtie2.bam;)


SRR5945809.bed: 
	$(foreach i, $(SAMPLES), bedtools bamtobed -i $(i).bowtie2.bam > $(i).bed;) 

#------------------------------------------------
#Run ribosomes for FP/RNA
#------------------------------------------------
WT_FP_WT_RNA.50.profile:         
	python ribosomes.py yeast.fa WT_FP.bed WT_RNA.bed --binsize 50 --mode 0

#------------------------------------------------
#Get SNPs: Quick variant call 
#------------------------------------------------

SRR5945809.bowtie2.sorted.bam:
	$(foreach i, $(XS), samtools sort -T $(i).tmp.sorted -o $(i).bowtie2.sorted.bam $(i).bowtie2.bam;)
	$(foreach i, $(XS), samtools index $(i).bowtie2.sorted.bam $(i).bowtie2.sorted.bam.bai;)


SRR5945809.vcf:
	$(foreach i, $(XS), samtools mpileup -uD -s -S -f yeast.fa $(i).bowtie2.sorted.bam | bcftools  call -m -v -Ob > $(i).bcf;)
	$(foreach i, $(XS), bcftools view $(i).bcf > $(i).vcf;)

WT_RNA_1snps.txt: 
	$(foreach i, $(SAMPLES), python countsnps.py yeast.fa $(i).vcf --cutoff 1;)  #Will generate sample_zero.snps  sample_1.snps for zeros snps and one or more snps respectively  
WT_FP_WT_RNA.50.0.WT_FP_zero.LogLinear.png :
	python ribosomes.py yeast.fa WT_FP.bed WT_RNA.bed --subset WT_FP_1snps.txt --binsize 50 --mode 0 
	python ribosomes.py yeast.fa WT_FP.bed WT_RNA.bed --subset WT_FP_zerosnps.txt --binsize 50 --mode 0 

#------------------------------------------------
#Prepare TPM files
#------------------------------------------------

SalmonYeastIndex: 
	salmon index --index SalmonYeastIndex --transcripts yeast.fa  

SalmonSRR5945809: 
	$(foreach i, $(SAMPLES), salmon quant -i SalmonYeastIndex --libType U  -r ${WHERE}/galore/$(i)_trimmed.fq.gz -o Salmon$(i) -p 3 --writeUnmappedNames --seqBias --gcBias --validateMappings;) 

wtfp_l50.tpm:
	salmon quant -i SalmonYeastIndex --libType U  -r ${WHERE}/galore/WT_FP_trimmed.fq.gz \
	-o SalmonWT_FP -p 3 --writeUnmappedNames --seqBias --gcBias --validateMappings
	@cmd1 = `cut -f1,4 SalmonWT_FP/quant.sf > wtfp.tpm` ; echo $$cmd1  
	@cmd2 = `awk '$$2 >50 {print $$0}' wtfp.tpm > wtfp_g50.tpm` ; echo $$cmd2 
	@cmd3 = `awk '$$2 <50 {print $$0}' wtfp.tpm > wtfp_l50.tpm` ; echo $$cmd3

wtrna_l50.tpm: 
	salmon quant -i SalmonYeastIndex --libType U  -r ${WHERE}/galore/WT_RNA_trimmed.fq.gz \
	-o SalmonWT_RNA -p 3 --writeUnmappedNames --seqBias --gcBias --validateMappings
	@cmd1 = `cut -f1,4 SalmonWT_RNA/quant.sf > wtrna.tpm` ; echo $$cmd1  
	@cmd2 = `awk '$$2 >50 {print $$0}' wtrna.tpm > wtrna_g50.tpm` ; echo $$cmd2 
	@cmd3 = `awk '$$2 >50 {print $$0}' wtrna.tpm > wtrna_l50.tpm` ; echo $$cmd3 


WT_FP_WT_RNA.50.0.wtfpl50tpm.LogLinear.png:
	 python ribosomes.py yeast.fa WT_FP.bed WT_RNA.bed --subset wtfp_g50.tpm --binsize 50 --mode 0 
	 python ribosomes.py yeast.fa WT_FP.bed WT_RNA.bed --subset wtfp_l50.tpm --binsize 50 --mode 0

#------------------------------------------------------------------------------------------------
#Get GO, GC, orthogolous, and Chromosomes from Ensembl Biomart
#------------------------------------------------------------------------------------------------
yeastchrXVI.txt:
	Rscript getBiomart.R  #This script will download the GO, Chromosomes, GC, and orthogolus data 

#------------------------------------------------
#Run Ribosomes for GO
#------------------------------------------------
WT_FP_WT_RNA.50.0.GO_0042254.regression.png:
	$(foreach i, $(GO), python ribosomes.py yeast.fa WT_FP.bed WT_RNA.bed --subset $(i).txt --binsize 50 --mode 0;)

#------------------------------------------------
#Run Ribosomes for Chromosomes
#------------------------------------------------
WT_FP_WT_RNA.50.0.yeastchrXV.LogLinear.png:
	$(foreach i, $(CHRs), python ribosomes.py yeast.fa WT_FP.bed WT_RNA.bed --subset yeast$(i).txt --binsize 50 --mode 0;) 

#-------------------------------------------
#Run Ribosomes for Orthogolus
#-------------------------------------------
WT_FP_WT_RNA.50.0.yeastortmouse.LogLinear.png:
	python ribosomes.py yeast.fa WT_FP.bed WT_RNA.bed --subset yeastorthuman.txt --binsize 50 --mode 0 
	python ribosomes.py yeast.fa WT_FP.bed WT_RNA.bed --subset yeastortmouse.txt --binsize 50 --mode 0 
#------------------------------------------------
#Prepare GC files
#------------------------------------------------

high.GC.txt:
	@cmd1 = `awk '$2 <40 {print $0}' yeastGC.tsv > lowGC.txt` ; echo $$cmd1 	
	@cmd2 = `cut  -f1 lowGC.txt > highGC.genes ` ; echo $$cmd2 
	@cmd3 = `awk '$2 >60 {print $0}' yeastGC.tsv > highGC.txt ` ; echo $$cmd3 
	@cmd4 = `cut  -f1 highGC.txt >highGC.genes ` ; echo $$cmd4 

lowGC.png:
	python ribosomes.py yeast.fa WT_FP.bed WT_RNA.bed --subset highGC.genes --binsize 50 --mode 0
	python ribosomes.py yeast.fa WT_FP.bed WT_RNA.bed --subset lowGC.genes --binsize 50 --mode 0 

#-------------------------------------------------
#Try out SV SUID
#-------------------------------------------------

S = WT_FP

yeast.fa.bwt:
	bwa index yeast.fa 

WT_FP.bwa.sorted.bam:
	#$(foreach i ,$(S), bwa aln yeast.fa  galore/$(i)_trimmed.fq.gz >  $(i).bwa.sai;)  
	$(foreach i ,$(S), bwa samse  yeast.fa  $(i).bwa.sai  galore/$(i)_trimmed.fq.gz  >  $(i).bwa.sam;)
	$(foreach i ,$(S), samtools import  yeast.fasta.fai  $(i).bwa.sam  $(i).bwa.bam;)
	$(foreach i, $(S), samtools sort -T $(i).bwa.tmp.sorted -o $(i).bwa.sorted.bam $(i).bwa.bam;) 
	$(foreach i, $(S), samtools index $(i).bwa.sorted.bam $(i).bwa.sorted.bam.bai;)

WT_FP.bwamem.sorted.bam:
	#bwa mem yeast.fa galore/WT_FP_trimmed.fq.gz > WT_FP.bwamem.sam
	#$(foreach i ,$(S), samtools import  yeast.fasta.fai  $(i).bwamem.sam  $(i).bwamem.bam;)
	$(foreach i, $(S), samtools sort -T $(i).bwamem.tmp.sorted -o $(i).bwamem.sorted.bam $(i).bwamem.bam;) 
	$(foreach i, $(S), samtools index $(i).bwamem.sorted.bam $(i).bwamem.sorted.bam.bai;)

SV:
	 squid --bwa -b WT_FP.bwa.sorted.bam -o wtfpsv #failed for segmentation fault, trying STAR 


yeastStarIndex:       
	STAR --runThreadN 4 --runMode genomeGenerate --genomeDir YeastStarIndex \
--genomeFastaFiles yeast.fa

WT_FP.star_Aligned.sortedByCoord.out.bam:
	$(foreach i, $(S), STAR --genomeDir YeastStarIndex --readFilesIn galore/$(i)_trimmed.fq.gz --readFilesCommand zcat --runThreadN 2 --runMode alignReads --chimOutType SeparateSAMold  --chimSegmentMin 15 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $(i).star_;) 
SV2:
	squid -b WT_FP.star_Aligned.sortedByCoord.out.bam -o wtfpsv 
