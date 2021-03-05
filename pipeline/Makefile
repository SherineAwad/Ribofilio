BBMAP = /home/sa941/software/bbmap

YEAST = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/annotations/Saccharomyces_cerevisiae/Ensembl/R64-1-1

DATASETS = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/RibosomalProfiling/ribofilio/datasets
SRC = /home/ec2-user/ribofilio/pipeline/src

CURRENT  =  /rds/project/yhbl2/rds-yhbl2-genehunter/SM/RibosomalProfiling/ribofilio/pipeline

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

#-----------------------
#Gene Length Experiment
#-----------------------
yeast_6000_100000.txt: yeast.fa 
	python ${SRC}/gcluster.py yeast.fa yeast_0_2000.txt --max 2000 --min 0
	python ${SRC}/gcluster.py yeast.fa yeast_2000_4000.txt --max 4000 --min 2000
	python ${SRC}/gcluster.py yeast.fa yeast_4000_6000.txt --max 6000 --min 4000
	python ${SRC}/gcluster.py yeast.fa yeast_6000_100000.txt --max 100000 --min 6000


#------------------------------------------------------------------------------------------
#Dataset: GSE91068
#------------------------------------------------------------------------------------------
s1 = SRR5090934 #GSM2420486
s2 = SRR5090935 #GSM2420487
s3 = SRR5090936 #GSM2420488
s4 = SRR5090937 #GSM2420489

GSE134152 = ${s1} ${s2} ${s3} ${s4}


#Adpaters reported by dataset author @jiashun@genome.ucsf.edu 
galore/SRR5090937_trimmed.fastq.gz: ${DATASETS}/GSE91068/SRR5090937.fastq
	$(foreach i, $(GSE134152), trim_galore --gzip --retain_unpaired --trim1 -a "CTGTAGGCACCATCAAT" --fastqc --fastqc_args "--outdir ${CURRENT}/fastqc" -o ${CURRENT}/galore  ${DATASETS}/GSE91068/$(i).fastq ; )

SRR5090937.bowtie2.bam: galore/SRR5090937_trimmed.fastq.gz
	$(foreach i, $(GSE134152), bowtie2 -x yeast -U galore/$(i)_trimmed.fq.gz  -S $(i).bowtie2.sam;)
	$(foreach i, $(GSE134152), samtools view -S -b $(i).bowtie2.sam > $(i).bowtie2.bam;)

SRR5090937.bed: 
	$(foreach i, $(GSE134152), bedtools bamtobed -i $(i).bowtie2.bam > $(i).bed;)

#--------------
#SHUFFLE
#--------------

SRR5090937_1000.bed:     
	$(foreach i, $(subset1), python ${SRC}/shuffle.py -t yeast.fa -b SRR5090937.bed  -o SRR5090937.$(i).bed;) 
	$(foreach i, $(subset2), python ${SRC}/shuffle.py -t yeast.fa -b SRR5090937.bed  -o SRR5090937.$(i).bed;) 
SRR5090936_1000.bed: 
	$(foreach i, $(subset1), python ${SRC}/shuffle.py -t yeast.fa -b SRR5090936.bed  -o SRR5090936.$(i).bed;) 
	$(foreach i, $(subset2), python ${SRC}/shuffle.py -t yeast.fa -b SRR5090936.bed  -o SRR5090936.$(i).bed;) 

#----------------------------------
#Main: SRR5090936 and SRR5090937
#----------------------------------
SRR5090936_regfile.log:	
	python  ${SRC}/ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -o SRR5090936
	$(foreach i, $(subset1), python  ${SRC}/ribofilio.py -t yeast.fa -f SRR5090936_$(i).bed -r SRR5090934.bed -o SRR5090936 -p 0;)
	$(foreach i, $(subset2), python  ${SRC}/ribofilio.py -t yeast.fa -f SRR5090936_$(i).bed -r SRR5090934.bed -o SRR5090936 -p 0;)

SRR5090937_regfile.log: 
	python  ${SRC}/ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed -o SRR5090937
	$(foreach i, $(subset1), python  ${SRC}/ribofilio.py -t yeast.fa -f SRR5090937_$(i).bed -r SRR5090935.bed -o SRR5090937 -p 0;)
	$(foreach i, $(subset2), python  ${SRC}/ribofilio.py -t yeast.fa -f SRR5090937_$(i).bed -r SRR5090935.bed -o SRR5090937 -p 0;)


#-------------------------------------
#Gene Length
#-------------------------------------
SRR5090936_yeast_0_2000_regfile.log:
	python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -o SRR5090936_yeast_4000_6000 -s yeast_4000_6000.txt  
	python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -o SRR5090936_yeast_6000_100000 -s yeast_6000_100000.txt 
	python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -o SRR5090936_yeast_2000_4000 -s yeast_2000_4000.txt  
	python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -o SRR5090936_yeast_0_2000 -s yeast_0_2000.txt

#---------------
#Gene Ontology
#---------------

#GO0031990
#----------
SRR5090936_GO0031990_regfile.log:
	python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed  -o SRR5090936_GO0031990 -s ${DATASETS}/GO_0031990.txt 
	$(foreach i, $(subset1), python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090936_$(i).bed -r SRR5090934.bed -o SRR5090936_GO0031990_$(i) -s ${DATASETS}/GO_0031990.txt -p 0 ;)
	$(foreach i, $(subset2), python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090936_$(i).bed -r SRR5090934.bed -o SRR5090936_GO0031990_$(i) -s ${DATASETS}/GO_0031990.txt -p 0 ;)

#GO0009651
#-----------
SRR5090936_GO0009651_regfile.log:
	python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed  -o SRR5090936_GO0009651 -s ${DATASETS}/GO_0009651.txt
	$(foreach i, $(subset1), python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090936_$(i).bed -r SRR5090934.bed  -o SRR5090936_GO0009651_$(i) -s ${DATASETS}/GO_0009651.txt -p 0;)
	$(foreach i, $(subset2), python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090936_$(i).bed -r SRR5090934.bed  -o SRR5090936_GO0009651_$(i) -s ${DATASETS}/GO_0009651.txt -p 0;)

#GO0006950
#-----------
SRR5090936_GO0006950_regfile.log:
	python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed  -o SRR5090936_GO0006950 -s ${DATASETS}/GO_0006950.txt
	$(foreach i, $(subset1), python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090936_$(i).bed -r SRR5090934.bed  -o SRR5090936_GO0006950_$(i) -s ${DATASETS}/GO_0006950.txt -p 0;)
	$(foreach i, $(subset2), python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090936_$(i).bed -r SRR5090934.bed  -o SRR5090936_GO0006950_$(i) -s ${DATASETS}/GO_0006950.txt -p 0;)


#Gene Length 
#---------------
SRR5090937_yeast_0_2000_regfile.log:
	python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed -o SRR5090937_yeast_4000_6000 -s yeast_4000_6000.txt 
	python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed -o SRR5090937_yeast_6000_100000 -s yeast_6000_100000.txt
	python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed -o SRR5090937_yeast_2000_4000 -s yeast_2000_4000.txt
	python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed -o SRR5090937_yeast_0_2000 -s yeast_0_2000.txt


#----
#GO
#-----
#GO0031990
#----------
SRR5090937_GO0031990_regfile.log:
	python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed  -o SRR5090937_GO0031990 -s ${DATASETS}/GO_0031990.txt
	$(foreach i, $(subset1), python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090937_$(i).bed -r SRR5090935.bed  -o SRR5090937_GO0031990_$(i) -s ${DATASETS}/GO_0031990.txt -p 0;)
	$(foreach i, $(subset2), python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090937_$(i).bed -r SRR5090935.bed  -o SRR5090937_GO0031990_$(i) -s ${DATASETS}/GO_0031990.txt -p 0;)

#GO0009651
#-----------
SRR5090937_GO0009651_regfile.log:
	python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed  -o SRR5090937_GO0009651 -s ${DATASETS}/GO_0009651.txt
	$(foreach i, $(subset1), python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090937_$(i).bed -r SRR5090935.bed  -o SRR5090937_GO0009651_$(i) -s ${DATASETS}/GO_0009651.txt -p 0;)
	$(foreach i, $(subset2), python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090937_$(i).bed -r SRR5090935.bed  -o SRR5090937_GO0009651_$(i) -s ${DATASETS}/GO_0009651.txt -p 0;)

#GO0006950
#-----------
SRR5090937_GO0006950_regfile.log:
	python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed  -o SRR5090937_GO0006950 -s ${DATASETS}/GO_0006950.txt
	$(foreach i, $(subset1), python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090937_$(i).bed -r SRR5090935.bed  -o SRR5090937_GO0006950_$(i) -s ${DATASETS}/GO_0006950.txt -p 0;)
	$(foreach i, $(subset2), python ${SRC}/ribofilio.py -t yeast.fa -f SRR5090937_$(i).bed -r SRR5090935.bed  -o SRR5090937_GO0006950_$(i) -s ${DATASETS}/GO_0006950.txt -p 0;)


#-------------------------------------------------------
#Dataset GSE134152
#-------------------------------------------------------

s21 = SRR9670814 #mRNA SD sample 1 #GSM3938051
s22 = SRR9670815 #mRNA GR sample 2 #GSM3938052
s23 = SRR9670816 #footprint SD sample 3-1 #GSM3938053
s24 = SRR9670817 #footprint SD sample 3-2 #GSM3938054
s25 = SRR9670818 #footprint GR sample 4-1 #GSM3938055
s26 = SRR9670819 #footprint GR sample 4-2 #GSM3938056
s27 = SRR9670820  #mRNA SD sample SD 5 #GSM3938057
s28 = SRR9670821  #mRNA sample GR  6 #GSM3938058
s29 = SRR9670822  #footprint SD  sample 7 #GSM3938059
s30 = SRR9670823  #footprint GR sample 8 #GSM3938060


GSE134152 = ${s21} ${s22} ${s23} ${s24} ${s25} ${s26} ${s27} ${s28} ${s29} ${s30}

#-------------------------------------
#Trim, align, and create bed file
#-------------------------------------
galore/SRR9670823_trimmed.fastq.gz:
	$(foreach i, $(GSE134152), trim_galore --gzip --retain_unpaired --trim1 -a "CTGTAGGCACCATCAAT" --fastqc --fastqc_args "--outdir ${WHERE}/fastqc" -o ${WHERE}/galore  ${WHERE}/$(i).fastq ; )

SRR9670823.bowtie2.bam: galore/SRR9670823_trimmed.fastq.gz
	$(foreach i, $(GSE134152), bowtie2 -x yeast -U galore/$(i)_trimmed.fq.gz  -S $(i).bowtie2.sam;)
	$(foreach i, $(GSE134152), samtools view -S -b $(i).bowtie2.sam > $(i).bowtie2.bam;)

SRR9670823.bed: SRR9670823.bowtie2.bam
	$(foreach i, $(GSE134152), bedtools bamtobed -i $(i).bowtie2.bam > $(i).bed;)


#---------
#Shuffle
#---------
SRR9670816_1000.bed:	
	$(foreach i, $(subset1), python src/shuffle.py -t yeast.fa -b SRR9670816.bed  -o SRR9670816_$(i).bed;)
	$(foreach i, $(subset2), python src/shuffle.py -t yeast.fa -b SRR9670816.bed  -o SRR9670816_$(i).bed;)

SRR9670817_1000.bed:     
	$(foreach i, $(subset1), python src/shuffle.py -t yeast.fa -b SRR9670817.bed  -o SRR9670817_$(i).bed;)
	$(foreach i, $(subset2), python src/shuffle.py -t yeast.fa -b SRR9670817.bed  -o SRR9670817_$(i).bed;)

SRR9670818_1000.bed:     
	$(foreach i, $(subset1), python src/shuffle.py -t yeast.fa -b SRR9670818.bed  -o SRR9670818_$(i).bed;)
	$(foreach i, $(subset2), python src/shuffle.py -t yeast.fa -b SRR9670818.bed  -o SRR9670818_$(i).bed;)

SRR9670822_1000.bed:     
	$(foreach i, $(subset1), python src/shuffle.py -t yeast.fa -b SRR9670822.bed  -o SRR9670822_$(i).bed;)
	$(foreach i, $(subset2), python src/shuffle.py -t yeast.fa -b SRR9670822.bed  -o SRR9670822_$(i).bed;)


SRR9670823_1000.bed:     
	$(foreach i, $(subset1), python src/shuffle.py -t yeast.fa -b SRR9670823.bed  -o SRR9670823_$(i).bed;)
	$(foreach i, $(subset2), python src/shuffle.py -t yeast.fa -b SRR9670823.bed  -o SRR9670823_$(i).bed;)


#Main SRR9670822
#------------------
SRR9670822_regfile.log:
	python  ${SRC}/ribofilio.py -t yeast.fa -f SRR9670822.bed -r SRR9670820.bed -o SRR9670822 
	$(foreach i, $(subset1), python  ${SRC}/ribofilio.py -t yeast.fa -f SRR9670822_$(i).bed -r SRR9670820.bed -o SRR9670822 -p 0;)
	$(foreach i, $(subset2), python  ${SRC}/ribofilio.py -t yeast.fa -f SRR9670822_$(i).bed -r SRR9670820.bed -o SRR9670822 -p 0;)

#Main SRR9670823
#------------------
SRR9670823_regfile.log:
	python  ${SRC}/ribofilio.py -t yeast.fa -f SRR9670823.bed -r SRR9670821.bed -o SRR9670823
	$(foreach i, $(subset1), python  ${SRC}/ribofilio.py -t yeast.fa -f SRR9670823_$(i).bed -r SRR9670821.bed -o SRR9670823 -p 0;)
	$(foreach i, $(subset2), python  ${SRC}/ribofilio.py -t yeast.fa -f SRR9670823_$(i).bed -r SRR9670821.bed -o SRR9670823 -p 0;)


#Yeast Genes 
#------------
SRR9670822_yeast_4000_6000_regfile.log:
	python ${SRC}/ribofilio.py -t yeast.fa -f SRR9670822.bed -r SRR9670820.bed -o SRR9670822_yeast_4000_6000 -s yeast_4000_6000.txt
	$(foreach i, $(subset1), python  ${SRC}/ribofilio.py  -t yeast.fa -f SRR9670822.bed -r SRR9670820.bed  -o SRR9670822_yeast_4000_6000 -s yeast_4000_6000.txt -p 0 ;) 
	$(foreach i, $(subset2), python  ${SRC}/ribofilio.py  -t yeast.fa -f SRR9670822.bed -r SRR9670820.bed  -o SRR9670822_yeast_4000_6000 -s yeast_4000_6000.txt -p 0 ;) 

SRR9670822_yeast_6000_10000_regfile.log:
	python ${SRC}/ribofilio.py -t yeast.fa -f SRR9670822.bed -r SRR9670820.bed -o SRR9670822_yeast_6000_100000 -s yeast_6000_100000.txt 
	$(foreach i, $(subset1), python  ${SRC}/ribofilio.py  -t yeast.fa -f SRR9670822.bed -r SRR9670820.bed  -o SRR9670822_yeast_6000_100000 -s yeast_6000_100000.txt -p 0;) 
	$(foreach i, $(subset2), python  ${SRC}/ribofilio.py  -t yeast.fa -f SRR9670822.bed -r SRR9670820.bed  -o SRR9670822_yeast_6000_100000 -s yeast_6000_100000.txt -p 0;) 

SRR9670822_yeast_2000_4000_regfile.log:
	python ${SRC}/ribofilio.py -t yeast.fa -f SRR9670822.bed -r SRR9670820.bed -o SRR9670822_yeast_2000_4000 -s yeast_2000_4000.txt
	$(foreach i, $(subset1), python  ${SRC}/ribofilio.py  -t yeast.fa -f SRR9670822.bed -r SRR9670820.bed -o SRR9670822_yeast_2000_4000 -s yeast_2000_4000.txt -p 0;)
	$(foreach i, $(subset2), python  ${SRC}/ribofilio.py  -t yeast.fa -f SRR9670822.bed -r SRR9670820.bed -o SRR9670822_yeast_2000_4000 -s yeast_2000_4000.txt -p 0;)

SRR9670822_yeast_0_2000_regfile.log:
	python ${SRC}/ribofilio.py -t yeast.fa -f SRR9670822.bed -r SRR9670820.bed -o SRR9670822_yeast_0_2000 -s yeast_0_2000.txt
	$(foreach i, $(subset1), python  ${SRC}/ribofilio.py  -t yeast.fa -f SRR9670822.bed -r SRR9670820.bed -o SRR9670822_yeast_0_2000 -s yeast_0_2000.txt -p 0;)
	$(foreach i, $(subset2), python  ${SRC}/ribofilio.py  -t yeast.fa -f SRR9670822.bed -r SRR9670820.bed -o SRR9670822_yeast_0_2000 -s yeast_0_2000.txt -p 0;)


#----
#GO
#-----
#GO0031990
#----------
SRR9670816_GO0031990_regfile.log: 
	python ${SRC}/ribofilio.py -t yeast.fa -f SRR9670822.bed -r SRR9670820.bed  -o SRR9670822_GO0031990 -s ${DATASETS}/GO_0031990.txt
	$(foreach i, $(subset1), python ${SRC}/ribofilio.py -t yeast.fa -f SRR9670822_$(i).bed -r SRR9670820.bed  -o SRR9670822_GO0031990 -s ${DATASETS}/GO_0031990.txt -p 0;)
	$(foreach i, $(subset2), python ${SRC}/ribofilio.py -t yeast.fa -f SRR9670822_$(i).bed -r SRR9670820.bed  -o SRR9670822_GO0031990 -s ${DATASETS}/GO_0031990.txt -p 0;)

#GO0009651
#-----------
SRR9670822_GO0009651_regfile.log: 
	python ${SRC}/ribofilio.py -t yeast.fa -f SRR9670822.bed -r SRR9670820.bed  -o SRR9670822_GO0009651 -s ${DATASETS}/GO_0009651.txt
	$(foreach i, $(subset1), python ${SRC}/ribofilio.py -t yeast.fa -f SRR9670816_$(i).bed -r SRR9670820.bed  -o SRR9670822_GO0009651 -s ${DATASETS}/GO_0009651.txt -p 0;)
	$(foreach i, $(subset2), python ${SRC}/ribofilio.py -t yeast.fa -f SRR9670816_$(i).bed -r SRR9670820.bed  -o SRR9670822_GO0009651 -s ${DATASETS}/GO_0009651.txt -p 0;)

