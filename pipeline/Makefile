BBMAP = /home/sa941/software/bbmap

YEAST = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/annotations/Saccharomyces_cerevisiae/Ensembl/R64-1-1

DATASETS = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/RibosomalProfiling/ribofilio/datasets
SCRIPTS = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/RibosomalProfiling/ribofilio/scripts
CURRENT  =  /rds/project/yhbl2/rds-yhbl2-genehunter/SM/RibosomalProfiling/ribofilio/pipeline

subset1 := $(shell seq 1 350)
subset2 := $(shell seq 351 700)
subset3 := $(shell seq 701 1000)


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
	python ${SCRIPTS}/gcluster.py yeast.fa yeast_0_2000.txt --max 2000 --min 0
	python ${SCRIPTS}/gcluster.py yeast.fa yeast_2000_4000.txt --max 4000 --min 2000
	python ${SCRIPTS}/gcluster.py yeast.fa yeast_4000_6000.txt --max 6000 --min 4000
	python ${SCRIPTS}/gcluster.py yeast.fa yeast_6000_100000.txt --max 100000 --min 6000


#------------------------------------------------------------------------------------------
#Dataset: GSE91068
#------------------------------------------------------------------------------------------
s1 = SRR5090934 #GSM2420486
s2 = SRR5090935 #GSM2420487
s3 = SRR5090936 #GSM2420488
s4 = SRR5090937 #GSM2420489

SAMPLES = ${s1} ${s2} ${s3} ${s4}


#Adpaters reported by dataset author @jiashun@genome.ucsf.edu 
galore/SRR5090937_trimmed.fastq.gz: ${DATASETS}/GSE91068/SRR5090937.fastq
	$(foreach i, $(SAMPLES), trim_galore --gzip --retain_unpaired --trim1 -a "CTGTAGGCACCATCAAT" --fastqc --fastqc_args "--outdir ${CURRENT}/fastqc" -o ${CURRENT}/galore  ${DATASETS}/GSE91068/$(i).fastq ; )

SRR5090937.bowtie2.bam: galore/SRR5090937_trimmed.fastq.gz
	$(foreach i, $(SAMPLES), bowtie2 -x yeast -U galore/$(i)_trimmed.fq.gz  -S $(i).bowtie2.sam;)
	$(foreach i, $(SAMPLES), samtools view -S -b $(i).bowtie2.sam > $(i).bowtie2.bam;)

SRR5090937.bed: 
	$(foreach i, $(SAMPLES), bedtools bamtobed -i $(i).bowtie2.bam > $(i).bed;)


SRR5090937_1000.bed:     
	$(foreach i, $(subset1), python ${SCRIPTS}/shuffle.py -t yeast.fa -b SRR5090937.bed  -o SRR5090937_$(i).bed;) 
	$(foreach i, $(subset2), python ${SCRIPTS}/shuffle.py -t yeast.fa -b SRR5090937.bed  -o SRR5090937_$(i).bed;) 


#SRR5090936
#-------------
SRR5090936_1000.bed: 
	$(foreach i, $(subset1), python ${SCRIPTS}/shuffle.py -t yeast.fa -b SRR5090936.bed  -o SRR5090936_$(i).bed;) 
	$(foreach i, $(subset2), python ${SCRIPTS}/shuffle.py -t yeast.fa -b SRR5090936.bed  -o SRR5090936_$(i).bed;) 


SRR5090936_WLRB_regfile.log:	
	python  ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -o SRR5090936_WLRB  -b 50 -n SRR5090936_WLRB_regfile.log
	$(foreach i, $(subset1), python  ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090936_$(i).bed -r SRR5090934.bed  -b 50 --plots 0 -n SRR5090936_WLRB_regfile.log;)
	$(foreach i, $(subset2), python  ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090936_$(i).bed -r SRR5090934.bed  -b 50 --plots 0 -n SRR5090936_WLRB_regfile.log;)
	$(foreach i, $(subset3), python  ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090936_$(i).bed -r SRR5090934.bed  -b 50 --plots 0 -n SRR5090936_WLRB_regfile.log;)

SRR5090936_WLRB_yeast_0_2000_regfile.log:
	python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -o SRR5090936_yeast_4000_6000 -s yeast_4000_6000.txt -b 50 -n SRR5090936_WLRB_yeast_4000_6000_regfile.log 
	python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -o SRR5090936_yeast_6000_100000 -s yeast_6000_100000.txt -b 50 -n SRR5090936_WLRB_yeast_6000_100000_regfile.log
	python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -o SRR5090936_yeast_2000_4000 -s yeast_2000_4000.txt -b 50 -n SRR5090936_WLRB_yeast_2000_4000_regfile.log 
	python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed -o SRR5090936_yeast_0_2000 -s yeast_0_2000.txt -b 50 -n SRR5090936_WLRB_yeast_0_2000_regfile.log

#----
#GO
#-----
#GO0031990
#----------
SRR5090936_GO0031990_WLRB_regfile.log:
	python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed  -o SRR5090936_GO0031990_WLRB -s ${DATASETS}/GO_0031990.txt  -b 50  -n SRR5090936_GO0031990_WLRB_regfile.log

SRR5090936_GO0031990_WLRB_350.png:
	$(foreach i, $(subset1), python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090936_$(i).bed -r SRR5090934.bed -s ${DATASETS}/GO_0031990.txt --plot 0 -b 50  -n SRR5090936_GO0031990_WLRB_regfile.log ;)

SRR5090936_GO0031990_WLRB_700.png:
	$(foreach i, $(subset2), python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090936_$(i).bed -r SRR5090934.bed -s ${DATASETS}/GO_0031990.txt --plot 0 -b 50  -n SRR5090936_GO0031990_WLRB_regfile.log ;)

SRR5090936_GO0031990_WLRB_1000.png:
	$(foreach i, $(subset3), python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090936_$(i).bed -r SRR5090934.bed  -s ${DATASETS}/GO_0031990.txt --plot 0 -b 50  -n SRR5090936_GO0031990_WLRB_regfile.log ;)

#GO0009651
#-----------
SRR5090936_GO0009651_WLRB_regfile.log:
	python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090936.bed -r SRR5090934.bed  -o SRR5090936_GO0009651_WLRB -s ${DATASETS}/GO_0009651.txt  -b 50  -n SRR5090936_GO0009651_WLRB_regfile.log

SRR5090936_GO0009651_WLRB_350.png:
	$(foreach i, $(subset1), python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090936_$(i).bed -r SRR5090934.bed  -o SRR5090936_GO0009651_WLRB_$(i) -s ${DATASETS}/GO_0009651.txt --plot 0 -b 50  -n SRR5090936_GO0009651_WLRB_regfile.log ;)

SRR5090936_GO0009651_WLRB_700.png:
	$(foreach i, $(subset2), python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090936_$(i).bed -r SRR5090934.bed  -o SRR5090936_GO0009651_WLRB_$(i) -s ${DATASETS}/GO_0009651.txt --plot 0 -b 50  -n SRR5090936_GO0009651_WLRB_regfile.log ;)

SRR5090936_GO0009651_WLRB_1000.png:
	$(foreach i, $(subset3), python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090936_$(i).bed -r SRR5090934.bed  -o SRR5090936_GO0009651_WLRB_$(i) -s ${DATASETS}/GO_0009651.txt --plot 0 -b 50  -n SRR5090936_GO0009651_WLRB_regfile.log ;)



#SRR5090937
#------------
SRR5090937_WLRB_regfile.log:
	python  ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed  -b 50 -o SRR5090937_WLRB -n SRR5090937_WLRB_regfile.log
	$(foreach i, $(subset1), python  ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090937_$(i).bed -r SRR5090935.bed  -b 50 --plots 0 -n SRR5090937_WLRB_regfile.log;)
	$(foreach i, $(subset2), python  ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090937_$(i).bed -r SRR5090935.bed  -b 50 --plots 0 -n SRR5090937_WLRB_regfile.log;)
	$(foreach i, $(subset3), python  ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090937_$(i).bed -r SRR5090935.bed  -b 50 --plots 0 -n SRR5090937_WLRB_regfile.log;)

#Gene Length 
#---------------
SRR5090937_WLRB_yeast_0_2000_regfile.log:
	python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed -o SRR5090937_yeast_4000_6000 -s yeast_4000_6000.txt -b 50 -n SRR5090937_WLRB_yeast_4000_6000_regfile.log 
	python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed -o SRR5090937_yeast_6000_100000 -s yeast_6000_100000.txt -b 50 -n SRR5090937_WLRB_yeast_6000_100000_regfile.log
	python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed -o SRR5090937_yeast_2000_4000 -s yeast_2000_4000.txt -b 50 -n SRR5090937_WLRB_yeast_2000_4000_regfile.log 
	python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed -o SRR5090937_yeast_0_2000 -s yeast_0_2000.txt -b 50 -n SRR5090937_WLRB_yeast_0_2000_regfile.log


#----
#GO
#-----
#GO0031990
#----------
SRR5090937_GO0031990_WLRB_regfile.log:
	python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed  -o SRR5090937_GO0031990_WLRB -s ${DATASETS}/GO_0031990.txt  -b 50  -n SRR5090937_GO0031990_WLRB_regfile.log

SRR5090937_GO0031990_WLRB_350.png:
	$(foreach i, $(subset1), python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090937_$(i).bed -r SRR5090935.bed  -o SRR5090937_GO0031990_WLRB_$(i) -s ${DATASETS}/GO_0031990.txt --plot 0 -b 50  -n SRR5090937_GO0031990_WLRB_regfile.log ;)

SRR5090937_GO0031990_WLRB_700.png:
	$(foreach i, $(subset2), python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090937_$(i).bed -r SRR5090935.bed  -o SRR5090937_GO0031990_WLRB_$(i) -s ${DATASETS}/GO_0031990.txt --plot 0 -b 50  -n SRR5090937_GO0031990_WLRB_regfile.log ;)

SRR5090937_GO0031990_WLRB_1000.png:
	$(foreach i, $(subset3), python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090937_$(i).bed -r SRR5090935.bed  -o SRR5090937_GO0031990_WLRB_$(i) -s ${DATASETS}/GO_0031990.txt --plot 0 -b 50  -n SRR5090937_GO0031990_WLRB_regfile.log ;)

#GO0009651
#-----------
SRR5090937_GO0009651_WLRB_regfile.log:
	python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed  -o SRR5090937_GO0009651_WLRB -s ${DATASETS}/GO_0009651.txt  -b 50  -n SRR5090937_GO0009651_WLRB_regfile.log

SRR5090937_GO0009651_WLRB_350.png:
	$(foreach i, $(subset1), python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090937_$(i).bed -r SRR5090935.bed  -o SRR5090937_GO0009651_WLRB_$(i) -s ${DATASETS}/GO_0009651.txt --plot 0 -b 50  -n SRR5090937_GO0009651_WLRB_regfile.log ;)

SRR5090937_GO0009651_WLRB_700.png:
	$(foreach i, $(subset2), python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090937_$(i).bed -r SRR5090935.bed  -o SRR5090937_GO0009651_WLRB_$(i) -s ${DATASETS}/GO_0009651.txt --plot 0 -b 50  -n SRR5090937_GO0009651_WLRB_regfile.log ;)

SRR5090937_GO0009651_WLRB_1000.png:
	$(foreach i, $(subset3), python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR5090937_$(i).bed -r SRR5090935.bed  -o SRR5090937_GO0009651_WLRB_$(i) -s ${DATASETS}/GO_0009651.txt --plot 0 -b 50  -n SRR5090937_GO0009651_WLRB_regfile.log ;)




#-----------
#Test mode 
#-----------

SRR5090937_testnoW_regfile.log:
	python ${SCRIPTS}/ribofilioWLR.py -t yeast.fa -f SRR5090937.bed -r SRR5090935.bed -o SRR5090937_testnoW -b 50 -n SRR5090937_testnoW_regfile.log

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



SRR9670816_1000.bed:	
	$(foreach i, $(subset1), python ${SCRIPTS}/shuffle.py -t yeast.fa -b SRR9670816.bed  -o SRR9670816_$(i).bed;)
	$(foreach i, $(subset2), python ${SCRIPTS}/shuffle.py -t yeast.fa -b SRR9670816.bed  -o SRR9670816_$(i).bed;)

SRR9670817_1000.bed:     
	$(foreach i, $(subset1), python ${SCRIPTS}/shuffle.py -t yeast.fa -b SRR9670817.bed  -o SRR9670817_$(i).bed;)
	$(foreach i, $(subset2), python ${SCRIPTS}/shuffle.py -t yeast.fa -b SRR9670817.bed  -o SRR9670817_$(i).bed;)

SRR9670818_1000.bed:     
	$(foreach i, $(subset1), python ${SCRIPTS}/shuffle.py -t yeast.fa -b SRR9670818.bed  -o SRR9670818_$(i).bed;)
	$(foreach i, $(subset2), python ${SCRIPTS}/shuffle.py -t yeast.fa -b SRR9670818.bed  -o SRR9670818_$(i).bed;)

SRR9670822_1000.bed:     
	$(foreach i, $(subset1), python ${SCRIPTS}/shuffle.py -t yeast.fa -b SRR9670822.bed  -o SRR9670822_$(i).bed;)
	$(foreach i, $(subset2), python ${SCRIPTS}/shuffle.py -t yeast.fa -b SRR9670822.bed  -o SRR9670822_$(i).bed;)


#Main SRR9670816
#------------------
SRR9670816_WLRB_regfile.log:
	python  ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR9670816.bed -r SRR9670814.bed -o SRR9670816_WLRB -b 50 -n SRR9670816_WLRB_regfile.log
	$(foreach i, $(subset1), python  ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR9670816_$(i).bed -r SRR9670814.bed  -b 50 --plots 0 -n SRR9670816_WLRB_regfile.log;)
	$(foreach i, $(subset2), python  ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR9670816_$(i).bed -r SRR9670814.bed  -b 50 --plots 0 -n SRR9670816_WLRB_regfile.log;)

#Yeast Genes 
#------------
SRR9670816_WLRB_yeast_0_2000_regfile.log:
	python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR9670816.bed -r SRR9670814.bed -o SRR9670816_yeast_4000_6000 -s yeast_4000_6000.txt -b 50 -n SRR9670814_WLRB_yeast_4000_6000_regfile.log
	python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR9670816.bed -r SRR9670814.bed -o SRR9670816_yeast_6000_100000 -s yeast_6000_100000.txt -b 50 -n SRR9670814_WLRB_yeast_6000_100000_regfile.log 
	python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR9670816.bed -r SRR9670814.bed -o SRR9670816_yeast_2000_4000 -s yeast_2000_4000.txt -b 50 -n SRR9670814_WLRB_yeast_2000_4000_regfile.log
	python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR9670816.bed -r SRR9670814.bed -o SRR9670816_yeast_0_2000 -s yeast_0_2000.txt -b 50 -n SRR9670814_WLRB_yeast_0_2000_regfile.log

#----
#GO
#-----
#GO0031990
#----------
SRR9670816_GO0031990_WLRB_regfile.log: 
	python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR9670816.bed -r SRR9670814.bed  -o SRR9670816_GO0031990_WLRB -s ${DATASETS}/GO_0031990.txt  -b 50  -n SRR9670816_GO0031990_WLRB_regfile.log

SRR9670816_GO0031990_WLRB_350.png: 
	$(foreach i, $(subset1), python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR9670816_$(i).bed -r SRR9670814.bed  -o SRR9670816_GO0031990_WLRB_$(i) -s ${DATASETS}/GO_0031990.txt --plot 0 -b 50  -n SRR9670816_GO0031990_WLRB_regfile.log ;)

SRR9670816_GO0031990_WLRB_700.png:
	$(foreach i, $(subset2), python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR9670816_$(i).bed -r SRR9670814.bed  -o SRR9670816_GO0031990_WLRB_$(i) -s ${DATASETS}/GO_0031990.txt --plot 0 -b 50  -n SRR9670816_GO0031990_WLRB_regfile.log ;)

SRR9670816_GO0031990_WLRB_1000.png:
	$(foreach i, $(subset3), python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR9670816_$(i).bed -r SRR9670814.bed  -o SRR9670816_GO0031990_WLRB_$(i) -s ${DATASETS}/GO_0031990.txt --plot 0 -b 50  -n SRR9670816_GO0031990_WLRB_regfile.log ;)

#GO0009651
#-----------
SRR9670816_GO0009651_WLRB_regfile.log: 
	python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR9670816.bed -r SRR9670814.bed  -o SRR9670816_GO0009651_WLRB -s ${DATASETS}/GO_0009651.txt  -b 50  -n SRR9670816_GO0009651_WLRB_regfile.log

SRR9670816_GO0009651_WLRB_350.png:
	$(foreach i, $(subset1), python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR9670816_$(i).bed -r SRR9670814.bed  -o SRR9670816_GO0009651_WLRB_$(i) -s ${DATASETS}/GO_0009651.txt --plot 0 -b 50  -n SRR9670816_GO0009651_WLRB_regfile.log ;)

SRR9670816_GO0009651_WLRB_700.png:
	$(foreach i, $(subset2), python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR9670816_$(i).bed -r SRR9670814.bed  -o SRR9670816_GO0009651_WLRB_$(i) -s ${DATASETS}/GO_0009651.txt --plot 0 -b 50  -n SRR9670816_GO0009651_WLRB_regfile.log ;)

SRR9670816_GO0009651_WLRB_1000.png:
	$(foreach i, $(subset3), python ${SCRIPTS}/ribofilioWLRB.py -t yeast.fa -f SRR9670816_$(i).bed -r SRR9670814.bed  -o SRR9670816_GO0009651_WLRB_$(i) -s ${DATASETS}/GO_0009651.txt --plot 0 -b 50  -n SRR9670816_GO0009651_WLRB_regfile.log ;)


#-----------------------------------------
#Dataset GSE81966
#-----------------------------------------

galore/SRR3591816_trimmed.fq.gz:
	trim_galore --gzip --retain_unpaired --trim1  --fastqc --fastqc_args "--outdir ${CURRENT}/fastqc" -o ${CURRENT}/galore  ${DATASETS}/GSE81966/SRR3591812.fastq
	trim_galore --gzip --retain_unpaired --trim1  --fastqc --fastqc_args "--outdir ${CURRENT}/fastqc" -o ${CURRENT}/galore  ${DATASETS}/GSE81966/SRR3591816.fastq 

SRR3591816.bam:	
	bowtie2 -x yeast -U galore/SRR3591812_trimmed.fq.gz -S SRR3591812.bowtie2.sam
	samtools view -S -b SRR3591812.bowtie2.sam > SRR3591812.bowtie2.bam
	bowtie2 -x yeast -U galore/SRR3591816_trimmed.fq.gz -S SRR3591816.bowtie2.sam
	samtools view -S -b SRR3591816.bowtie2.sam > SRR3591816.bowtie2.bam

