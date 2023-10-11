# APAIQ

## generate input files

APAIQ use bedGraph files as input, it can either take two bedGraph files from the forward strand and reverse strand, or single bedGraph file without strand information

#### generate from mapping tools

The bedGraph file could be generated during the alighment using STAR with the option '--outWigType bedGraph --outWigNorm RPM'

#### generate from bam file

The bedGraph could be also obtained from bam file by using genomeCoverageBed from BEDTOOLS 

#### examples to generate bedGraph from bam file 

`genomeCoverageBed -ibam file.bam -split -bg -scale scaled.factor`\
To get the value of reads per million (RPM), instead of raw count, a scaled factor should be specified.\
The factor is equal to 1000000/(total unique mapped reads), and it should be caution to use '-du' if the data is paired-end.

## install APAIQ 
To install the compiled version of APAIQ from conda, we recommend to creat a conda enviroment firstly using\ 
`conda create -n apaiq_env` and then\
`conda install -c joshuachou apaiq`.\
All the enviroment files and required libraries would be installed automatically   


To run APAIQ using source code from Github, please create a enviroment using the provided env files.\
`conda create --name apaiq_env --file apaiq.env.txt`

To install each dependency package manually\
`conda create --name apaiq_env python=3.7`\
`conda install -c bioconda pybedtools`\
`pip install tensorflow`\
`conda install -c anaconda biopython`

Fast install dependency manually, works for python 3.7 to 3.10 
`conda create --name apaiq_env`\
`pip install pybedtools`\
`pip install tensorflow`\
`pip install biopython`\
`pip install pandas`\

To run with GPU, additional dependency should be installed by using the below code
`conda install -c conda-forge tensorflow-gpu`
to replace 
`pip install tensorflow`\

## run APAIQ

`apaiq --input_file=RNAseq.depth.bedGraph --out_dir=out_dir/ --fa_file=genome_fa --name=sample_id --DB_file polyA.bed --model $model`
a test data, pre-trained model and annotation db_file could be found through the link below:
https://drive.google.com/drive/folders/1KNj-dsh5hCmKI3dyhIsi_OuU6-3mLpBW?usp=sharing

Run APAIQ with source code:
`python src_v2/APAIQ.v.1.2.py --input_file=RNAseq.depth.bedGraph --out_dir=out_dir/ --fa_file=genome_fa --name=sample_id --DB_file polyA.bed --model $model --t 30`
Please use the prefix of the model files as the input of `--model`. For instance for the provided model in Google Drive, the option
should be `--model snu398_model.ckpt`

### Options
	--input_file <bedGraph file>			input bedGraph file from strandless data 

	--input_plus <bedGraph file>			input bedGraph file from forward strand

	--input_minus <bedGraph file>			input bedGraph file from reverse strand 

	--fa_file				fasta file of the genome 

	--name					sample name

	--model					model

	--depth					default=1. For unnormalized input. use --depth=10 for 10 millions total single-end mapped reads, 5 millions total paired-end mapped reads

	--out_dir				default='out_dir'. output directory. 
	
	--RNASeqRCThreshold			default=0.05. Minimum RPM threshold for scaning item

	--window				default=201. Window sizes of the scaning item. Please do not change this value if you used the pre-trained model

	--keep_temp				use --keep_temp='yes', if you want to keep the temporary files.
	

## Running time and memory 
Parallelization was implanted in APAIQ using multiple-processing in Python, for which DNA sequence and RNA-seq coverage across the whole genome were divided  into hundreds of blocks. Thus, the run time and memory usage could be variable using different number of cores/CPU and different size of block (default is 100k bp). 
Usually, 2.5-4 G memory were required for each parallelized core/CPU.

## Examples of input and output 
Examples of input and output can be found under the 'demo' directory. A example of input was shown below.

chrX	20505	20601	0.00673982
chrX	20601	20605	0.0134796
chrX	20605	20701	0.00673982
chrX	20908	20980	0.00673982
chrX	20980	21008	0.0134796
chrX	21008	21080	0.00673982
chrX	21948	22047	0.00673982
chrX	22047	22048	0.0134796

The first column is the chromosome. The 2nd and 3rd column are start and end of a genomic region. The fourth column indicate RNAseq coverage in the corresponding region.
A exmaple of output was shown below.

\#chromosme	start	end	score	id	strand	anno_id	anno_source	distance
chrX	303355	303356	34.55842697620392	chrX:+:16	+	chrX:303356:+:PLCXD1	Gencode	0
chrX	1602523	1602524	44.11139154434204	chrX:+:26	+	chrX:1602520:+:AKAP17A	Gencode	4
chrX	2741311	2741312	35.39809334278107	chrX:+:33	+	chrX:2741309:+:CD99	Gencode	3
chrX	2882817	2882818	43.826815128326416	chrX:+:40	+	chrX:2882818:+:GYG2	3'UTR(M)	0
chrX	9717311	9717312	15.25328278541565	chrX:+:47	+	chrX:9717314:+:TBL1X	3'UTR(M)	-2
chrX	9719649	9719650	16.97154474258423	chrX:+:53	+	chrX:9719652:+:TBL1X	3'UTR(M)	-2
chrX	9719739	9719740	35.48719763755798	chrX:+:54	+	chrX:9719740:+:TBL1X	Gencode	0
chrX	9948362	9948363	16.15381270647049	chrX:+:55	+	chrX:9948359:+:SHROOM2	3'UTR(M)	4

The first three columns indicate the genomic coordinate of the identified PAS. The 4th column indicates the score of the PAS and usually 12 was used as the threshhold to filter reliable PAS, the 5th column is the PAS id and the 6th column is the strand. The 7th column is the cloest PAS from the DB file and the 8th column is the information of PAS from the DB file. The 9th is the distance. 
 


