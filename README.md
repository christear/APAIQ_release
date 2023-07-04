# APAIQ

## generate input files

APAIQ use bedGraph files as input, it can either take two bedGraph files from the forward strand and reverse strand, or single bedGraph file without strand information

#### generate from mapping tools

The bedGraph file could be generated during the alighment using STAR with the option '--outWigType bedGraph'

#### generate from bam file

The bedGraph could be also obtained from bam file by using genomeCoverageBed from BEDTOOLS 

#### examples to generate bedGraph from bam file 

`genomeCoverageBed -ibam file.bam -split -bg -scale scaled.factor`
To get the value of reads per million (RPM), instead of raw count, a scaled factor should be specified.
The factor is equal to 1000000/(total unique mapped reads) 

## install APAIQ 
To install the compiled version of APAIQ from conda, we recommend to creat a conda enviroment firstly using 
`conda create -n apaiq_env` and then
`conda install -c joshuachou apaiq`.
All the enviroment files and required libraries would be installed automatically   


To run APAIQ using source code from Github, please create a enviroment using the provided env files.
`conda create --name apaiq_env --file apaiq.env.txt`

To install each dependency package manually
`conda create --name apaiq_env python=3.7`\
`conda install -c bioconda pybedtools`\
`pip install tensorflow`\
`conda install -c anaconda biopython`\

## run APAIQ

`apaiq --input_file=RNAseq.depth.bedGraph --out_dir=out_dir/ --fa_file=genome_fa --name=sample_id --DB_file polyA.bed --model $model`
a test data, pre-trained model and annotation db_file could be found through the link below:
https://drive.google.com/drive/folders/1D-I_LN1DXQno8BUXUEEIVu7QUMZuEkWQ?usp=sharing

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
