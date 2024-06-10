Follow the same instructions for installing the required files as the new Nanopore Sequencing Amplicon Pipeline

Prior to running the snakefile you need to manually demultiplex the fastq_pass files

Firstly, unqip the fastq_pass files with ##cd to the fastq_pass folder
- gunzip *.gz

Concatenate the fastq files
- cat *fastq > merged.fastq

Demultiplex and remove sequencing adapters
- porechop -i "merged.fastq" -b "/fastq_trimmed/" --threads 8       ##create a new folder in the data directory and name it fastq_trimmed. Specify this directory in the porechop command

Porechop will by default name the demultiplexed fastq files BC01, BC02 etc. You need to now manually rename the the files in the format "trimmed_barcode01.fastq"

Now run the snakemake file contained in this repository (after configuring the config.json file. This is the same as for the new snakefile). I.E.
- snakemake --cores 2
