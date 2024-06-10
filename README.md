Follow the same instructions for installing the required files as the new Nanopore Sequencing Amplicon Pipeline

Prior to running the snakefile you need to manually demultiplex the fastq_pass files

Firstly, unqip the fastq_pass files with ##cd to the fastq_pass folder
-gunzip *.gz

Concatenate the fastq files
-cat *fastq > merged.fastq

Demultiplex and remove sequencing adapters
-porechop -i "merged.fastq" -b "/demultiplex_fastq/" --threads 8       ##create a new folder in the data directory and name it demultiplex_fastq. Specify this directoy in the porechop command

Now run the snakemake file contained in this page (after configuring the config.json file. This is the same as for the new snakefile). I.E.
--snakemake --cores 2
