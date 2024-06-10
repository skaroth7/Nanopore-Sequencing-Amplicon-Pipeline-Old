import os
import shutil
count = 0

configfile:
    "config.json"


SAMPLES, = glob_wildcards(config['data']+"/fastq_pass/demultiplexed_fastq/BC{id}.fastq")

rule all:
    input:
        expand(config['data']+"/vcf_files/BC{sample}_site_qual.lqual", sample=SAMPLES)
        
    params:
        out = config['data']+"/vcf_files/merged_phased_variants.vcf",
        imp = config['data']+"/vcf_files/"
        
    shell:
        "bcftools merge {params.imp}*.gz > {params.out} --force-samples | cat {params.imp}*.lqual > {params.imp}_merged_qual_score.lqual | cat {params.imp}*.INFO > {params.imp}_merged_total_reads.INFO"
  
        
rule get_tot_reads:
    input:
        config['data']+"/vcf_files/BC{sample}.phased.vcf.gz.tbi",

    params:
        pvcf = config['data']+"/vcf_files/BC{sample}.phased.vcf",
        tot_read = config['data']+"/vcf_files/BC{sample}_total_reads"

    output:
        config['data']+"/vcf_files/BC{sample}_total_reads.INFO",

    shell:
        "vcftools --vcf {params.pvcf} --get-INFO TotalReads --out {params.tot_read}"

rule get_info_site_qual:
    input:
        config['data']+"/vcf_files/BC{sample}_total_reads.INFO",
    
    params:
        pvcf = config['data']+"/vcf_files/BC{sample}.phased.vcf",
        qual = config['data']+"/vcf_files/BC{sample}_site_qual",
        
    output:
        config['data']+"/vcf_files/BC{sample}_site_qual.lqual",
        
    shell:
        "vcftools --vcf {params.pvcf} --site-quality --out {params.qual}" 
        
        
  

rule make_directories:
    output:
        vcf = config['data']+"/vcf_files/",
        align = config['data']+"/alignment_files/"

    shell:
        "mkdir -p {output.vcf}"
        "mkdir -p {output.align}"

rule Alignement:
    input:
        fastq = config['data']+"/fastq_pass/demultiplexed_fastq/BC{sample}.fastq",

    output:
        config['data']+"/alignment_files/BC{sample}.bam"

    params:
        ref = config['reference'],
        qual = "BC{sample}.bam.pdf"

    shell:
        "minimap2 -ax map-ont {params.ref} {input.fastq} --secondary=no | samtools sort -o {output} | qualimap bamqc -bam {output} -outfile {params.qual}"

rule samtools_index:
    input:
        config['data']+"/alignment_files/BC{sample}.bam"

    output:
        config['data']+"/alignment_files/BC{sample}.bam.bai"

    shell:
        "samtools index {input}"


rule nanopolish_index:
    input:
        fastq = config['data']+"/fastq_pass/demultiplexed_fastq/BC{sample}.fastq",
        fast5 = config['fast5'],
    params:
        summary = config['sequencing_summary']
    output:
        config['data']+"/fastq_pass/demultiplexed_fastq/BC{sample}.fastq.index",
    shell:
        "nanopolish index -d {input.fast5} -s {params.summary} {input.fastq}"

rule variant_calling:
    input:
        fastq = config['data']+"/fastq_pass/demultiplexed_fastq/BC{sample}.fastq",
        bam = config['data']+"/alignment_files/BC{sample}.bam",
        bambai = config['data']+"/alignment_files/BC{sample}.bam.bai",
        nanindex = config['data']+"/fastq_pass/demultiplexed_fastq/BC{sample}.fastq.index",

    params:
        ref = config['reference'],
        nanref = config['nanopolish_reference']

    output:
        config['data']+"/vcf_files/BC{sample}.vcf",

    shell:
        "nanopolish variants -g {params.ref} -r {input.fastq} -b {input.bam} --min-candidate-frequency 0.2 -w {params.nanref} --ploidy 2 -o {output} -d 10"

rule phasing:
    input:
        vcf = config['data']+"/vcf_files/BC{sample}.vcf",
        bam = config['data']+"/alignment_files/BC{sample}.bam"
    output:
        config['data']+"/vcf_files/BC{sample}.phased.vcf"
    params:
        ref = config['reference'],
    shell:
        "whatshap phase --reference {params.ref} -o {output} {input.vcf} {input.bam} --tag=PS --ignore-read-groups"

rule zip_vcf:
    input:
        config['data']+"/vcf_files/BC{sample}.phased.vcf"
    output:
        config['data']+"/vcf_files/BC{sample}.phased.vcf.gz",
    shell:
        "bgzip -c {input} > {output}"


rule index_VCF:
    input:
        config['data']+"/vcf_files/BC{sample}.phased.vcf.gz"
    output:
        tbi = config['data']+"/vcf_files/BC{sample}.phased.vcf.gz.tbi",
    shell:
        "tabix -p vcf {input}"

            
        
        
