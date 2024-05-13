import os

data_dir = "data/"


SAMPLES = set()
for filename in os.listdir(data_dir):
    if filename.endswith("_1.fastq.gz") or filename.endswith("_2.fastq.gz"):
        sample_name = filename.split('_')[0]
        SAMPLES.add(sample_name)

SAMPLES = sorted(list(SAMPLES))

READ1_FILES = expand("data/{sample}_1.fastq.gz", sample=SAMPLES)

rule all:
    input:
        "tabla_de_conteos.txt",
        expand("fastqc_reports/{sample}_1_fastqc.html", sample=SAMPLES),
        expand("fastqc_reports/{sample}_1_fastqc.zip", sample=SAMPLES)


rule fastqc:
    input:
        read1="data/{sample}_1.fastq.gz"
    output:
        html="fastqc_reports/{sample}_1_fastqc.html",
        zip="fastqc_reports/{sample}_1_fastqc.zip"
    params:
        outdir="fastqc_reports/"
    shell:
        "fastqc {input.read1} --outdir={params.outdir}"


rule mapear_secuencias:
    input:
        read1="data/{sample}_1.fastq.gz",
        read2="data/{sample}_2.fastq.gz"
    output:
        temp("mapeado/{sample}.sam")
    params:
        index="Index/L_vannamei",
        maxins=800,
        threads=workflow.cores
    shell:
        "bowtie2 --maxins {params.maxins} -t -x {params.index} -p {params.threads} --fr \
            -1 {input.read1} -2 {input.read2} > {output}"

rule convertir_sam_a_bam:
    input:
        "mapeado/{sample}.sam"
    output:
        temp("mapeado/{sample}.bam")
    params:
        threads=workflow.cores/2
    shell:
        "samtools view -Sb -@ {params.threads} {input} > {output}"

rule ordenar_bam:
    input:
        "mapeado/{sample}.bam"
    output:
        bam="mapeado/{sample}.sorted.bam"
    params:
        threads=workflow.cores/2
    shell:
        "samtools sort -@ {params.threads} {input} -o {output}"

rule obtener_cobertura:
    input:
        "mapeado/{sample}.sorted.bam"
    output:
        cov="cobertura/{sample}.All_transcripts"
    shell:
        "samtools coverage {input} | awk '{{print $1\"\\t\"$4}}' > {output}"

rule generar_tabla_conteos:
    input:
        expand("cobertura/{sample}.All_transcripts", sample=SAMPLES)
    output:
        "tabla_de_conteos.txt"
    params:
        path="cobertura/",
        ext=".All_transcripts",
        script="/scratch/yahelgb/Project_Lvannamei/scripts/genera_tabla_de_conteos.R"
    shell:
        "Rscript {params.script} {params.path} {params.ext} {output}"
