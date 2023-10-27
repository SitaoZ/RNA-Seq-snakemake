configfile: "config.yaml"

WORKDIR = config["workdir"]

# samples
FILES = config["samples"]
SAMPLE = FILES.keys()
READ = [1, 2]

# database
REF_BOWTIE_INDEX = config["ref_bowtie_index"]
rRNA_DATABASE = config["rRNA_databases"]
silva_euk_18s = rRNA_DATABASE['silva-euk-18s']
silva_euk_28s = rRNA_DATABASE['silva-euk-28s']
rfam_5_8s = rRNA_DATABASE['rfam-5.8s']
rfam_5s = rRNA_DATABASE['rfam-5s']

USED_GTF = config["USED_GTF"]

# comparison format string for difference gene expression
Ctrl = config['diff_group']['Control']
Treat = config['diff_group']['Treat']
print("Ctrl:", type(Ctrl))
print("Treat:", Treat)
GROUP = "Ctrl:" + ",".join(Ctrl.split(" ")) + "&" + "Treat:" + ",".join(Treat.split(" "))

# adapter
ADAPTER_FORWARD = config['adapter']['forward']
ADAPTER_BACKWARD = config['adapter']['backward']

# software
CUTADAPT = config['software']['cutadapt']
BOWTIE2 = config['software']['bowtie2']
SAMTOOLS = config['software']['samtools']
RSCRIPT = config['software']['Rscript']

rule all:
    input:
        "08.diffplots/diffgene/volcano.diffexp.pdf",
        "08.diffplots/diffgene/volcano.diffexp.filter.pdf",
        "08.diffplots/diffgene/pheatmap.diffexp.filter.pdf",
        "09.enrichment/go.pdf",
        "09.enrichment/kegg.pdf"

def get_fastqc_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]

rule prepare_fq:
    input:
        fq1 = lambda wildcards: FILES[wildcards.sample].split(",")[0],
        fq2 = lambda wildcards: FILES[wildcards.sample].split(",")[1]
    output:
        fq1 = "data/samples/{sample}/{sample}_1.fq.gz",
        fq2 = "data/samples/{sample}/{sample}_2.fq.gz"
    log:
        stdout="logs/prepare_fq.{sample}.stdout",
        stderr="logs/prepare_fq.{sample}.stderr"
    shell:
        "if [ ! -d data/samples ]; then mkdir -p data/samples;fi; ln -s {input.fq1} {output.fq1}; ln -s {input.fq2} {output.fq2};"

rule fastqc:
    input:
        fastq="data/samples/{sample}/{sample}_{read}.fq.gz",
    output:
        out1="01.fastqc/{sample}_{read}_fastqc.html",
        out2="01.fastqc/{sample}_{read}_fastqc.zip"
    threads: 16
    resources:
        tmpdir="tmp"
    message: "fastqc {input.fastq}: {threads} threads"
    benchmark:
        "benchmarks/fastqc.{sample}.{read}.csv"
    log:
        stdout="logs/fastqc.{sample}.{read}.stdout",
        stderr="logs/fastqc.{sample}.{read}.stderr"
    shell:
        "fastqc -t {threads} -o 01.fastqc {input.fastq} > {log.stdout} 2>{log.stderr}"

rule cutadapt:
    input:
        qc1a="01.fastqc/{sample}_1_fastqc.html",
        qc1b="01.fastqc/{sample}_1_fastqc.zip",
        qc2a="01.fastqc/{sample}_2_fastqc.html",
        qc2b="01.fastqc/{sample}_2_fastqc.zip",
        r1="data/samples/{sample}/{sample}_1.fq.gz",
        r2="data/samples/{sample}/{sample}_2.fq.gz"
    output:
        clean1="02.cutadapt/{sample}_1_clean.fq.gz",
        clean2="02.cutadapt/{sample}_2_clean.fq.gz" 
    threads: 16
    resources:
        tmpdir="tmp"
    message: "cutadapt {input.r1} {input.r2} : {threads} threads"
    benchmark:
        "benchmarks/cutadapt.{sample}.csv"
    log:
        stdout="logs/cutadapt.{sample}.stdout",
        stderr="logs/cutadapt.{sample}.stderr"
    params:
        "-m 20 -q 20 --max-n 0.2 -e 0.08 -Z --trim-n"
    shell:
        "{CUTADAPT} {params} -j {threads} -a {ADAPTER_FORWARD} -A {ADAPTER_BACKWARD} -o {output.clean1} -p {output.clean2} {input.r1} {input.r2} > {log.stdout} 2>{log.stderr}"

rule sortmerna:
    input:
        r1="02.cutadapt/{sample}_1_clean.fq.gz",
        r2="02.cutadapt/{sample}_2_clean.fq.gz"
    output:
        out1="03.sortmerna/{sample}/{sample}_non_rRNA_fwd.fq.gz",
        out2="03.sortmerna/{sample}/{sample}_non_rRNA_rev.fq.gz"
    threads: 16
    resources:
        tmpdir = "tmp"
    message: "sortmerna {input.r1} {input.r2} : {threads} threads"
    benchmark:
        "benchmarks/sortmerna.{sample}.csv"
    log:
        stdout="logs/sortmerna.{sample}.stdout",
        stderr="logs/sortmerna.{sample}.stderr"
    params:
        fixed = "--paired_out --out2 --sam --num_alignments 1 -fastx",
        aligned = "03.sortmerna/{sample}/{sample}_rRNA",
        other = "03.sortmerna/{sample}/{sample}_non_rRNA",
        workdir = "03.sortmerna/{sample}",
        idx = "03.sortmerna/{sample}/idx",
        kvdb = "03.sortmerna/{sample}/kvdb",
        readb = "03.sortmerna/{sample}/readb"
    shell:
        "sortmerna {params.fixed} --threads {threads} --workdir {params.workdir} \
         --ref {silva_euk_18s} \
         --ref {silva_euk_28s} \
         --ref {rfam_5_8s} \
         --ref {rfam_5s} \
         --reads {input.r1} --reads {input.r2} --aligned {params.aligned} --other {params.other} -v > {log.stdout} 2>{log.stderr};"
         "rm -rf {params.idx} {params.kvdb} {params.readb}"


rule bowtie2:
    input:
        r1="03.sortmerna/{sample}/{sample}_non_rRNA_fwd.fq.gz",
        r2="03.sortmerna/{sample}/{sample}_non_rRNA_rev.fq.gz"
    output:
        "04.bowtie2/{sample}/{sample}.bam"
    threads: 16 
    resources:
        tmpdir = "tmp"
    message: "bowtie2 {input.r1} {input.r2} : {threads} threads"
    benchmark:
        "benchmarks/bowtie2.{sample}.csv"
    log:
        stdout="logs/bowtie2.{sample}.stdout",
        stderr="logs/bowtie2.{sample}.stderr"
    params:
        "-q --phred33 --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 \
         --score-min L,0,-0.1 -I 1 -X 1000 --no-mixed --no-discordant -k 200"
    shell:
        "{BOWTIE2} {params} -p {threads} -x {REF_BOWTIE_INDEX} \
         -1 {input.r1} -2 {input.r2} 2>{log.stderr} | {SAMTOOLS} view -b -o {output} -"

rule rsem:
    input:
        "04.bowtie2/{sample}/{sample}.bam"
    output:
        gene="05.rsem/{sample}/{sample}.genes.results",
        isoform="05.rsem/{sample}/{sample}.isoforms.results"
    threads: 16
    resources:
        tmpdir = "tmp"
    log:
        stdout="logs/rsem.{sample}.stdout",
        stderr="logs/rsem.{sample}.stderr"
    shell:
        "rsem-calculate-expression --paired-end \
                               --alignments \
                               -p {threads} \
                               {input} \
                               {REF_BOWTIE_INDEX} \
                               05.rsem/{wildcards.sample}/{wildcards.sample} > {log.stdout} 2>{log.stderr}"

rule gene_exp:
    input:
        "05.rsem/{sample}/{sample}.genes.results"
    output:
        "06.genexp/{sample}/{sample}.gene.xls"
    threads: 1
    resources:
        tmpdir = "tmp"
    log:
        stdout="logs/genexp.{sample}.stdout",
        stderr="logs/genexp.{sample}.stderr"
    shell:
        """awk '{{ if($7!=0.00) print $1"\t"$2"\t"$3"\t"$5"\t"$7}}'  {input} > {output}"""

rule diff_exp:
    input:
        expand("06.genexp/{sample}/{sample}.gene.xls", sample=SAMPLE)
    output:
        "07.diffgene/diffgene/FPKM.xls",
        "07.diffgene/diffgene/Ctrl-vs-Treat_DESeq2.diffexp.xls",
        "07.diffgene/diffgene/Ctrl-vs-Treat_DESeq2.diffexp.filter.xls"
    params:
        group={GROUP}
    shell:
        "python scripts/deseq2.py -i 06.genexp -d '{params.group}' -l 1 -p 0.05 -w {WORKDIR} -o 07.diffgene/diffgene "

rule diff_stat:
    input:
        expression="07.diffgene/diffgene/FPKM.xls",
        total="07.diffgene/diffgene/Ctrl-vs-Treat_DESeq2.diffexp.xls",
        filter="07.diffgene/diffgene/Ctrl-vs-Treat_DESeq2.diffexp.filter.xls"
    output:
        vol_total="08.diffplots/diffgene/volcano.diffexp.pdf",
        vol_filter="08.diffplots/diffgene/volcano.diffexp.filter.pdf",
        phe_filter="08.diffplots/diffgene/pheatmap.diffexp.filter.pdf"
    log:
        v_stdout="logs/diff_stat.volcano.stdout",
        v_stderr="logs/diff_stat.volcano.stderr",
        p_stdout="logs/diff_stat.pheatmap.stdout",
        p_stderr="logs/diff_stat.pheatmap.stderr"
    shell:
        "{RSCRIPT} scripts/volcano.R -i {input.total} -o {output.vol_total} >{log.v_stdout} 2>{log.v_stderr};"
        "{RSCRIPT} scripts/volcano.R -i {input.filter} -o {output.vol_filter} >{log.v_stdout} 2>{log.v_stderr};"
        "{RSCRIPT} scripts/pheatmap.R -i {input.expression} -d {input.filter} -o {output.phe_filter} >{log.p_stdout} 2>{log.p_stderr}"

rule enrichment:
    input:
        "07.diffgene/diffgene/Ctrl-vs-Treat_DESeq2.diffexp.filter.xls"
    output:
        go="09.enrichment/go.pdf",
        kegg="09.enrichment/kegg.pdf"
    log:
        go_stdout="logs/enrichment.go.stdout",
        go_stderr="logs/enrichment.go.stderr",
        kegg_stdout="logs/enrichment.kegg.stdout",
        kegg_stderr="logs/enrichment.kegg.stderr"
    shell:
        "{RSCRIPT} scripts/go.R -i {input} -o {output.go} >{log.go_stdout} 2>{log.go_stderr};"
        "{RSCRIPT} scripts/kegg.R -i {input} -o {output.kegg} >{log.kegg_stdout} 2>{log.kegg_stderr}"
