rule trinity:
    input:
        fwd=RESULTS + "/fastq_trimmed/{sample}.1.trimmed.fastq",
        rev=RESULTS + "/fastq_trimmed/{sample}.2.trimmed.fastq"
    output:
        assembly=RESULTS + "/assemblies/{sample}_trinity.Trinity.fasta",
        map=RESULTS + "/assemblies/{sample}_trinity.Trinity.fasta.gene_trans_map"
    log:
        RESULTS + "/logs/trinity/{sample}.trinity.log"
    params:
        prefix=RESULTS + "/assemblies/{sample}_trinity"
    threads: 12
    shell:
        "/home/progs/trinityrnaseq-v2.11.0/Trinity --seqType fq --max_memory 150G --left {input.fwd} --right {input.rev} "
        "--CPU {threads} --full_cleanup --output {params.prefix} > {log}"
