rule trinity:
    input:
        fwd=RESULTS + "/fastq_trimmed/{sample}.1.trimmed.fastq",
        rev=RESULTS + "/fastq_trimmed/{sample}.2.trimmed.fastq"
    output:
        assembly=RESULTS + "/assemblies/{sample}_trinity/Trinity.fasta",
        map=RESULTS + "/assemblies/{sample}_trinity/Trinity.fasta.gene_trans_map"
    log:
        RESULTS + "/logs/trinity/{sample}.trinity.log"
    params:
        cleanup="--full_cleanup",
        mem="--max_memory 50G"
    threads: 8
    wrapper:
        "0.63.0/bio/trinity"
