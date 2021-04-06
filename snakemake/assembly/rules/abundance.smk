rule abundance:
    input:
        assembly=RESULTS + "/assemblies/{sample}_trinity.Trinity.fasta",
        fwd=RESULTS + "/fastq_trimmed/{sample}.1.trimmed.fastq",
        rev=RESULTS + "/fastq_trimmed/{sample}.2.trimmed.fastq",
        map=RESULTS + "/assemblies/{sample}_trinity.Trinity.fasta.gene_trans_map"
    output:
        check=RESULTS + "/quantification/{sample}.ok",
        #salmon_dir=directory(RESULTS + "/quantification/{sample}"),
        quant=RESULTS + "/quantification/{sample}/quant.sf"
    log:
        RESULTS + "/logs/quant/{sample}.quant.log"
    params:
        method="--est_method salmon",
        type="--seqType fq",
        prep="--prep_reference",
        #options="--validateMappings",
        outdir=RESULTS + "/quantification/{sample}"
    shell:
        "/home/progs/trinityrnaseq-v2.11.0/util/align_and_estimate_abundance.pl "
        "--transcripts {input.assembly} {params.type} --left {input.fwd} --right {input.rev} "
        "{params.method} --gene_trans_map {input.map} {params.prep} "
        "--output {params.outdir} 2> {log} && touch {output.check}"
