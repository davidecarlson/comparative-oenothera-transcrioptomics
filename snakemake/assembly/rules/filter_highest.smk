rule filter_transcript_highest:
    input:
        matrix=RESULTS + "/matrices/{sample}_out.isoform.TPM.not_cross_norm",
        map=RESULTS + "/assemblies/{sample}_trinity.Trinity.fasta.gene_trans_map",
        assembly=RESULTS + "/assemblies/{sample}_trinity.Trinity.fasta"
    output:
        filtered=RESULTS + "/assemblies/{sample}.highestiso.Trinity.fasta"
    log:
        RESULTS + "/logs/filter/{sample}.filter.log"
    params:
        method="--highest_iso_only"
    shell:
        "/home/progs/trinityrnaseq-v2.11.0/util/filter_low_expr_transcripts.pl "
        "--matrix {input.matrix} --transcripts {input.assembly} "
        "--gene_to_trans_map {input.map} {params.method} > {output.filtered} 2> {log}"
