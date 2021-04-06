rule make_matrix:
    input:
        quant=RESULTS + "/quantification/{sample}/quant.sf",
        map=RESULTS + "/assemblies/{sample}_trinity.Trinity.fasta.gene_trans_map"
    output:
        check=RESULTS + "/matrices/{sample}.ok",
        matrix=RESULTS + "/matrices/{sample}_out.isoform.TPM.not_cross_norm"
    log:
        RESULTS + "/logs/matrix/{sample}.matrix.log"
    params:
        method="--est_method salmon",
        prefix=RESULTS + "/matrices/{sample}_out"
    shell:
        "/home/progs/trinityrnaseq-v2.11.0/util/abundance_estimates_to_matrix.pl "
        "{params.method} --gene_trans_map {input.map} --out_prefix {params.prefix} {input.quant} 2> {log} && touch {output.check}"
