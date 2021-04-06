TRANSDECODER=config['TRANSDECODER']

rule transdecoder_long:
    input:
        filtered=RESULTS + "/assemblies/{sample}.longestiso.Trinity.fasta",
        map=RESULTS + "/assemblies/{sample}_trinity.Trinity.fasta.gene_trans_map"
    output:
        pep=RESULTS + "/transdecoder/{sample}/longest_orfs.pep"
    threads:
        1
    params:
        prog=TRANSDECODER + "/TransDecoder.LongOrfs",
        outdir=RESULTS + "/transdecoder/{sample}"
    log:
        RESULTS + "/logs/transdecoder/{sample}.longestorfs.log"
    shell:
        "rm -r {params.outdir} "
        "&& {params.prog} -t {input.filtered} --gene_trans_map {input.map} "
        "--output_dir {params.outdir} 2> {log}"
