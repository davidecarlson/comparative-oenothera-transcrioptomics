TRANSDECODER=config['TRANSDECODER']

rule transdecoder_predict:
    input:
        filtered=RESULTS + "/assemblies/{sample}.longestiso.Trinity.fasta",
        pep=RESULTS + "/transdecoder/{sample}/longest_orfs.pep",
        pfam=RESULTS + "/hmmer/{sample}.PFAM.out",
        blastp=RESULTS + "/blastp/{sample}.blastp.outfmt6"
    output:
        final_pep=RESULTS + "/transdecoder/{sample}/{sample}.longestiso.Trinity.fasta.transdecoder.pep"
    threads:
        1
    params:
        prog=TRANSDECODER + "/TransDecoder.Predict",
        outdir=RESULTS + "/transdecoder/{sample}",
        indir=RESULTS + "/transdecoder/{sample}",
        pfam="--retain_pfam_hits",
        blastp="--retain_blastp_hits"
    threads:
        1
    log:
        RESULTS + "/logs/transdecoder/{sample}.predict.log"
    shell:
        "cd {params.indir} && "
        "{params.prog} -t {input.filtered} --output_dir {params.outdir} {params.pfam} {input.pfam} {params.blastp} {input.blastp} 2> {log}"
        " && cd -"
