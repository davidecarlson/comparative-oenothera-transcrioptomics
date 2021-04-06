DB=config['DB']

rule blastp:
    input:
        filtered=RESULTS + "/assemblies/{sample}.longestiso.Trinity.fasta",
        pep=RESULTS + "/transdecoder/{sample}/longest_orfs.pep"
    output:
        results=RESULTS + "/blastp/{sample}.blastp.outfmt6"
    threads:
        5
    params:
        uniprot=DB +"/uniprot_sprot.pep",
        max=1,
        fmt=6,
        evalue=1e-5
    log:
        RESULTS + "/logs/blastp/{sample}.blastp.log"
    shell:
        "blastp -query {input.pep} -db {params.uniprot} "
        "-num_threads {threads} -max_target_seqs {params.max} "
        "-outfmt {params.fmt} -evalue {params.evalue} > {output.results} 2> {log}"
