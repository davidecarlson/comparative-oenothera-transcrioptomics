rule hmmer:
    input:
        filtered=RESULTS + "/assemblies/{sample}.longestiso.Trinity.fasta",
        pep=RESULTS + "/transdecoder/{sample}/longest_orfs.pep"
    output:
        results=RESULTS + "/hmmer/{sample}.PFAM.out"
    threads:
        2
    params:
        pfam=DB +"/Pfam-A.hmm"
    log:
        RESULTS + "/logs/hmmer/{sample}.hmmer.log"
    shell:
        "hmmscan --cpu {threads} --domtblout {output.results} "
        " {params.pfam} {input.pep} > {log}"
