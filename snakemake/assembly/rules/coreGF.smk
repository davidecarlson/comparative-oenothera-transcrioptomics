DB=config['coreGF']

rule blastx:
    input:
        transcripts=RESULTS + "/assemblies/{sample}_trinity.Trinity.fasta",
    output:
        results=RESULTS + "/coreGF/{sample}.blastx.outfmt6"
    threads:
        5
    params:
        blastDB=DB,
        max=20,
        fmt=6,
        evalue=1e-05
    log:
        RESULTS + "/logs/coreGF/{sample}.blastx.log"
    shell:
        "blastx -query {input.transcripts} -db {params.blastDB} "
        "-num_threads {threads} -max_target_seqs {params.max} "
        "-outfmt {params.fmt} -evalue {params.evalue} > {output.results} 2> {log}"

rule coreGF:
    input:
        blast=RESULTS + "/coreGF/{sample}.blastx.outfmt6"
    output:
        coreGFs=RESULTS + "/coreGF/{sample}.coreGFscore.txt"
    params:
        table="/datahome/public/PLAZA_coreGF/coreGF_plaza2.5_rosids.txt"
    log:
        RESULTS + "/logs/coreGF/{sample}.coreGFs.log"
    shell:
        "python /datahome/public/PLAZA_coreGF/coreGF_plaza2.5_geneset.py {params.table} {input.blast} > {output.coreGFs} &> {log}"
