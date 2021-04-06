rule busco:
    input:
        assembly=RESULTS + "/assemblies/{sample}_trinity.Trinity.fasta"
    output:
        check=RESULTS + "/busco/{sample}_out/{sample}.busco.ok",
        summary=RESULTS + "/busco/{sample}_out/run_{sample}_out/short_summary_{sample}_out.txt"
    log:
        RESULTS + "/logs/busco/{sample}.busco.log"
    params:
        db="/datahome/public/dammit_databases/busco2db/embryophyta_odb9",
        outdir=RESULTS + "/busco/{sample}_out",
        prefix="{sample}_out"
    threads: 1
    shell:
        "cd {params.outdir} && "
        "run_busco --in {input.assembly} --force  --cpu {threads} --mode tran "
        "--lineage {params.db} --out {params.prefix} > {log} "
        "&& touch {output.check} "
        "&& cd -"
