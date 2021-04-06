
rule fastp:
    input:
        fwd=INPUT + "/{sample}-read_1.fq",
        rev=INPUT + "/{sample}-read_2.fq"
    output:
        fwd=RESULTS + "/fastq_trimmed/{sample}.1.trimmed.fastq",
        rev=RESULTS + "/fastq_trimmed/{sample}.2.trimmed.fastq",
        html=RESULTS + "/fastq_trimmed/{sample}.html",
        json=RESULTS + "/fastq_trimmed/{sample}.json"
    threads:
        5
    log:
        RESULTS + "/logs/fastp/{sample}.preprocess.log"
    shell:
        "fastp --in1 {input.fwd} --in2 {input.rev} "
        "--out1 {output.fwd} --out2 {output.rev} --thread {threads} --cut_tail --html {output.html} --json {output.json} 2> {log}"
