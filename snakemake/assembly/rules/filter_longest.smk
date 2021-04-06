rule filter_transcript_longest:
    input:
        assembly=RESULTS + "/assemblies/{sample}_trinity.Trinity.fasta"
    output:
        filtered=RESULTS + "/assemblies/{sample}.longestiso.Trinity.fasta"
    log:
        RESULTS + "/logs/filter/{sample}.filter_longest.log"

    shell:
        "/home/progs/trinityrnaseq-v2.11.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl "
        "{input.assembly} > {output.filtered} 2> {log}"
