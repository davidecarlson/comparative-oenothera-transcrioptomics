rule stats:
    input:
        filtered=RESULTS + "/assemblies/{sample}.highestiso.Trinity.fasta",
        assembly=RESULTS + "/assemblies/{sample}_trinity.Trinity.fasta"
    output:
        stats_highestiso=RESULTS + "/assemblies/{sample}.highestiso.stats.txt",
        stats_all=RESULTS + "/assemblies/{sample}.stats.txt"
    shell:
        "/home/progs/trinityrnaseq-v2.11.0/util/TrinityStats.pl {input.assembly} > {output.stats_all} "
        "&& /home/progs/trinityrnaseq-v2.11.0/util/TrinityStats.pl {input.filtered} > {output.stats_highestiso}"
