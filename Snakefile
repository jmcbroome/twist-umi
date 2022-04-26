configfile: "config.yaml"
rule all:
    input:
        "{sample}_consensus.bam"

rule index:
    output:
        "{config[dependencies][bwa]} index {config[reference][fasta]}"
