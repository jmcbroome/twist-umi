configfile: "config.yaml"
rule all:
    input:
        "{sample}_final.bam"
rule index:
    output:
        expand("{reference}.bwt",reference=config['reference']['fasta'])
    shell:
        "{config[dependencies][bwa]} index {config[reference][fasta]}"

rule convert_to_bam:
    input:
        "{sample}_1.fq.gz",
        "{sample}_2.fq.gz"
    output:
        "{sample}_unmapped.bam"
    shell:
        "java -jar {config[dependencies][picard]} FastqToSam -F1 {input[0]} -F2 {input[1]} -O {output} -SM {wildcards.sample}"

rule extract_umi:
    input:
        "{sample}_unmapped.bam"
    output:
        "{sample}_tagged.bam"
    shell:
        "java -jar {config[dependencies][fgbio]} ExtractUmisFromBam -i {input} -o {output} -r {config[read-structure]}"

rule convert_to_fastq:
    input: 
        "{sample}_unmapped.bam"
    output: 
        "{sample}_undone_1.fastq",
        "{sample}_undone_2.fastq"
    shell: 
        "java -jar {config[dependencies][picard]} SamToFastq -F {output[0]} -F2 {output[1]} -I {input}"

rule bwa_map:
    input:
        expand("{reference}",reference=config['reference']['fasta']),
        expand("{reference}.bwt",reference=config['reference']['fasta']),
        "{sample}_undone_1.fastq",
        "{sample}_undone_2.fastq"
    output:
        "{sample}_aligned.bam"
    shell:
        "{config[dependencies][bwa]} mem {input[0]} {input[2]} {input[3]} | {config[dependencies][sambamba]}-view -b > {output}"

rule merge_bams:
    input:
        "{sample}_tagged.bam",
        "{sample}_aligned.bam"
    output:
        "{sample}_merged.bam"
    shell:
        "java -jar {config[dependencies][picard]} MergeBamAlignment -UNMAPPED {input[0]} -ALIGNED {input[1]} -O {output}"

rule markdup:
    input:
        "{sample}_merged.bam"
    output:
        "{sample}_markdup.bam"
    shell:
        "{config[dependencies][sambamba]}-markdup {input} {output} -p"

rule group_reads:
    input:
        "{sample}_markdup.bam"
    output:
        "{sample}_grouped.bam"
    shell:
        "java -jar {config[dependencies][fgbio]} GroupReadsByUmi -i {input} -o {output} -f {wildcards.sample}_family_size.txt -s paired "

rule call_consensus:
    input:
        "{sample}_grouped.bam"
    output:
        "{sample}_consensus.bam"
    shell:
        "java -jar {config[dependencies][fgbio]} CallDuplexConsensusReads -i {input} -o {output}"

rule convert_consensus_to_fastq:
    input: 
        "{sample}_consensus.bam"
    output: 
        "{sample}_consensus.fastq"
    shell: 
        "java -jar {config[dependencies][picard]} SamToFastq -F {output} -I {input}"

rule bwa_map_consensus:
    input:
        expand("{reference}",reference=config['reference']['fasta']),
        expand("{reference}.bwt",reference=config['reference']['fasta']),
        "{sample}_consensus.fastq"
    output:
        "{sample}_final.bam"
    shell:
        "{config[dependencies][bwa]} mem {input[0]} {input[2]} | {config[dependencies][sambamba]}-view -b > {output}"
