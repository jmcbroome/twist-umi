configfile: "config.yaml"
rule all:
    input:
        "{sample}_final.bam"

rule bwa_index:
    output:
        expand("{reference}.bwt",reference=config['reference'] + ".fa")
    shell:
        "{config[dependencies][bwa]} index {config[reference]}.fa"

rule picard_index:
    output:
        expand("{reference}.dict",reference=config['reference'])
    shell:
        "java -jar ./picard.jar CreateSequenceDictionary -R {config[reference]}.fa -O {config[reference]}.dict"

rule trim_adapters:
    input:
        "{sample}_1.fq.gz",
        "{sample}_2.fq.gz"    
    output:
        "{sample}_1.trimmed.fq.gz",
        "{sample}_2.trimmed.fq.gz",
        "{sample}_fastp.txt"
    threads:
        24
    shell:
        "{config[dependencies][fastp]} --in1 {input[0]} --in2 {input[1]}  --out1 {output[0]} --out2 {output[1]} --thread {threads} --detect_adapter_for_pe -j /dev/null -h /dev/null 2> {output[2]}"

rule convert_to_bam:
    input:
        "{sample}_1.trimmed.fq.gz",
        "{sample}_2.trimmed.fq.gz"
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
        "java -jar {config[dependencies][fgbio]} ExtractUmisFromBam -i {input} -o {output} -r {config[read_structure_r1]} {config[read_structure_r2]} -t RX"

rule convert_to_fastq:
    input: 
        "{sample}_tagged.bam"
    output: 
        "{sample}_undone_1.fastq",
        "{sample}_undone_2.fastq"
    shell: 
        "java -jar {config[dependencies][picard]} SamToFastq -F {output[0]} -F2 {output[1]} -I {input}"

rule bwa_map:
    input:
        expand("{reference}.fa",reference=config['reference']),
        expand("{reference}.fa.bwt",reference=config['reference']),
        "{sample}_undone_1.fastq",
        "{sample}_undone_2.fastq"
    output:
        "{sample}_aligned.bam"
    threads: 24
    shell:
        "{config[dependencies][bwa]} mem -t {threads} {input[0]} {input[2]} {input[3]} | {config[dependencies][samtools]} view -b > {output}"

rule merge_bams:
    input:
        "{sample}_tagged.bam",
        "{sample}_aligned.bam"
    output:
        "{sample}_merged.bam"
    shell:
        "java -jar {config[dependencies][picard]} MergeBamAlignment -UNMAPPED {input[0]} -ALIGNED {input[1]} -O {output} -REFERENCE_SEQUENCE {config[reference]}.fa"

rule markdup:
    input:
        "{sample}_merged.bam"
    output:
        "{sample}_markdup.bam"
    threads: 24
    log:
        "{sample}_markdup.log"
    shell:
        "java -jar {config[dependencies][picard]} MarkDuplicates -I {input} -O {output} -M {log} --DUPLEX_UMI"

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
    threads: 24
    log:
        "{sample}_consensus.log"
    shell:
        "java -jar {config[dependencies][fgbio]} CallDuplexConsensusReads -i {input} -o {output} --threads {threads} -M 3 0 0 2> {log}"

rule convert_consensus_to_fastq:
    input: 
        "{sample}_consensus.bam"
    output: 
        "{sample}_consensus_r1.fastq",
        "{sample}_consensus_r2.fastq"
    shell: 
        "java -jar {config[dependencies][picard]} SamToFastq -F {output[0]} -F2 {output[1]} -I {input}"

rule bwa_map_consensus:
    input:
        expand("{reference}.fa",reference=config['reference']),
        expand("{reference}.fa.bwt",reference=config['reference']),
        "{sample}_consensus_r1.fastq",
        "{sample}_consensus_r2.fastq"
    output:
        "{sample}_final.bam"
    threads: 24
    shell:
        "{config[dependencies][bwa]} mem -t {threads} {input[0]} {input[2]} {input[3]} | {config[dependencies][samtools]} view -b > {output}"
