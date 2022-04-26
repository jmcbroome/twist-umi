configfile: "config.yaml"
rule all:
    input:
        "{sample}_consensus_aligned.bam"

rule index:
    output:
        "{config[reference][fasta]}.bwt"
    shell:
        "{config[dependencies][bwa]} index {config[reference][fasta]}"

rule convert_to_bam:
    input:
        "{sample}_1.fastq",
        "{sample}_2.fastq"
    output:
        "{sample}_unmapped.bam"
    shell:
        "java -jar {config[dependencies][picard]} FastqToSam -F1={input[0]} -F2={input[1]} -O={output}"

rule extract_umi:
    input:
        "{sample}_unmapped.bam"
    output:
        "{sample}_tagged.bam"
    shell:
        "java -jar {config[dependencies][fgbio]} ExtractUmisFromBam -i {input} -o {output}"

rule convert_to_fastq:
    input: 
        "{sample}_unmapped.bam"
    output: 
        "{sample}_undone_1.fastq",
        "{sample}_undone_2.fastq"
    shell: 
        "java -jar {config[dependencies][picard]} SamToFastq -F={output[0]} -F2={output[1]} -I={input}"

rule bwa_map:
    input:
        "{config[reference][fasta]}",
        "{config[reference][fasta]}.bwt",
        "{sample}_undone_1.fastq",
        "{sample}_undone_2.fastq"
    output:
        "{sample}_aligned.bam"
    shell:
        "{config[dependencies][bwa]} mem {input[0]} {input[2]} {input[3]} | {config[dependencies][samtools]} view -b > {output}"

rule merge_bams:
    input:
        "{sample}_tagged.bam",
        "{sample}_aligned.bam"
    output:
        "{sample}_merged.bam"
    shell:
        "java -jar {config[dependencies][picard]} MergeBamAlignment -UNMAPPED={input[0]} -ALIGNED={input[1]} -O {output}"

rule markdup:
    input:
        "{sample}_merged.bam"
    output:
        "{sample}_markdup.bam",
        "{sample}_dup_metrics.txt"
    shell:
        "java -jar {config[dependencies][picard]} MarkDuplicates -I={input} -O={output[0]} -M={output[1]}"

rule group_reads:
    input:
        "{sample}_markdup.bam"
    output:
        "{sample}_grouped.bam"
    shell:
        "java -jar {config[dependencies][fgbio]} GroupReadsByUmi -i {input} -o {output} -f {sample}_family_size.txt -s paired"

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
        "java -jar {config[dependencies][picard]} SamToFastq -F={output} -I={input}"

rule bwa_map_consensus:
    input:
        "{config[reference][fasta]}",
        "{config[reference][fasta]}.bwt",
        "{sample}_consensus.fastq"
    output:
        "{sample}_consensus_aligned.bam"
    shell:
        "{config[dependencies][bwa]} mem {input[0]} {input[2]} | {config[dependencies][samtools]} view -b > {output}"
