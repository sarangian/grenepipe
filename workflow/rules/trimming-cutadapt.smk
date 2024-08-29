# =================================================================================================
#     Trimming
# =================================================================================================

rule trim_reads_se:
    input:
        unpack(get_fastq),
    output:
        fastq=(
            "trimming/{sample}-{unit}.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimming/{sample}-{unit}.fastq.gz")
        ),
        qc="trimming/{sample}-{unit}.qc-se.txt",
        done=touch("trimming/{sample}-{unit}.se.done"),
    params:
        adapters=config["params"]["cutadapt"]["se"]["adapters"],
        extra=config["params"]["cutadapt"]["se"]["extra"],
    threads: config["params"]["cutadapt"]["threads"]
    log:
        "logs/trimming/cutadapt/{sample}-{unit}.log",
    benchmark:
        "benchmarks/trimming/cutadapt/{sample}-{unit}.log"
    shell:
        """
        cutadapt {params.extra} -a {params.adapters} -o {output.fastq} {input} > {output.qc} 2> {log}
        """

rule trim_reads_pe:
    input:
        unpack(get_fastq),
    output:
        fastq1=(
            "trimming/{sample}-{unit}.1.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimming/{sample}-{unit}.1.fastq.gz")
        ),
        fastq2=(
            "trimming/{sample}-{unit}.2.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimming/{sample}-{unit}.2.fastq.gz")
        ),
        qc="trimming/{sample}-{unit}.qc-pe.txt",
        done=touch("trimming/{sample}-{unit}.pe.done"),
    params:
        adapters=config["params"]["cutadapt"]["pe"]["adapters"],
        extra=config["params"]["cutadapt"]["pe"]["extra"],
    threads: config["params"]["cutadapt"]["threads"]
    log:
        "logs/trimming/cutadapt/{sample}-{unit}.log",
    benchmark:
        "benchmarks/trimming/cutadapt/{sample}-{unit}.log"
    shell:
        """
        cutadapt {params.extra} -a {params.adapters[0]} -A {params.adapters[1]} -o {output.fastq1} -p {output.fastq2} {input[0]} {input[1]} > {output.qc} 2> {log}
        """
