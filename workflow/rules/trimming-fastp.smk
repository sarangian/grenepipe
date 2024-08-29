# =================================================================================================
#     Trimming
# =================================================================================================

def unpack_fastp_files(wildcards):
    return list(get_fastq(wildcards).values())

rule trim_reads_se:
    input:
        sample=unpack_fastp_files,
    output:
        trimmed=(
            "trimming/{sample}-{unit}.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimming/{sample}-{unit}.fastq.gz")
        ),
        html="trimming/{sample}-{unit}-se-fastp.html",
        json="trimming/{sample}-{unit}-se-fastp.json",
        done=touch("trimming/{sample}-{unit}-se.done"),
    log:
        "logs/trimming/fastp/{sample}-{unit}.log",
    benchmark:
        "benchmarks/trimming/fastp/{sample}-{unit}.log",
    params:
        extra=config["params"]["fastp"]["se"],
    threads: config["params"]["fastp"]["threads"]
    shell:
        """
        fastp -i {input.sample} -o {output.trimmed} -h {output.html} -j {output.json} {params.extra} > {log} 2>&1
        """

rule trim_reads_pe:
    input:
        sample=unpack_fastp_files,
    output:
        trimmed=(
            ["trimming/{sample}-{unit}.1.fastq.gz", "trimming/{sample}-{unit}.2.fastq.gz"]
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp(
                ["trimming/{sample}-{unit}.1.fastq.gz", "trimming/{sample}-{unit}.2.fastq.gz"]
            )
        ),
        html="trimming/{sample}-{unit}-pe-fastp.html",
        json="trimming/{sample}-{unit}-pe-fastp.json",
        done=touch("trimming/{sample}-{unit}-pe.done"),
    log:
        "logs/trimming/fastp/{sample}-{unit}.log",
    benchmark:
        "benchmarks/trimming/fastp/{sample}-{unit}.log",
    params:
        extra=config["params"]["fastp"]["pe"],
    threads: config["params"]["fastp"]["threads"]
    shell:
        """
        fastp -i {input.sample[0]} -I {input.sample[1]} -o {output.trimmed[0]} -O {output.trimmed[1]} -h {output.html} -j {output.json} {params.extra} > {log} 2>&1
        """

rule trim_reads_pe_merged:
    input:
        sample=unpack_fastp_files,
    output:
        merged=(
            "trimming/{sample}-{unit}-merged.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimming/{sample}-{unit}-merged.fastq.gz")
        ),
        html="trimming/{sample}-{unit}-pe-merged-fastp.html",
        json="trimming/{sample}-{unit}-pe-merged-fastp.json",
        done=touch("trimming/{sample}-{unit}-pe-merged.done"),
    log:
        "logs/trimming/fastp/{sample}-{unit}.log",
    benchmark:
        "benchmarks/trimming/fastp/{sample}-{unit}.log",
    params:
        extra=config["params"]["fastp"]["pe"]
        + " --merge --merged_out trimming/{sample}-{unit}-merged.fastq.gz"
        + " --out1 trimming/{sample}-{unit}-unmerged.pass-1.fastq.gz"
        + " --out2 trimming/{sample}-{unit}-unmerged.pass-2.fastq.gz"
        + " --unpaired1 trimming/{sample}-{unit}-unmerged.unpaired-1.fastq.gz"
        + " --unpaired2 trimming/{sample}-{unit}-unmerged.unpaired-2.fastq.gz",
    threads: config["params"]["fastp"]["threads"]
    shell:
        """
        fastp -i {input.sample[0]} -I {input.sample[1]} {params.extra} -h {output.html} -j {output.json} > {log} 2>&1
        """

# =================================================================================================
#     Trimming Results
# =================================================================================================

def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if is_single_end(wildcards.sample, wildcards.unit):
        # single end sample
        return [
            "trimming/{sample}-{unit}.fastq.gz".format(sample=wildcards.sample, unit=wildcards.unit)
        ]
    elif config["settings"]["merge-paired-end-reads"]:
        # merged paired-end samples
        return [
            "trimming/{sample}-{unit}-merged.fastq.gz".format(
                sample=wildcards.sample, unit=wildcards.unit
            )
        ]
    else:
        # paired-end sample
        return expand(
            "trimming/{sample}-{unit}.{pair}.fastq.gz",
            pair=[1, 2],
            sample=wildcards.sample,
            unit=wildcards.unit,
        )


def get_trimming_report(sample, unit):
    """Get the report needed for MultiQC."""
    if is_single_end(sample, unit):
        # single end sample
        return "trimming/" + sample + "-" + unit + "-se-fastp.json"
    elif config["settings"]["merge-paired-end-reads"]:
        # merged paired-end samples
        return "trimming/" + sample + "-" + unit + "-pe-merged-fastp.json"
    else:
        # paired-end sample
        return "trimming/" + sample + "-" + unit + "-pe-fastp.json"
