shell("set -e;")

import os
import sys
outputs = {}
if not os.path.exists(config['sample_ids']):
    print(f"Specified sample id path {config['sample_ids']} does not exists")
sample_ids = open(config['sample_ids']).read().strip().split('\n')


def get_output(config):
    # get output requested in configure file
    if config['quality_control']:
        outputs["qc-raw"] = expand(config['output_dir'] + "/qc-raw/{sample_id}_{mate}_fastqc.html",
                                   sample_id = sample_ids,mate=["1","2"])

    if config['trim_adapter']:
        outputs["trimmed"] = expand(config['output_dir'] + "/trimmed/{sample_id}_{mate}.fastq.gz", 
                                   sample_id = sample_ids,mate=["1","2"])
        if config['quality_control']:
            outputs["qc-trimmed"] = expand(config['output_dir'] + "/qc-trimmed/{sample_id}_{mate}_fastqc.html",
                                   sample_id = sample_ids,mate=["1","2"])

    if config['remove_unwanted_sequences']:
        outputs["cleaned"] = expand(config['output_dir'] + "/cleaned/{sample_id}_{mate}.fastq.gz", 
                                   sample_id = sample_ids,mate=["1","2"])
        if config['quality_control']:
            outputs["qc-cleaned"] = expand(config['output_dir'] + "/qc-cleaned/{sample_id}_{mate}_fastqc.html",
                                   sample_id = sample_ids,mate=["1","2"])

    if config['genome_mapping']:
        outputs["genome"] = expand(config['output_dir'] + "/bam/{sample_id}/genome.bam", sample_id = sample_ids)

    if config['circRNA_mapping']:
        outputs["circRNA"] = expand(config['output_dir'] + "/bam/{sample_id}/circRNA.bam", sample_id = sample_ids)

    if config['metagenomic_classification']:
        outputs["microbe-counts"] = expand(config['output_dir'] + "/microbe/report/{sample_id}.txt", sample_id = sample_ids)
    if config['count_gene']:
        #outputs["gene-counts"] = expand(config['output_dir'] + "/counts/gene/{sample_id}.txt", sample_id = sample_ids)
        outputs["gene-counts"] = config['output_dir'] + '/counts/matrix/gene.txt'
    if config['count_circRNA']:
        #outputs["circRNA-counts"] = expand(config['output_dir'] + "/counts/circRNA/{sample_id}.txt", sample_id = sample_ids)
        outputs["circRNA-counts"] = config['output_dir'] + '/counts/matrix/circRNA.txt'
    if config["count_editing"]:
        #outputs["editing"] = expand(config['output_dir'] + "/editing/coverage/{sample_id}.txt", sample_id = sample_ids)
        outputs["editing-edited"] = config['output_dir'] + '/editing/matrix/editing-coverage.txt',
        outputs["editing-ref"] =  config['output_dir'] + '/editing/matrix/reference-coverage.txt'
    if config['splicing']:
        outputs['splicing'] = expand(config['output_dir'] + '/splicing/events/{event}.MATS.JC.txt', event = ["A3SS","A5SS","MXE","RI","SE"])
    #outputs["coverage"] = expand(config['output_dir'] + '/bigwig/{sample_id}.bigwig',sample_id = sample_ids)
    if config["APA"]:
        outputs["PDUI"] = config["output_dir"] + '/APA/PDUI.txt'
    return list(outputs.values())
        

rule all:
    input:
        get_output(config)


### Quality control before reads trimming
rule qc_raw_pe:
    input:
        fastq_1 = config["input_dir"] + '/{sample_id}_1.fastq.gz',
        fastq_2 = config["input_dir"] + '/{sample_id}_2.fastq.gz'
    output:
        report_1 = '{outdir}/qc-raw/{sample_id}_1_fastqc.html',
        report_2 = '{outdir}/qc-raw/{sample_id}_2_fastqc.html'
    shell:
        """
        fastqc -o {wildcards.outdir}/qc-raw {input.fastq_1} 
        fastqc -o {wildcards.outdir}/qc-raw {input.fastq_2}
        """


### Reads trimming
rule trimming_pe:
    input:
        fastq_1 = config["input_dir"] + '/{sample_id}_1.fastq.gz',
        fastq_2 = config["input_dir"] + '/{sample_id}_2.fastq.gz'
    output:
        fastq_1 = '{outdir}/trimmed/{sample_id}_1.fastq.gz',
        fastq_2 = '{outdir}/trimmed/{sample_id}_2.fastq.gz',
        report_1 = '{outdir}/log/{sample_id}/trimming_statistics_1.txt',
        report_2 = '{outdir}/log/{sample_id}/trimming_statistics_2.txt',
    params:
        quality = 30,
        adapter_1 = "" if config["adapter_1"] == "auto" else f"-a {config['adapter_1']}", 
        adapter_2 = "" if config["adapter_2"] == "auto" else f"-a2 {config['adapter_2']}"
    threads: config["trimming_threads"]
    log: '{outdir}/log/{sample_id}/trimming.txt'
    shell:
        """
        trim_galore --phred33 --paired  {params.adapter_1} {params.adapter_2} \
        --cores {threads} --quality {params.quality} -o {wildcards.outdir}/trimmed \
         --basename {wildcards.sample_id} {input.fastq_1} {input.fastq_2} > {log} 2>&1
        mv {wildcards.outdir}/trimmed/{wildcards.sample_id}_val_1.fq.gz {output.fastq_1}
        mv {wildcards.outdir}/trimmed/{wildcards.sample_id}_val_2.fq.gz {output.fastq_2}
        mv {wildcards.outdir}/trimmed/{wildcards.sample_id}_1.fastq.gz_trimming_report.txt {output.report_1}
        mv {wildcards.outdir}/trimmed/{wildcards.sample_id}_2.fastq.gz_trimming_report.txt {output.report_2}
        """

### Quality control for cleaned reads
rule qc_trimmed_pe:
    input:
        fastq_1 = '{outdir}/trimmed/{sample_id}_1.fastq.gz',
        fastq_2 = '{outdir}/trimmed/{sample_id}_2.fastq.gz'
    output:
        report_1 = '{outdir}/qc-trimmed/{sample_id}_1_fastqc.html',
        report_2 = '{outdir}/qc-trimmed/{sample_id}_2_fastqc.html'
    shell:
        """
        fastqc -o {wildcards.outdir}/qc-trimmed {input.fastq_1} 
        fastqc -o {wildcards.outdir}/qc-trimmed {input.fastq_2}
        """

rule remove_unwanted_sequences:
    input:
        fastq_1 = '{outdir}/trimmed/{sample_id}_1.fastq.gz',
        fastq_2 = '{outdir}/trimmed/{sample_id}_2.fastq.gz',
    output:
        fastq_1 = '{outdir}/cleaned/{sample_id}_1.fastq.gz',
        fastq_2 = '{outdir}/cleaned/{sample_id}_2.fastq.gz',
        bam = '{outdir}/bam/{sample_id}/unwanted.bam'
    threads: config["mapping_threads"]
    params:
        prefix = config['unwanted_sequences']
    log: '{outdir}/log/{sample_id}/remove-unwanted-sequences.log'
    shell:
        """
        bowtie2 --no-unal -p {threads} -1 {input.fastq_1} -2 {input.fastq_2}  --un-conc-gz {wildcards.outdir}/cleaned/{wildcards.sample_id}_%.fastq.gz --no-discordant --end-to-end -x {params.prefix} 2> {log} | samtools view -b > {output.bam}
        """

### Quality control for cleaned reads
rule qc_cleaned_pe:
    input:
        fastq_1 = '{outdir}/cleaned/{sample_id}_1.fastq.gz', 
        fastq_2 = '{outdir}/cleaned/{sample_id}_2.fastq.gz'
    output:
        report_1 = '{outdir}/qc-cleaned/{sample_id}_1_fastqc.html',
        report_2 = '{outdir}/qc-cleaned/{sample_id}_2_fastqc.html'
    shell:
        """
        fastqc -o {wildcards.outdir}/qc-cleaned {input.fastq_1}
        fastqc -o {wildcards.outdir}/qc-cleaned {input.fastq_2}
        """

### Mapping to human genome 
rule genome_mapping:
    input:
        fastq_1 = '{outdir}/cleaned/{sample_id}_1.fastq.gz',
        fastq_2 = '{outdir}/cleaned/{sample_id}_2.fastq.gz'
    output:
        bam = '{outdir}/bam/{sample_id}/genome.bam',
        fastq_1 = '{outdir}/unmapped/genome/{sample_id}_1.fastq.gz',
        fastq_2 = '{outdir}/unmapped/genome/{sample_id}_2.fastq.gz'
    threads: 4
    params:
        genome = config['genome']
    log: '{outdir}/log/{sample_id}/genome-mapping/Log.out'
    shell:
        """
        STAR --genomeDir {params.genome} \
            --readFilesIn {input.fastq_1} {input.fastq_2} \
            --runThreadN {threads} \
            --outFileNamePrefix {wildcards.outdir}/log/{wildcards.sample_id}/genome-mapping/ \
            --outSAMtype BAM Unsorted \
            --outReadsUnmapped Fastx \
            --readFilesCommand gzip -d -c \
            --outSAMmultNmax 1
        mv {wildcards.outdir}/log/{wildcards.sample_id}/genome-mapping/Aligned.out.bam {output.bam}
        gzip {wildcards.outdir}/log/{wildcards.sample_id}/genome-mapping/Unmapped.out.mate1
        gzip {wildcards.outdir}/log/{wildcards.sample_id}/genome-mapping/Unmapped.out.mate2
        mv {wildcards.outdir}/log/{wildcards.sample_id}/genome-mapping/Unmapped.out.mate1.gz {output.fastq_1}
        mv {wildcards.outdir}/log/{wildcards.sample_id}/genome-mapping/Unmapped.out.mate2.gz {output.fastq_2}
        """

### Mapping to circRNA back-spliced junctions

rule circRNA_mapping:
    input:
        fastq_1 = '{outdir}/unmapped/genome/{sample_id}_1.fastq.gz',
        fastq_2 = '{outdir}/unmapped/genome/{sample_id}_2.fastq.gz',
    output:
        fastq_1 = '{outdir}/unmapped/circRNA/{sample_id}_1.fastq.gz',
        fastq_2 = '{outdir}/unmapped/circRNA/{sample_id}_2.fastq.gz',
        bam = '{outdir}/bam/{sample_id}/circRNA.bam'
    params:
        prefix = config['circRNA']
    threads: 4
    log: '{outdir}/log/{sample_id}/map-circRNA.log'
    shell:
        """
        bowtie2 --no-unal -p {threads} -1 {input.fastq_1} -2 {input.fastq_2}  --un-conc-gz {wildcards.outdir}/unmapped/circRNA/{wildcards.sample_id}_%.fastq.gz --no-discordant --end-to-end -x {params.prefix} 2> {log} | samtools view -b > {output.bam}
        """

### Count gene expression with featurecount
rule count_gene:
    input:
        bam = '{outdir}/bam/{sample_id}/genome.bam',
        annotation = config["genome_annotation"] 
    output:  
        count = '{outdir}/counts/gene/{sample_id}.txt'
    threads:
        2
    params:
        strandness = {'no': 0, 'forward': 1, 'reverse': 2}.get(config['strandness'], 0)
    shell:
        """
        featureCounts -p -s {params.strandness} -t exon -T {threads} -g gene_id -a {input.annotation} -o {output.count} {input.bam}
        """

rule summarize_gene_count:
    input:
        counts = expand('{outdir}/counts/gene/{sample_id}.txt',outdir=config["output_dir"],sample_id=sample_ids),
        sample_ids = config["sample_ids"]
    output:
        matrix = '{outdir}/counts/matrix/gene.txt'
    params:
        indir = '{outdir}/counts/gene',
        sample_ids = config["sample_ids"]
    shell:
        """
         scripts/summarize-table.py -i {params.indir} -o {output.matrix} -s {input.sample_ids} --value-field 6 --row-field 0 --row-name gene_id --first-line 2 --fillna 0 --integer 
        """

### Count reads spanning back-spliced junctions
rule count_circRNA:
    input:
        bam = '{outdir}/bam/{sample_id}/circRNA.bam',
    output:
        count = '{outdir}/counts/circRNA/{sample_id}.txt',
        stats = '{outdir}/counts/circRNA/{sample_id}.stats'
    params:
        strandness = config['strandness']
    shell:
        """
        scripts/count-circRNA.py -b {input.bam} -s {params.strandness} -c {output.count} --stats {output.stats}
        """ 

rule summarize_circRNA_count:
    input:
        counts = expand('{outdir}/counts/circRNA/{sample_id}.txt',outdir=config["output_dir"],sample_id=sample_ids),
        sample_ids = config["sample_ids"] 
    output:
        matrix = '{outdir}/counts/matrix/circRNA.txt'
    params:
        indir = '{outdir}/counts/circRNA'
    shell:
        """
        scripts/summarize-table.py -i {params.indir} -o {output.matrix} -s {input.sample_ids} --value-field 1 --row-field 0 --row-name circrna_id --first-line 1 --fillna 0 --integer
        """

rule sort_genome_bam:
    input:
        bam = '{outdir}/bam/{sample_id}/genome.bam'
    output:
        bam = '{outdir}/bam/{sample_id}/genome.sorted.bam',
        bai = '{outdir}/bam/{sample_id}/genome.sorted.bam.bai'
    threads: 1
    shell:
        """
        samtools sort --threads {threads} {input.bam} > {output.bam}
        samtools index  {output.bam}
        """

rule get_coverage:
    input:
        bam = '{outdir}/bam/{sample_id}/genome.sorted.bam',
        chromsize = "reference/fasta/chrom.size" 
    output:
        bedgraph = "{outdir}/bedgraph/{sample_id}.bedgraph",
        bigwig = "{outdir}/bigwig/{sample_id}.bigwig",
    shell:
        """
        bedtools genomecov -ibam {input.bam} -bg -split | LC_ALL=C sort -k1,1 -k2,2n > {output.bedgraph}
        bedGraphToBigWig {output.bedgraph} {input.chromsize} {output.bigwig} 
        """


### Pileup bam files at kown RNA editing sites

rule pile_up:
    input:
        bam = '{outdir}/bam/{sample_id}/genome.sorted.bam',
        genome = config['genome_fasta'],
        positions = 'reference/editing-sites/REDIportal.txt'
    output:
        pileup = '{outdir}/editing/pileup/{sample_id}.txt'
    shell:
        """
        samtools mpileup -l {input.positions} --reference {input.genome} {input.bam}  | awk '$5~/[ACGTacgt]/{{print}}' > {output.pileup}
        """
        

### Get coverage and A to G conversion at known RNA editing sites
rule count_RNA_editing:
    input:
        pileup = '{outdir}/editing/pileup/{sample_id}.txt'
    output:
        coverage = '{outdir}/editing/coverage/{sample_id}.txt'
    shell:
        """  
        scripts/parse-pileup.py -i {input.pileup} -o {output.coverage}
        """

rule summarize_RNA_editing:
    input:
        counts = expand('{outdir}/editing/coverage/{sample_id}.txt',outdir=config["output_dir"],sample_id=sample_ids),
        sample_ids = config["sample_ids"]
    output:
        reference = '{outdir}/editing/matrix/reference-coverage.txt',
        editing = '{outdir}/editing/matrix/editing-coverage.txt'
    params:
        indir = '{outdir}/editing/coverage'
    shell:
        """
        scripts/summarize-table.py -i {params.indir} -o {output.reference} -s {input.sample_ids} --value-field 4 --row-field 0 --row-name site_id --first-line 1 --fillna 0 --integer
        scripts/summarize-table.py -i {params.indir} -o {output.editing} -s {input.sample_ids} --value-field 5 --row-field 0 --row-name site_id --first-line 1 --fillna 0 --integer
        """


rule splicing_analysis:
    input:
        bam = '{outdir}/bam/{sample_id}/genome.sorted.bam',
        gtf = config['genome_annotation'] 
    output:
        count = '{outdir}/splicing/counts/{sample_id}.rmats'
    params:
        read_length = 150
    shell:
        """
        echo {input.bam} > {wildcards.outdir}/splicing/counts/{wildcards.sample_id}.path.txt 
        rmats.py --b1 {wildcards.outdir}/splicing/counts/{wildcards.sample_id}.path.txt -t paired --gtf {input.gtf} --readLength {params.read_length} --variable-read-length --od {wildcards.outdir}/splicing/counts/{wildcards.sample_id} --task prep --tmp {wildcards.outdir}/splicing/counts/{wildcards.sample_id}-tmp
        mv {wildcards.outdir}/splicing/counts/{wildcards.sample_id}-tmp/*.rmats {output.count}
        rm -r {wildcards.outdir}/splicing/counts/{wildcards.sample_id}
        rm -r {wildcards.outdir}/splicing/counts/{wildcards.sample_id}-tmp
        """

rule summarize_splicing:
    input:
        counts = expand('{outdir}/splicing/counts/{sample_id}.rmats', outdir=config["output_dir"], sample_id=sample_ids),
        gtf = config['genome_annotation']
    output:
        events = expand('{outdir}/splicing/events/{event}.MATS.JC.txt', outdir=config["output_dir"], event = ["A3SS","A5SS","MXE","RI","SE"])
    params:
        read_length = 150,
        outdir = config["output_dir"]
    shell:
        """
        cat {params.outdir}/splicing/counts/*path.txt | awk 'BEGIN{{ORS=",";}}{{print}}' | sed 's/,$/\n/' > {params.outdir}/splicing/bam.path.txt
        rmats.py --b1 {params.outdir}/splicing/bam.path.txt --gtf {input.gtf} -t paired --readLength {params.read_length} --nthread 1 --tmp {params.outdir}/splicing/counts --od {params.outdir}/splicing/events  --task post --statoff 
        """


rule APA_analysis:
    input:
        coverages = expand('{outdir}/bigwig/{sample_id}.bigwig', outdir=config["output_dir"], sample_id=sample_ids),
        utr = "reference/bed/gencode.v38.annotation.3putr.subtracted.selected.bed"
    output:
        PDUI = config["output_dir"] + '/APA/PDUI.txt'
    params:
        sample_ids = " ".join(sample_ids),
        wigdir = config["output_dir"] + '/bigwig',
        outdir = config["output_dir"] + '/APA'
    threads: 4
    shell:
        """
        [ -f {params.outdir}/bigwig.paths.txt ] && rm {params.outdir}/bigwig.paths.txt
        for sample_id in {params.sample_ids} ;do
          echo "{params.wigdir}/${{sample_id}}.bigwig" >> {params.outdir}/bigwig.paths.txt
        done
        scripts/infer-APA.py -bws {params.outdir}/bigwig.paths.txt -b {input.utr} --PDUI {output.PDUI} --njobs {threads}
        """


rule metagenomic_classification:
    input:
        fastq_1 = '{outdir}/unmapped/circRNA/{sample_id}_1.fastq.gz',
        fastq_2 = '{outdir}/unmapped/circRNA/{sample_id}_2.fastq.gz',
        database = config['kraken_database'] 
    output:
        report = '{outdir}/microbe/report/{sample_id}.txt',
        assignment = '{outdir}/microbe/assignment/{sample_id}.txt.gz', 
        fastq_1 = '{outdir}/unmapped/unclassified/{sample_id}_1.fastq.gz',
        fastq_2 = '{outdir}/unmapped/unclassified/{sample_id}_2.fastq.gz'
    threads: 4
    log: '{outdir}/log/{sample_id}/kraken2-classification.txt'
    shell:
        """
        kraken2 --db {input.database} --threads {threads} --unclassified-out {wildcards.outdir}/unmapped/unclassified/{wildcards.sample_id}#.fastq --report {output.report} --paired --use-names {input.fastq_1} {input.fastq_2} 2> {log} | gzip -c > {output.assignment}
        gzip {wildcards.outdir}/unmapped/unclassified/{wildcards.sample_id}_1.fastq
        gzip {wildcards.outdir}/unmapped/unclassified/{wildcards.sample_id}_2.fastq 
        """
