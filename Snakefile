'''
Aim: A Snakemake workflow to process SPIDR paired end data (in-progress)
'''

import os 
import sys
import numpy as np
import datetime
from pathlib import Path


##############################################################################
#Initialize settings
##############################################################################


#Copy config file into logs
v = datetime.datetime.now()
run_date = v.strftime('%Y.%m.%d.')

try:
    config_path = config["config_path"]
except:
    config_path = 'config.yaml'

configfile: config_path

try:
    email = config['email']
except:
    email = None
    print("Will not send email on error")


##############################################################################
#Location of scripts
##############################################################################

barcode_id_jar = "scripts/java/BarcodeIdentification_v1.2.0.jar"
lig_eff = "scripts/python/get_ligation_efficiency.py"
split_bpm_rpm = "scripts/python/split_rpm_bpm_fq.py"
add_chr = "scripts/python/ensembl2ucsc.py"
get_clusters = "scripts/python/get_clusters.py"
merge_clusters = "scripts/python/merge_clusters.py"
label_clusters = "scripts/python/label_clusters.py"
fq_to_bam = "scripts/python/fastq_to_bam.py"
add_RG_to_bam = "scripts/python/add_tag_to_bam.py"
split_fastq = "scripts/split_fastq.sh"
add_chr_bt2 = "scripts/python/add_chr_bt2.py"

##############################################################################
#General settings
##############################################################################

try:
    bid_config = config['bID']
    print('Using BarcodeID config', bid_config)
except:
    bid_config = 'config.txt'
    print('Config "bID" not specified, looking for config at:', bid_config)

try:
    num_tags = config['num_tags']
    print('Using', num_tags, 'tags')
except:
    num_tags = "6"
    print('Config "num_tags" not specified, using:', num_tags)

try:
    assembly = config['assembly']
    assert assembly in ['mm10', 'hg38'], 'Only "mm10" or "hg38" currently supported'
    print('Using', assembly)
except:
    print('Config "assembly" not specified, defaulting to "mm10"')
    assembly = 'mm10'

try:
    samples = config['samples']
    print('Using samples file:', samples)
except:
    samples = './samples.json'
    print('Defaulting to working directory for samples json file')

try:
    out_dir = config['output_dir']
    print('All data will be written to:', out_dir)
except:
    out_dir = ''
    print('Defaulting to working directory as output directory')

try:
    temp_dir = config['temp_dir']
    print("Using temporary directory:", temp_dir)
except:
    temp_dir = '/central/scratch/'
    print('Defaulting to central scratch as temporary directory')

try:
    num_chunks = config['num_chunks']
except:
    num_chunks = 2

##############################################################################
##Trimming Sequences
##############################################################################

try:
    adapters = "-g file:" + config['cutadapt_dpm']
    print('Using cutadapt sequence file', adapters)
except:
    adapters = "-g GGTGGTCTTT -g GCCTCTTGTT \
        -g CCAGGTATTT -g TAAGAGAGTT -g TTCTCCTCTT -g ACCCTCGATT"
    print("No file provided for cutadapt. Using standard cutadapt sequences")

try:
    oligos = "-g file:" + config['cutadapt_oligos']
    print('Using bead oligo file', oligos)
except:
    oligos = "-g GGTGGTCTTT -g GCCTCTTGTT"
    print("Using junk oligos. FIX ME")

################################################################################
#Aligner Indexes
################################################################################

#RNA aligner

try:
    bowtie2_index = config['bowtie2_index'][config['assembly']]
except:
    print('Bowtie2 index not specified in config.yaml')
    sys.exit() #no default, exit

try:
    star_index = config['star_index'][config['assembly']]
except:
    print('Star index not specified in config.yaml')
    sys.exit() #no default, exit

################################################################################
#make output directories (aren't created automatically on cluster)
################################################################################

Path(out_dir + "workup/logs/cluster").mkdir(parents=True, exist_ok=True)
out_created = os.path.exists(out_dir + "workup/logs/cluster")
print('Output logs path created:', out_created)

################################################################################
#Get sample files
###############################################################################

#Prep samples from fastq directory using fastq2json_updated.py, now load json file
FILES = json.load(open(samples))
ALL_SAMPLES = sorted(FILES.keys())

ALL_FASTQ = []
for SAMPLE, file in FILES.items():
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R1')])
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R2')])

CONFIG = [out_dir + "workup/logs/config_" + run_date + "yaml"]

NUM_CHUNKS = [f"{i:03}" for i in np.arange(0, num_chunks)]

################################################################################
#Trimming
################################################################################

SPLIT_FQ = expand(out_dir + "workup/splitfq/{sample}_{read}.part_{splitid}.fastq.gz", sample=ALL_SAMPLES, read = ["R1", "R2"], splitid=NUM_CHUNKS)

TRIM = expand([out_dir + "workup/trimmed/{sample}_R1.part_{splitid}_val_1.fq.gz", out_dir + "workup/trimmed/{sample}_R2.part_{splitid}_val_2.fq.gz"], sample = ALL_SAMPLES, splitid=NUM_CHUNKS)
TRIM_LOG = expand(out_dir + "workup/trimmed/{sample}_{read}.part_{splitid}.fastq.gz_trimming_report.txt", sample = ALL_SAMPLES, read = ["R1", "R2"], splitid = NUM_CHUNKS)
TRIM_RD = expand([out_dir + "workup/trimmed/{sample}_R1.part_{splitid}.barcoded_rpm.RDtrim.fastq.gz",
                  out_dir + "workup/trimmed/{sample}_R1.part_{splitid}.barcoded_bpm.RDtrim.fastq.gz",
                  out_dir + "workup/trimmed/{sample}_R2.part_{splitid}.barcoded_rpm.RDtrim.fastq.gz"],
                  sample = ALL_SAMPLES, splitid=NUM_CHUNKS)

R2_RPM = expand(out_dir + "workup/fastqs/{sample}_R2.part_{splitid}.barcoded_rpm.fastq.gz", sample = ALL_SAMPLES, splitid = NUM_CHUNKS)

################################################################################
#Logging
################################################################################

LE_LOG_ALL = [out_dir + "workup/ligation_efficiency.txt"]

MULTI_QC = [out_dir + "workup/qc/multiqc_report.html"]

################################################################################
#Barcoding
################################################################################

BARCODEID = expand(out_dir + "workup/fastqs/{sample}_{read}.part_{splitid}.barcoded.fastq.gz", sample = ALL_SAMPLES, read = ["R1", "R2"], splitid=NUM_CHUNKS)

SPLIT_RPM_BPM = expand([out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_bpm.fastq.gz",
                    out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_rpm.fastq.gz"], sample=ALL_SAMPLES, splitid=NUM_CHUNKS)

################################################################################
#RNA workup
################################################################################

BT2_RNA_ALIGN = expand([out_dir + "workup/alignments/{sample}.part_{splitid}.bowtie2.sorted.mapped.bam",
                        out_dir + "workup/alignments/{sample}.part_{splitid}.bowtie2.sorted.unmapped.bam",
                        out_dir + "workup/alignments/{sample}.part_{splitid}.bowtie2.sorted.mapped.bam.bai"],
                        sample=ALL_SAMPLES, splitid=NUM_CHUNKS)

BAM_TO_FQ = expand([out_dir + "workup/fastqs/{sample}.part_{splitid}.bowtie2.unmapped_R1.fq.gz",
                    out_dir + "workup/fastqs/{sample}.part_{splitid}.bowtie2.unmapped_R2.fq.gz"],
                    sample=ALL_SAMPLES, splitid=NUM_CHUNKS)

STAR_ALIGN = expand(out_dir + "workup/alignments/{sample}.part_{splitid}.Aligned.out.sorted.bam", sample=ALL_SAMPLES, splitid=NUM_CHUNKS)

UNALIGNED = expand([out_dir + "workup/unmapped/{sample}_R1.part_{splitid}.unaligned.fastq.gz",
                    out_dir + "workup/unmapped/{sample}_R2.part_{splitid}.unaligned.fastq.gz"],
                    sample=ALL_SAMPLES, splitid=NUM_CHUNKS)

MERGE_RNA = expand(out_dir + "workup/alignments/{sample}.merged.RPM.bam", sample = ALL_SAMPLES)

CHR_RPM = expand([out_dir + "workup/alignments/{sample}.part_{splitid}.Aligned.out.sorted.chr.bam",
                  out_dir + "workup/alignments/{sample}.part_{splitid}.bowtie2.sorted.mapped.chr.bam"], 
                  sample=ALL_SAMPLES, splitid=NUM_CHUNKS)

################################################################################
#Bead workup
################################################################################

FQ_TO_BAM = expand(out_dir + "workup/alignments/{sample}.part_{splitid}.BPM.bam", sample=ALL_SAMPLES, splitid=NUM_CHUNKS)

MERGE_BEAD = expand(out_dir + "workup/alignments/{sample}.merged.BPM.bam", sample=ALL_SAMPLES)

################################################################################
#Clustering
################################################################################

CLUSTERS = expand(out_dir + "workup/clusters/{sample}.part_{splitid}.clusters", sample=ALL_SAMPLES, splitid=NUM_CHUNKS)
CLUSTERS_MERGED = expand(out_dir + "workup/clusters/{sample}.clusters", sample=ALL_SAMPLES)

################################################################################
################################################################################
#RULE ALL
################################################################################
################################################################################

rule all:
    input: CLUSTERS + CLUSTERS_MERGED + MERGE_RNA

#Send and email if an error occurs during execution
onerror:
    shell('mail -s "an error occurred" ' + email + ' < {log}')

wildcard_constraints:
    sample = "[^\.]+"


##################################################################################
#Trimming and barcode identification
##################################################################################

rule splitfq:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample]['R2']
    output:
        temp(expand([(out_dir + "workup/splitfq/{{sample}}_R1.part_{splitid}.fastq"), (out_dir + "workup/splitfq/{{sample}}_R2.part_{splitid}.fastq")],  splitid=NUM_CHUNKS))
    params:
        dir = out_dir + "workup/splitfq",
        prefix_r1 = "{sample}_R1.part_0",
        prefix_r2 = "{sample}_R2.part_0"
    log:
        out_dir + "workup/logs/{sample}.splitfq.log"
    conda:
        "envs/sprite.yaml"
    threads: 
        8
    shell:
        '''
        mkdir -p {params.dir}
        bash {split_fastq} {input.r1} {num_chunks} {params.dir} {params.prefix_r1}
        bash {split_fastq} {input.r2} {num_chunks} {params.dir} {params.prefix_r2}
        '''


rule compress_fastq:
    input:
        r1 = out_dir + "workup/splitfq/{sample}_R1.part_{splitid}.fastq",
        r2 = out_dir + "workup/splitfq/{sample}_R2.part_{splitid}.fastq"
    output:
        r1 = out_dir + "workup/splitfq/{sample}_R1.part_{splitid}.fastq.gz",
        r2 = out_dir + "workup/splitfq/{sample}_R2.part_{splitid}.fastq.gz"
    conda:
        "envs/pigz.yaml"
    threads:
        8
    shell:
        '''
        pigz -p {threads} {input.r1}
        pigz -p {threads} {input.r2}
        '''        

#Trim adaptors
#multiple cores requires pigz to be installed on the system
rule adaptor_trimming_pe:
    input:
        [out_dir + "workup/splitfq/{sample}_R1.part_{splitid}.fastq.gz", 
         out_dir + "workup/splitfq/{sample}_R2.part_{splitid}.fastq.gz"]
    output:
         out_dir + "workup/trimmed/{sample}_R1.part_{splitid}_val_1.fq.gz",
         out_dir + "workup/trimmed/{sample}_R1.part_{splitid}.fastq.gz_trimming_report.txt",
         out_dir + "workup/trimmed/{sample}_R2.part_{splitid}_val_2.fq.gz",
         out_dir + "workup/trimmed/{sample}_R2.part_{splitid}.fastq.gz_trimming_report.txt"
    threads:
        10
    log:
        out_dir + "workup/logs/{sample}.{splitid}.trim_galore.log"
    conda:
        "envs/sprite.yaml"
    shell:
        '''
        if [[ {threads} -gt 8 ]]
        then
            cores=2
        else
            cores=1
        fi

        trim_galore \
        --paired \
        --gzip \
        --cores $cores \
        --quality 20 \
        --fastqc \
        -o {out_dir}workup/trimmed/ \
        {input} &> {log}
        '''

#Identify barcodes using BarcodeIdentification_v1.2.0.jar
rule barcode_id:
    input:
        r1 = out_dir + "workup/trimmed/{sample}_R1.part_{splitid}_val_1.fq.gz",
        r2 = out_dir + "workup/trimmed/{sample}_R2.part_{splitid}_val_2.fq.gz"
    output:
        r1_barcoded = out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded.fastq.gz",
        r2_barcoded = out_dir + "workup/fastqs/{sample}_R2.part_{splitid}.barcoded.fastq.gz"
    log:
        out_dir + "workup/logs/{sample}.{splitid}.bID.log"
    shell:
        "java -jar {barcode_id_jar} \
        --input1 {input.r1} --input2 {input.r2} \
        --output1 {output.r1_barcoded} --output2 {output.r2_barcoded} \
        --config {bid_config} &> {log}"


#Get ligation efficiency
rule get_ligation_efficiency:
    input:
        r1 = out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded.fastq.gz" 
    output:
        temp(out_dir + "workup/{sample}.part_{splitid}.ligation_efficiency.txt")
    conda:
        "envs/sprite.yaml"
    shell:
        "python {lig_eff} {input.r1} > {output}"


rule cat_ligation_efficiency:
    input:
        expand(out_dir + "workup/{sample}.part_{splitid}.ligation_efficiency.txt", sample=ALL_SAMPLES, splitid=NUM_CHUNKS)
    output:
        out_dir + "workup/ligation_efficiency.txt"
    shell:
        "tail -n +1 {input} > {output}"


rule split_bpm_rpm:
    '''
    split bpm and rpm will also remove incomplete barcodes
    '''
    input:
        out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded.fastq.gz"
    output:
        out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_rpm.fastq.gz",
        out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_bpm.fastq.gz",
        out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_short.fastq.gz"
    log:
        out_dir + "workup/logs/{sample}.{splitid}.BPM_RPM.log"
    conda:
       "envs/sprite.yaml"
    shell:
        "python {split_bpm_rpm} --r1 {input} &> {log}"

rule get_barcoded_rpm_r2:
    input:
        r1 = out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_rpm.fastq.gz",
        r2 = out_dir + "workup/fastqs/{sample}_R2.part_{splitid}.barcoded.fastq.gz"
    output:
        r2 = out_dir + "workup/fastqs/{sample}_R2.part_{splitid}.barcoded_rpm.fastq.gz",
        id = out_dir + "workup/fastqs/{sample}_R2.part_{splitid}_ID.txt"
    log:
        out_dir + "workup/logs/{sample}.{splitid}.RPM.r2barcoded.log"
    threads: 10
    conda:
        "envs/seqkit.yaml"
    shell:
        '''
        gzip -d -c {input.r1} {input.r2} | \
        seqkit seq --name --only-id | \
        sort | uniq -d > {output.id} 

        gzip -d -c {input.r2} | \
        seqkit grep --pattern-file {output.id} | \
        gzip -c > {output.r2}
        '''

################################################################################
#Cutadapt
################################################################################

rule cutadapt_rpm:
    input:
        read1 = out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_rpm.fastq.gz",
        read2 = out_dir + "workup/fastqs/{sample}_R2.part_{splitid}.barcoded_rpm.fastq.gz"
    output:
        r1=out_dir + "workup/trimmed/{sample}_R1.part_{splitid}.barcoded_rpm.RDtrim.fastq.gz",
        r2=out_dir + "workup/trimmed/{sample}_R2.part_{splitid}.barcoded_rpm.RDtrim.fastq.gz",
        qc=out_dir + "workup/trimmed/{sample}_R1.part_{splitid}.barcoded_rpm.RDtrim.qc.txt"
    params:
        adapters_r1 = "-a ATCAGCACTTAGCGTCAG",
        adapters_r2 = "-G CTGACGCTAAGTGCTGAT",
        others = "--minimum-length 20"
    log:
        out_dir + "workup/logs/{sample}.{splitid}.RPM.cutadapt.log"
    threads: 10
    conda:
        "envs/sprite.yaml"
    shell:
        '''
        (cutadapt \
         {params.adapters_r1} \
         {params.adapters_r2} \
         {params.others} \
         -o {output.r1} \
         -p {output.r2} \
         -j {threads} \
         {input.read1} {input.read2} > {output.qc}) &> {log}
        
         #fastqc {output.r1}
         #fastqc {output.r2}
        '''

rule cutadapt_oligo:
    '''
    Trim 9mer oligo sequence from bead barcode
    '''
    input:
        out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_bpm.fastq.gz"
    output:
        fastq=out_dir + "workup/trimmed/{sample}_R1.part_{splitid}.barcoded_bpm.RDtrim.fastq.gz",
        qc=out_dir + "workup/trimmed/{sample}_R1.part_{splitid}.barcoded_bpm.RDtrim.qc.txt"
    params:
        adapters_r1 = oligos
    log:
        out_dir + "workup/logs/{sample}.{splitid}.BPM.cutadapt.log"
    threads: 10
    conda:
        "envs/sprite.yaml"
    shell:
        '''
        (cutadapt \
         {params.adapters_r1} \
         -o {output.fastq} \
         -j {threads} \
         {input} > {output.qc}) &> {log}
        '''

################################################################################
#RNA alignment
################################################################################

rule bowtie2_align:
    input:
        fq1 = out_dir + "workup/trimmed/{sample}_R1.part_{splitid}.barcoded_rpm.RDtrim.fastq.gz",
        fq2 = out_dir + "workup/trimmed/{sample}_R2.part_{splitid}.barcoded_rpm.RDtrim.fastq.gz"
    output:
        bam = temp(out_dir + "workup/alignments/{sample}.part_{splitid}.bowtie2.bam"),
        mapped = out_dir + "workup/alignments/{sample}.part_{splitid}.bowtie2.sorted.mapped.bam",
        unmapped = out_dir + "workup/alignments/{sample}.part_{splitid}.bowtie2.sorted.unmapped.bam",
        index = out_dir + "workup/alignments/{sample}.part_{splitid}.bowtie2.sorted.mapped.bam.bai"
    log:
        out_dir + "workup/logs/{sample}.{splitid}.bt2.log"
    threads: 
        10
    conda:
        "envs/sprite.yaml"
    shell:
        '''
        (bowtie2 \
        -p 10 \
        -t \
        -x {bowtie2_index} \
        -1 {input.fq1} -2 {input.fq2} | \
        samtools view -bS -> {output.bam}) &> {log}
        samtools view -b -f 4 {output.bam} | samtools sort -n -o {output.unmapped}
        samtools view -b -F 4 {output.bam} | samtools sort > {output.mapped}
        samtools index {output.mapped}
        '''

rule bam_to_fq:
    input: out_dir + "workup/alignments/{sample}.part_{splitid}.bowtie2.sorted.unmapped.bam"
    output: 
        r1 = out_dir + "workup/fastqs/{sample}.part_{splitid}.bowtie2.unmapped_R1.fq.gz",
        r2 = out_dir + "workup/fastqs/{sample}.part_{splitid}.bowtie2.unmapped_R2.fq.gz"
    log:
        out_dir + "workup/logs/{sample}.{splitid}.bam2fq.log"
    threads:
        10
    conda:
        "envs/sprite.yaml"
    shell:
        '''
        samtools fastq -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n {input} 
        '''

rule star_align:
    input:
        r1 = out_dir + "workup/fastqs/{sample}.part_{splitid}.bowtie2.unmapped_R1.fq.gz",
        r2 = out_dir + "workup/fastqs/{sample}.part_{splitid}.bowtie2.unmapped_R2.fq.gz"
    output:
        sam = temp(out_dir + "workup/alignments/{sample}.part_{splitid}.Aligned.out.sam"),
        sorted = out_dir + "workup/alignments/{sample}.part_{splitid}.Aligned.out.sorted.bam",
        filtered = temp(out_dir + "workup/alignments/{sample}.part_{splitid}.Aligned.out.bam")
    params:
        STAR_OPTIONS = "--readFilesCommand zcat --alignEndsType EndToEnd --outFilterScoreMin 10 --outFilterMultimapNmax 1 --outFilterMismatchNmax 10 --alignIntronMax 100000 --alignMatesGapMax 1300 --alignIntronMin 80 --alignSJDBoverhangMin 5 --alignSJoverhangMin 8 --chimSegmentMin 20 --alignSJstitchMismatchNmax 5 -1 5 5 --outSAMunmapped Within --outReadsUnmapped Fastx", prefix = out_dir + "workup/alignments/{sample}.part_{splitid}."
    log:
        out_dir + "workup/logs/{sample}.{splitid}.star.log"
    threads:
        10
    conda:
        "envs/sprite.yaml"
    shell:
        '''
        (STAR \
        --genomeDir {star_index} \
        --readFilesIn {input.r1} {input.r2} \
        --runThreadN {threads} {params.STAR_OPTIONS} \
        --outFileNamePrefix {params.prefix}) &> {log}

        samtools view -@ {threads} -bS {output.sam} > {output.filtered}
        samtools sort -@ {threads} -o {output.sorted} {output.filtered}
        samtools index {output.sorted}
        '''

rule compress_unaligned:
    input:
        r1 = out_dir + "workup/alignments/{sample}.part_{splitid}.Unmapped.out.mate1",
        r2 = out_dir + "workup/alignments/{sample}.part_{splitid}.Unmapped.out.mate2"
    output:
        fq1 = out_dir + "workup/unmapped/{sample}_R1.part_{splitid}.unaligned.fastq.gz",
        fq2 = out_dir + "workup/unmapped/{sample}_R2.part_{splitid}.unaligned.fastq.gz"
    conda:
        "envs/pigz.yaml"
    threads:
        8
    shell:
        '''
        pigz -p {threads} {input.r1} {input.r2}
       
        mv {input.r1}.gz {output.fq1}
        mv {input.r2}.gz {output.fq2}
        '''

rule add_chr:
    input:
        star = out_dir + "workup/alignments/{sample}.part_{splitid}.Aligned.out.sorted.bam",
        bt2 = out_dir + "workup/alignments/{sample}.part_{splitid}.bowtie2.sorted.mapped.bam"
    output:
        star = out_dir + "workup/alignments/{sample}.part_{splitid}.Aligned.out.sorted.chr.bam",
        bt2 = out_dir + "workup/alignments/{sample}.part_{splitid}.bowtie2.sorted.mapped.chr.bam"
    log:
        out_dir + "workup/logs/{sample}.{splitid}.add_chr.log",
    conda:
        "envs/sprite.yaml"
    shell:
        '''
        python {add_chr} -i {input.star} -o {output.star} --assembly {assembly} &> {log}
        python {add_chr_bt2} -i {input.bt2} -o {output.bt2} --assembly {assembly} 
        '''

rule merge_rna:
    input:
        bt2 = expand(out_dir + "workup/alignments/{{sample}}.part_{splitid}.bowtie2.sorted.mapped.chr.bam", splitid=NUM_CHUNKS),
        star = expand(out_dir + "workup/alignments/{{sample}}.part_{splitid}.Aligned.out.sorted.chr.bam", splitid=NUM_CHUNKS)
    output:
        out_dir + "workup/alignments/{sample}.merged.RPM.bam"
    conda:
        "envs/sprite.yaml"
    threads:
        8
    log:
        out_dir + "workup/logs/{sample}.merge_bams.log"
    shell:
        '''
        (samtools merge -@ {threads} {output} {input.bt2} {input.star}) &> {log}
        samtools index {output}
        '''

##############################################################################
#Workup Bead Oligo
##############################################################################

rule fastq_to_bam:
    input:
        out_dir + "workup/trimmed/{sample}_R1.part_{splitid}.barcoded_bpm.RDtrim.fastq.gz"
    output:
        sorted = out_dir + "workup/alignments/{sample}.part_{splitid}.BPM.bam",
        bam = temp(out_dir + "workup/alignments/{sample}.part_{splitid}.BPM.unsorted.bam")
    log:
        out_dir + "workup/logs/{sample}.{splitid}.make_bam.log"
    conda:
        "envs/sprite.yaml"
    threads:
        8
    shell:
        '''
        python {fq_to_bam} --input {input} --output {output.bam} --config {bid_config} &> {log}
        samtools sort -@ {threads} -o {output.sorted} {output.bam}

        '''

rule merge_beads:
    input:
        expand(out_dir + "workup/alignments/{{sample}}.part_{splitid}.BPM.bam", splitid = NUM_CHUNKS)
    output:
        out_dir + "workup/alignments/{sample}.merged.BPM.bam"
    conda:
        "envs/sprite.yaml"
    log:
        out_dir + "workup/logs/{sample}.merge_beads.log"
    threads:
        8
    shell:
        '''
        (samtools merge -@ {threads} {output} {input}) >& {log}
        '''

###########################################################################
#Make clusters
###########################################################################

rule make_clusters:
    input:
        rpm=out_dir + "workup/alignments/{sample}.part_{splitid}.Aligned.out.sorted.chr.bam",
        bpm=out_dir + "workup/alignments/{sample}.part_{splitid}.BPM.bam",
        bt2=out_dir + "workup/alignments/{sample}.part_{splitid}.bowtie2.sorted.mapped.chr.bam"
    output:
        unsorted = temp(out_dir + "workup/clusters/{sample}.part_{splitid}.unsorted.clusters"),
        sorted = out_dir + "workup/clusters/{sample}.part_{splitid}.clusters"
    log:
        out_dir + "workup/logs/{sample}.{splitid}.make_clusters.log"
    conda:
        "envs/sprite.yaml"
    shell:
        '''
        (python {get_clusters} \
        -i {input.bpm} {input.bt2} {input.rpm}\
        -o {output.unsorted} \
        -n {num_tags})  &> {log}

        sort -k 1 -T {temp_dir} {output.unsorted} > {output.sorted}
        '''

rule merge_clusters:
    input:
        expand(out_dir + "workup/clusters/{{sample}}.part_{splitid}.clusters", splitid=NUM_CHUNKS)
    output:
        mega = temp(out_dir + "workup/clusters/{sample}.duplicated.clusters"),
        final = out_dir + "workup/clusters/{sample}.clusters"
    conda:
       "envs/sprite.yaml"
    log:
        out_dir + "workup/logs/{sample}.merge_clusters.log"
    shell:
        '''
         sort -k 1 -T {temp_dir} -m {input} > {output.mega}
        (python {merge_clusters} -i {output.mega} -o {output.final}) &> {log}
        '''        

################################################################################
#MultiQC
################################################################################

rule multiqc:
    input:
        expand([out_dir + "workup/clusters/{sample}.clusters"], sample=ALL_SAMPLES) 
    output:
        out_dir + "workup/qc/multiqc_report.html"
    log:
        out_dir + "workup/logs/multiqc.log"
    conda:
        "envs/multiqc.yaml"
    shell:
        "multiqc {out_dir}workup -o {out_dir}workup/qc"

