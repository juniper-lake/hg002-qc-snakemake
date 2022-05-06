import os
import pandas as pd
from collections import defaultdict


# assign workflow directory name based on location of snakefile
workflow_dir = os.path.dirname(workflow.snakefile)
shell.prefix(f"set -o pipefail; umask 002; ")  # set g+w


# configuration
configfile: f"{workflow_dir}/config.yaml"
sample_sheet = pd.read_csv(config['sample_sheet'], sep='\t', header=0)
output_dir = config['output_dir']
workflow_targets = config['targets']
reference_name = config['reference_name']
reference = config['reference']
reference_index = config['reference_index']
chromosome_lengths = config['chromosome_lengths']
tandem_repeats_bed = config['tandem_repeats_bed']
stratifications = config['stratifications']
sdf = config['sdf']
all_chroms = config['chroms']
deepvariant_version = config['deepvariant_version']
n_shards = config['n_shards']
shards = [f"{x:05}" for x in range(config['n_shards'])]

# build a list of targets
targets = []

                                                            
# create dictionaries for multiple inputs per sample/condition (e.g. multiple movies, benchmarks, etc.)
ubam_dict = defaultdict(dict)
truvari_vcf_dict = defaultdict(dict)
truvari_bed_dict = defaultdict(dict)
happy_vcf_dict = defaultdict(dict)
happy_bed_dict = defaultdict(dict)
whatshap_vcf_dict = defaultdict(dict)
sr_yak_dict = defaultdict(dict)


# populate dictionaries and resource-dependent targets from sample sheet
for index, row in sample_sheet.iterrows():
    for movie in row['movies'].split(','):
        ubam_dict[row['condition_name']][os.path.basename(movie).split('.')[0]] = movie
    # truvari
    if pd.notnull(row['truvari_vcf']):
        for item in row['truvari_vcf'].split(','):
            truvari_vcf_dict[row['condition_name']][item.split(':')[0]] = item.split(':')[1]
            # truvari targets
            targets.append(f"{output_dir}/{row['condition_name']}/truvari/{item.split(':')[0]}/summary.txt")
    if pd.notnull(row['truvari_bed']):
        for item in row['truvari_bed'].split(','):
            truvari_bed_dict[row['condition_name']][item.split(':')[0]] = item.split(':')[1]
    # happy
    if pd.notnull(row['happy_vcf']):
        for item in row['happy_vcf'].split(','):
            happy_vcf_dict[row['condition_name']][item.split(':')[0]] = item.split(':')[1]
            # happy targets
            targets.append(f"{output_dir}/{row['condition_name']}/happy/{item.split(':')[0]}.summary.csv")
    if pd.notnull(row['happy_bed']):
        for item in row['happy_bed'].split(','):
            happy_bed_dict[row['condition_name']][item.split(':')[0]] = item.split(':')[1]
    # whatshap
    if pd.notnull(row['whatshap_vcf']):
        for item in row['whatshap_vcf'].split(','):
            whatshap_vcf_dict[row['condition_name']][item.split(':')[0]] = item.split(':')[1]
            # whatshap targets
            targets.append(f"{output_dir}/{row['condition_name']}/whatshap/{item.split(':')[0]}.phase.eval.tsv")
    # short-read yak db
    if pd.notnull(row['sr_yak']):
        for item in row['sr_yak'].split(','):
            sr_yak_dict[row['condition_name']][item.split(':')[0]] = item.split(':')[1]
            # yak qv targets
            targets.append(f"{output_dir}/{row['condition_name']}/hifiasm/{item.split(':')[0]}.asm.p_ctg.qv.txt")
            targets.append(f"{output_dir}/{row['condition_name']}/hifiasm/{item.split(':')[0]}.asm.a_ctg.qv.txt")

print(targets)

# targets that don't depend on resources
targets.extend([f"{output_dir}/{condition}/{filename}"
                    for condition in ubam_dict.keys()
                    for filename in [
                        "smrtcell_stats/all_movies.read_length_and_quality.tsv",
                        "mosdepth/coverage.mosdepth.summary.txt",
                        "mosdepth/M2_ratio.txt",
                        "mosdepth/gc_coverage.summary.txt",
                        "mosdepth/coverage.thresholds.summary.txt"
                                    ]])
if 'assembly' in workflow_targets:
    targets.extend([f"{output_dir}/{condition}/{filename}"
                        for condition in ubam_dict.keys()
                        for filename in [
                            "hifiasm/asm.p_ctg.fasta.stats.txt",
                            "hifiasm/asm.a_ctg.fasta.stats.txt"
                                        ]])
if 'small_variants' in workflow_targets:
    targets.extend([f"{output_dir}/{condition}/{filename}"
                        for condition in ubam_dict.keys()
                        for filename in [
                        "deepvariant/deepvariant.vcf.stats.txt",
                        "whatshap/deepvariant.phased.tsv",
                        "whatshap/deepvariant.phased.blocklist"
                                        ]])
if 'structural_variants' in workflow_targets:
    targets.extend([f"{output_dir}/{condition}/{filename}"
                        for condition in ubam_dict.keys()
                        for filename in [
                        "pbsv/all_chroms.pbsv.vcf.gz"
                                        ]])

    
# rules
localrules: all, 
            calculate_m2_ratio, 
            calculate_gc_coverage, 
            calculate_coverage_thresholds,
            bcftools_concat_pbsv_vcf, 
            split_deepvariant_vcf_round1, 
            split_deepvariant_vcf_round2, 
            whatshap_bcftools_concat_round1, 
            whatshap_bcftools_concat_round2,
            bgzip_vcf, 
            tabix_vcf


ruleorder: pbmm2_align > samtools_index_bam


rule all:
    input: targets


rule smrtcell_stats:
    input: lambda wc: ubam_dict[wc.condition][wc.movie]
    output: f"{output_dir}/{{condition}}/smrtcell_stats/{{movie}}.read_length_and_quality.tsv"
    log: f"{output_dir}/{{condition}}/logs/smrtcell_stats.{{movie}}.log"
    conda: "envs/smrtcell_stats.yaml"
    shell: f"(python3 {workflow_dir}/scripts/extract_read_length_and_qual.py {{input}} > {{output}}) > {{log}} 2>&1"


rule combine_smrtcell_stats:
    input: lambda wc: expand(f"{output_dir}/{wc.condition}/smrtcell_stats/{{movie}}.read_length_and_quality.tsv", movie=ubam_dict[wc.condition].keys())
    output: f"{output_dir}/{{condition}}/smrtcell_stats/all_movies.read_length_and_quality.tsv"
    log: f"{output_dir}/{{condition}}/logs/combine_smrtcell_stats.log"
    shell: "(cat {input} > {output}) > {log} 2>&1"    


rule pbmm2_align:
    input:
        ref = reference,
        ref_index = reference_index,
        query = lambda wc: ubam_dict[wc.condition][wc.movie]
    output:
        bam = f"{output_dir}/{{condition}}/aligned/{{movie}}.{reference_name}.bam",
        bai = f"{output_dir}/{{condition}}/aligned/{{movie}}.{reference_name}.bam.bai"
    log: f"{output_dir}/{{condition}}/logs/pbmm2_align.{{movie}}.log"
    params:
        condition = lambda wc: wc.condition,
        preset = "CCS",
        extra = "--sort --unmapped -c 0 -y 70",
        loglevel = "INFO"
    threads: 24
    conda: "envs/pbmm2.yaml"
    shell:
        """
        (pbmm2 align --num-threads {threads} \
            --preset {params.preset} \
            --sample {params.condition} \
            --log-level {params.loglevel} \
            {params.extra} \
            {input.ref} \
            {input.query} \
            {output.bam}) > {log} 2>&1
        """


rule samtools_fasta:
    input: lambda wc: ubam_dict[wc.condition][wc.movie]
    output: f"{output_dir}/{{condition}}/fasta/{{movie}}.fasta"
    log: f"{output_dir}/{{condition}}/logs/samtools_fasta.{{movie}}.log"
    threads: 4
    conda: "envs/samtools.yaml"
    shell: "(samtools fasta -@ 3 {input} > {output}) > {log} 2>&1"


rule hifiasm_assemble:
    input: lambda wc: expand(f"{output_dir}/{wc.condition}/fasta/{{movie}}.fasta", movie=ubam_dict[wc.condition].keys())
    output:
        f"{output_dir}/{{condition}}/hifiasm/asm.p_ctg.gfa",
        f"{output_dir}/{{condition}}/hifiasm/asm.a_ctg.gfa",
    log: f"{output_dir}/{{condition}}/logs/hifiasm.log"
    conda: "envs/hifiasm.yaml"
    params: prefix = f"{output_dir}/{{condition}}/hifiasm/asm"
    threads: 48
    shell: "(hifiasm -o {params.prefix} --primary -t {threads} {input}) > {log} 2>&1"


rule gfa2fa:
    input: f"{output_dir}/{{condition}}/hifiasm/asm.{{infix}}.gfa"
    output: f"{output_dir}/{{condition}}/hifiasm/asm.{{infix}}.fasta"
    log: f"{output_dir}/{{condition}}/logs/gfa2fa.{{infix}}.log"
    conda: "envs/gfatools.yaml"
    shell: "(gfatools gfa2fa {input} > {output}) 2> {log}"


rule bgzip_fasta:
    input: f"{output_dir}/{{condition}}/hifiasm/asm.{{infix}}.fasta"
    output: f"{output_dir}/{{condition}}/hifiasm/asm.{{infix}}.fasta.gz"
    log: f"{output_dir}/{{condition}}/logs/bgzip_fasta.{{infix}}.log"
    threads: 4
    conda: "envs/htslib.yaml"
    shell: "(bgzip --threads {threads} {input}) > {log} 2>&1"


rule asm_stats:
    input: 
        fasta = f"{output_dir}/{{condition}}/hifiasm/asm.{{infix}}.fasta.gz",
        ref_index = reference_index
    output: f"{output_dir}/{{condition}}/hifiasm/asm.{{infix}}.fasta.stats.txt"
    log: f"{output_dir}/{{condition}}/logs/asm_stats.{{infix}}.log"
    conda: "envs/k8.yaml"
    shell: f"(k8 {workflow_dir}/scripts/calN50.js -f {{input.ref_index}} {{input.fasta}} > {{output}}) > {{log}} 2>&1"


rule yak_qv:
    input:
        asm = f"{output_dir}/{{condition}}/hifiasm/asm.{{infix}}.fasta.gz",
        sr = lambda wc: sr_yak_dict[wc.condition][wc.version]
    output: f"{output_dir}/{{condition}}/hifiasm/{{version}}.asm.{{infix}}.qv.txt"
    log: f"{output_dir}/{{condition}}/logs/yak_qv.{{infix}}.{{version}}.log"
    threads: 24
    conda: "envs/yak.yaml"
    shell: "(yak qv -t{threads} -p -K3.2g -l100k {input.sr} <(zcat {input.asm}) > {output}) > {log} 2>&1"



rule bgzip_vcf:
    input: f"{output_dir}/{{condition}}/{{prefix}}.vcf"
    output: f"{output_dir}/{{condition}}/{{prefix}}.vcf.gz"
    log: f"{output_dir}/{{condition}}/logs/bgzip_vcf/{{prefix}}.log"
    threads: 2
    conda: "envs/htslib.yaml"
    shell: "(bgzip --threads {threads} {input}) > {log} 2>&1"


rule tabix_vcf:
    input: f"{output_dir}/{{condition}}/{{prefix}}.vcf.gz"
    output: f"{output_dir}/{{condition}}/{{prefix}}.vcf.gz.tbi"
    log: f"{output_dir}/{{condition}}/logs/tabix_vcf/{{prefix}}.log"
    params: "-p vcf"
    conda: "envs/htslib.yaml"
    shell: "tabix {params} {input} > {log} 2>&1"


rule pbsv_discover:
    input:
        bam = rules.pbmm2_align.output.bam,
        bai = rules.pbmm2_align.output.bai,
        tr_bed = tandem_repeats_bed
    output: f"{output_dir}/{{condition}}/pbsv/svsig/{{movie}}.{{chromosome}}.svsig.gz"
    log: f"{output_dir}/{{condition}}/logs/pbsv_discover.{{movie}}.{{chromosome}}.log"
    params:
        extra = "--hifi",
        chromosome = lambda wc: wc.chromosome,
        loglevel = "INFO"
    conda: "envs/pbsv.yaml"
    shell:
        """
        (pbsv discover {params.extra} \
            --log-level {params.loglevel} \
            --region {wildcards.chromosome} \
            --tandem-repeats {input.tr_bed} \
            {input.bam} {output}) > {log} 2>&1
        """


rule pbsv_call:
    input:
        svsigs = lambda wc: [f"{output_dir}/{wc.condition}/pbsv/svsig/{movie}.{wc.chromosome}.svsig.gz" for movie in ubam_dict[wc.condition].keys()],
        ref = reference
    output: f"{output_dir}/{{condition}}/pbsv/chrom_vcfs/{{chromosome}}.pbsv.vcf"
    log: f"{output_dir}/{{condition}}/logs/pbsv_call.{{chromosome}}.log"
    params:
        extra = "--hifi -m 20",
        loglevel = "INFO"
    threads: 8
    conda: "envs/pbsv.yaml"
    shell:
        """
        (pbsv call {params.extra} \
            --log-level {params.loglevel} \
            --num-threads {threads} \
            {input.ref} {input.svsigs} {output}) > {log} 2>&1
        """


rule bcftools_concat_pbsv_vcf:
    input:
        calls = lambda wc: expand(f"{output_dir}/{wc.condition}/pbsv/chrom_vcfs/{{chromosome}}.pbsv.vcf.gz", chromosome=all_chroms),
        indices = lambda wc: expand(f"{output_dir}/{wc.condition}/pbsv/chrom_vcfs/{{chromosome}}.pbsv.vcf.gz.tbi", chromosome=all_chroms)
    output: f"{output_dir}/{{condition}}/pbsv/all_chroms.pbsv.vcf"
    log: f"{output_dir}/{{condition}}/logs/bcftools_concat_pbsv_vcf.log"
    conda: "envs/bcftools.yaml"
    shell: "(bcftools concat -a -o {output} {input.calls}) > {log} 2>&1"


rule truvari_benchmark:
    input:
        ref = reference,
        query_vcf = f"{output_dir}/{{condition}}/pbsv/all_chroms.pbsv.vcf.gz",
        query_tbi = f"{output_dir}/{{condition}}/pbsv/all_chroms.pbsv.vcf.gz.tbi",
        bench_vcf = lambda wc: truvari_vcf_dict[wc.condition][wc.version],
    output: 
        f"{output_dir}/{{condition}}/truvari/{{version}}/tp-call.vcf",
        f"{output_dir}/{{condition}}/truvari/{{version}}/tp-base.vcf",
        f"{output_dir}/{{condition}}/truvari/{{version}}/fn.vcf",
        f"{output_dir}/{{condition}}/truvari/{{version}}/fp.vcf",
        f"{output_dir}/{{condition}}/truvari/{{version}}/base-filter.vcf",
        f"{output_dir}/{{condition}}/truvari/{{version}}/call-filter.vcf",
        f"{output_dir}/{{condition}}/truvari/{{version}}/summary.txt",
        f"{output_dir}/{{condition}}/truvari/{{version}}/log.txt",
        f"{output_dir}/{{condition}}/truvari/{{version}}/giab_report.txt"
    log: f"{output_dir}/{{condition}}/logs/truvari_benchmark.{{version}}.log"
    params: 
        prefix = f"{output_dir}/{{condition}}/truvari/{{version}}",
        bed = lambda wc: "--includebed "+truvari_bed_dict[wc.condition][wc.version] if wc.version in truvari_bed_dict[wc.condition].keys() else ""
    conda: "envs/truvari.yaml"
    shell:
        """
        (rmdir {params.prefix} && \
        truvari \
            -f {input.ref} \
            -b {input.bench_vcf} {params.bed} \
            -o {params.prefix} \
            --passonly --giabreport \
            -r 1000 -p 0.00 \
            -c {input.query_vcf}) > {log} 2>&1
        """


rule deepvariant_make_examples_round1:
    input:
        bams = lambda wc: [f"{output_dir}/{{condition}}/aligned/{movie}.{reference_name}.bam" for movie in ubam_dict[wc.condition].keys()],
        bais = lambda wc: [f"{output_dir}/{{condition}}/aligned/{movie}.{reference_name}.bam.bai" for movie in ubam_dict[wc.condition].keys()],
        ref = reference
    output:
        tfrecord = f"{output_dir}/{{condition}}/deepvariant_intermediate/examples/examples.tfrecord-{{shard}}-of-{n_shards:05}.gz"
    log: f"{output_dir}/{{condition}}/logs/deepvariant_make_examples_round1.{{shard}}-of-{n_shards:05}.log"
    container: f"docker://google/deepvariant:{deepvariant_version}"
    params:
        vsc_min_fraction_indels = "0.12",
        pileup_image_width = 199,
        shard = lambda wc: wc.shard,
        reads = lambda wc: ','.join([f"{output_dir}/{wc.condition}/aligned/{movie}.{reference_name}.bam" for movie in ubam_dict[wc.condition].keys()])
    resources:
        extra = '--constraint=avx512'
    shell:
        f"""
        (/opt/deepvariant/bin/make_examples \
            --norealign_reads \
            --vsc_min_fraction_indels {{params.vsc_min_fraction_indels}} \
            --pileup_image_width {{params.pileup_image_width}} \
            --alt_aligned_pileup=diff_channels \
            --add_hp_channel \
            --sort_by_haplotypes \
            --parse_sam_aux_fields \
            --mode calling \
            --ref {{input.ref}} \
            --reads {{params.reads}} \
            --examples {output_dir}/{{wildcards.condition}}/deepvariant_intermediate/examples/examples.tfrecord@{n_shards}.gz \
            --task {{wildcards.shard}}) > {{log}} 2>&1
        """


rule deepvariant_call_variants_gpu_round1:
    input: lambda wc: expand(f"{output_dir}/{wc.condition}/deepvariant_intermediate/examples/examples.tfrecord-{{shard}}-of-{n_shards:05}.gz", shard=shards)
    output: f"{output_dir}/{{condition}}/deepvariant_intermediate/call_variants_output.tfrecord.gz"
    log: f"{output_dir}/{{condition}}/logs/deepvariants_call_variants_gpu_round1.log"
    container: f"docker://google/deepvariant:{deepvariant_version}"
    params: model = "/opt/models/pacbio/model.ckpt"
    threads: 8
    resources:
        partition = 'ml',
        extra = '--gpus=1'
    shell:
        f"""
        (/opt/deepvariant/bin/call_variants \
            --outfile {{output}} \
            --examples {output_dir}/{{wildcards.condition}}/deepvariant_intermediate/examples/examples.tfrecord@{n_shards}.gz \
            --checkpoint {{params.model}}) > {{log}} 2>&1
        """


rule deepvariant_postprocess_variants_round1:
    input:
        tfrecord = rules.deepvariant_call_variants_gpu_round1.output,
        ref = reference
    output:
        vcf = f"{output_dir}/{{condition}}/deepvariant_intermediate/deepvariant.vcf.gz",
        vcf_index = f"{output_dir}/{{condition}}/deepvariant_intermediate/deepvariant.vcf.gz.tbi",
        report = f"{output_dir}/{{condition}}/deepvariant_intermediate/deepvariant.visual_report.html"
    log: f"{output_dir}/{{condition}}/logs/deepvariant_postprocess_variants_round1.log"
    container: f"docker://google/deepvariant:{deepvariant_version}"
    threads: 4
    resources: 
        extra = '--constraint=avx512'
    shell:
        """
        (/opt/deepvariant/bin/postprocess_variants \
            --ref {input.ref} \
            --infile {input.tfrecord} \
            --outfile {output.vcf}) > {log} 2>&1
        """


rule split_deepvariant_vcf_round1:
    input: rules.deepvariant_postprocess_variants_round1.output.vcf,
    output: f"{output_dir}/{{condition}}/whatshap_intermediate/{{chromosome}}.deepvariant.vcf"
    log: f"{output_dir}/{{condition}}/logs/split_deepvariant_vcf_round1.{{chromosome}}.log"
    params: extra = '-h'
    conda: "envs/htslib.yaml"
    shell: "tabix {params.extra} {input} {wildcards.chromosome} > {output} 2> {log}"


rule whatshap_phase_round1:
    input:
        ref = reference,
        vcf = f"{output_dir}/{{condition}}/whatshap_intermediate/{{chromosome}}.deepvariant.vcf.gz",
        tbi = f"{output_dir}/{{condition}}/whatshap_intermediate/{{chromosome}}.deepvariant.vcf.gz.tbi",
        phaseinput = lambda wc: [f"{output_dir}/{{condition}}/aligned/{movie}.{reference_name}.bam" for movie in ubam_dict[wc.condition].keys()],
        phaseinputindex = lambda wc: [f"{output_dir}/{{condition}}/aligned/{movie}.{reference_name}.bam.bai" for movie in ubam_dict[wc.condition].keys()],
    output: f"{output_dir}/{{condition}}/whatshap_intermediate/{{chromosome}}.deepvariant.phased.vcf.gz"
    log: f"{output_dir}/{{condition}}/logs/whatshap_phase_round1.{{chromosome}}.log"
    params: chromosome = lambda wc: wc.chromosome
    conda: "envs/whatshap.yaml"
    shell:
        """
        (whatshap phase \
            --chromosome {wildcards.chromosome} \
            --output {output} \
            --reference {input.ref} \
            {input.vcf} {input.phaseinput}) > {log} 2>&1
        """


rule whatshap_bcftools_concat_round1:
    input:
        calls = lambda wc: expand(f"{output_dir}/{wc.condition}/whatshap_intermediate/{{chromosome}}.deepvariant.phased.vcf.gz", chromosome=all_chroms),
        indices = lambda wc: expand(f"{output_dir}/{wc.condition}/whatshap_intermediate/{{chromosome}}.deepvariant.phased.vcf.gz.tbi", chromosome=all_chroms)
    output: f"{output_dir}/{{condition}}/whatshap_intermediate/deepvariant.phased.vcf.gz"
    log: f"{output_dir}/{{condition}}/logs/whatshap_bcftools_concat_round1.log"
    params: "-a -Oz"
    conda: "envs/bcftools.yaml"
    shell: "bcftools concat {params} -o {output} {input.calls} > {log} 2>&1"


rule whatshap_haplotag_round1:
    input:
        ref = reference,
        vcf = rules.whatshap_bcftools_concat_round1.output,
        tbi = f"{output_dir}/{{condition}}/whatshap_intermediate/deepvariant.phased.vcf.gz.tbi",
        bam = rules.pbmm2_align.output.bam
    output: f"{output_dir}/{{condition}}/whatshap_intermediate/{{movie}}.deepvariant.haplotagged.bam"
    log: f"{output_dir}/{{condition}}/logs/whatshap_haplotag_round1.{{movie}}.log"
    params: "--tag-supplementary"
    conda: "envs/whatshap.yaml"
    shell:
        """
        (whatshap haplotag {params} \
            --output {output} \
            --reference {input.ref} \
            {input.vcf} {input.bam}) > {log} 2>&1
        """


rule samtools_index_bam:
    input: f"{output_dir}/{{condition}}/{{folder}}/{{prefix}}.bam"
    output: f"{output_dir}/{{condition}}/{{folder}}/{{prefix}}.bam.bai"
    log: f"{output_dir}/{{condition}}/logs/samtools_index_bam/{{folder}}/{{prefix}}.log"
    threads: 4
    conda: "envs/samtools.yaml"
    shell: "(samtools index -@ 3 {input}) > {log} 2>&1"


rule deepvariant_make_examples_round2:
    input:
        bams = lambda wc: [f"{output_dir}/{{condition}}/whatshap_intermediate/{movie}.deepvariant.haplotagged.bam" for movie in ubam_dict[wc.condition].keys()],
        bais = lambda wc: [f"{output_dir}/{{condition}}/whatshap_intermediate/{movie}.deepvariant.haplotagged.bam.bai" for movie in ubam_dict[wc.condition].keys()],
        ref = reference
    output:
        tfrecord = f"{output_dir}/{{condition}}/deepvariant/examples/examples.tfrecord-{{shard}}-of-{n_shards:05}.gz",
        nonvariant_site_tfrecord = f"{output_dir}/{{condition}}/deepvariant/examples/gvcf.tfrecord-{{shard}}-of-{n_shards:05}.gz"
    log: f"{output_dir}/{{condition}}/logs/deepvariant_make_examples_round2.{{shard}}-of-{n_shards:05}.log"
    container: f"docker://google/deepvariant:{deepvariant_version}"
    params:
        vsc_min_fraction_indels = "0.12",
        pileup_image_width = 199,
        shard = lambda wc: wc.shard,
        reads = lambda wc: ','.join([f"{output_dir}/{wc.condition}/whatshap_intermediate/{movie}.deepvariant.haplotagged.bam" for movie in ubam_dict[wc.condition].keys()])
    resources: 
        extra = '--constraint=avx512'
    shell:
        f"""
        (/opt/deepvariant/bin/make_examples \
            --norealign_reads \
            --vsc_min_fraction_indels {{params.vsc_min_fraction_indels}} \
            --pileup_image_width {{params.pileup_image_width}} \
            --alt_aligned_pileup=diff_channels \
            --add_hp_channel \
            --sort_by_haplotypes \
            --parse_sam_aux_fields \
            --mode calling \
            --ref {{input.ref}} \
            --reads {{params.reads}} \
            --examples {output_dir}/{{wildcards.condition}}/deepvariant/examples/examples.tfrecord@{n_shards}.gz \
            --gvcf {output_dir}/{{wildcards.condition}}/deepvariant/examples/gvcf.tfrecord@{n_shards}.gz \
            --task {{wildcards.shard}}) > {{log}} 2>&1
        """


rule deepvariant_call_variants_gpu_round2:
    input: lambda wc: expand(f"{output_dir}/{wc.condition}/deepvariant/examples/examples.tfrecord-{{shard}}-of-{n_shards:05}.gz", shard=shards)
    output: f"{output_dir}/{{condition}}/deepvariant/call_variants_output.tfrecord.gz"
    log: f"{output_dir}/{{condition}}/logs/deepvariant_call_variants_gpu_round2.log"
    container: f"docker://google/deepvariant:{deepvariant_version}"
    params: model = "/opt/models/pacbio/model.ckpt"
    threads: 8
    resources:
        partition = 'ml',
        extra = '--gpus=1'
    shell:
        f"""
        (echo "CUDA_VISIBLE_DEVICES=" $CUDA_VISIBLE_DEVICES; \
        /opt/deepvariant/bin/call_variants \
            --outfile {{output}} \
            --examples {output_dir}/{{wildcards.condition}}/deepvariant/examples/examples.tfrecord@{n_shards}.gz \
            --checkpoint {{params.model}}) > {{log}} 2>&1
        """


rule deepvariant_postprocess_variants_round2:
    input:
        tfrecord = f"{output_dir}/{{condition}}/deepvariant/call_variants_output.tfrecord.gz",
        nonvariant_site_tfrecord = lambda wc: expand(f"{output_dir}/{wc.condition}/deepvariant/examples/gvcf.tfrecord-{{shard:05}}-of-{n_shards:05}.gz",
                                          shard=range(n_shards)),
        ref = reference
    output:
        vcf = f"{output_dir}/{{condition}}/deepvariant/deepvariant.vcf.gz",
        vcf_index = f"{output_dir}/{{condition}}/deepvariant/deepvariant.vcf.gz.tbi",
        gvcf = f"{output_dir}/{{condition}}/deepvariant/deepvariant.g.vcf.gz",
        gvcf_index = f"{output_dir}/{{condition}}/deepvariant/deepvariant.g.vcf.gz.tbi",
        report = f"{output_dir}/{{condition}}/deepvariant/deepvariant.visual_report.html"
    log: f"{output_dir}/{{condition}}/logs/deepvariant_postprocess_variants_round2.log"
    container: f"docker://google/deepvariant:{deepvariant_version}"
    threads: 4
    resources: 
        extra = '--constraint=avx512'
    shell:
        f"""
        (/opt/deepvariant/bin/postprocess_variants \
            --ref {{input.ref}} \
            --infile {{input.tfrecord}} \
            --outfile {{output.vcf}} \
            --nonvariant_site_tfrecord_path {output_dir}/{{wildcards.condition}}/deepvariant/examples/gvcf.tfrecord@{n_shards}.gz \
            --gvcf_outfile {{output.gvcf}}) > {{log}} 2>&1
        """


rule happy_benchmark_deepvariant:
    input:
        ref = reference,
        query_vcf = rules.deepvariant_postprocess_variants_round2.output.vcf,
        query_tbi = rules.deepvariant_postprocess_variants_round2.output.vcf_index,
        bench_vcf = lambda wc: happy_vcf_dict[wc.condition][wc.version],
        strats = stratifications,
        sdf = sdf
    output:
        f"{output_dir}/{{condition}}/happy/{{version}}.extended.csv",
        f"{output_dir}/{{condition}}/happy/{{version}}.metrics.json.gz",
        f"{output_dir}/{{condition}}/happy/{{version}}.roc.all.csv.gz",
        f"{output_dir}/{{condition}}/happy/{{version}}.roc.Locations.INDEL.csv.gz",
        f"{output_dir}/{{condition}}/happy/{{version}}.roc.Locations.INDEL.PASS.csv.gz",
        f"{output_dir}/{{condition}}/happy/{{version}}.roc.Locations.SNP.csv.gz",
        f"{output_dir}/{{condition}}/happy/{{version}}.roc.Locations.SNP.PASS.csv.gz",
        f"{output_dir}/{{condition}}/happy/{{version}}.runinfo.json",
        f"{output_dir}/{{condition}}/happy/{{version}}.summary.csv",
        f"{output_dir}/{{condition}}/happy/{{version}}.vcf.gz",
        f"{output_dir}/{{condition}}/happy/{{version}}.vcf.gz.tbi"
    log: f"{output_dir}/{{condition}}/logs/happy_benchmark_deepvariant.{{version}}.log"
    container: "docker://pkrusche/hap.py:latest"
    params:
        prefix = f"{output_dir}/{{condition}}/happy/{{version}}",
        bed = lambda wc: "-f "+happy_bed_dict[wc.condition][wc.version] if wc.version in happy_bed_dict[wc.condition].keys() else ""
    threads: 12
    shell:
        """
        (/opt/hap.py/bin/hap.py \
            --threads {threads} \
            -r {input.ref} {params.bed} \
            -o {params.prefix} \
            --engine=vcfeval --engine-vcfeval-template {input.sdf} \
            --stratification {input.strats} \
            {input.bench_vcf} {input.query_vcf}) > {log} 2>&1
        """


rule deepvariant_bcftools_stats:
    input: rules.deepvariant_postprocess_variants_round2.output.vcf
    output: f"{output_dir}/{{condition}}/deepvariant/deepvariant.vcf.stats.txt"
    log: f"{output_dir}/{{condition}}/logs/deepvariant_bcftools_stats.log"
    params: f"--fasta-ref {reference} --apply-filters PASS -s {{condition}}"
    threads: 4
    conda: "envs/bcftools.yaml"
    shell: "(bcftools stats --threads 3 {params} {input} > {output}) > {log} 2>&1"


rule split_deepvariant_vcf_round2:
    input: rules.deepvariant_postprocess_variants_round2.output.vcf
    output: f"{output_dir}/{{condition}}/whatshap/{{chromosome}}.deepvariant.vcf"
    log: f"{output_dir}/{{condition}}/logs/split_deepvariant_vcf_round2.{{chromosome}}.log"
    params: region = lambda wc: wc.chromosome, extra = '-h'
    conda: "envs/htslib.yaml"
    shell: "tabix {params.extra} {input} {params.region} > {output} 2> {log}"


rule whatshap_phase_round2:
    input:
        ref = reference,
        vcf = f"{output_dir}/{{condition}}/whatshap/{{chromosome}}.deepvariant.vcf.gz",
        tbi = f"{output_dir}/{{condition}}/whatshap/{{chromosome}}.deepvariant.vcf.gz.tbi",
        phaseinput = lambda wc: [f"{output_dir}/{wc.condition}/aligned/{movie}.{reference_name}.bam" for movie in ubam_dict[wc.condition].keys()],
        phaseinputindex = lambda wc: [f"{output_dir}/{wc.condition}/aligned/{movie}.{reference_name}.bam.bai" for movie in ubam_dict[wc.condition].keys()]
    output: f"{output_dir}/{{condition}}/whatshap/{{chromosome}}.deepvariant.phased.vcf.gz"
    log: f"{output_dir}/{{condition}}/logs/whatshap_phase_round2.{{chromosome}}.log"
    params: chromosome = lambda wc: wc.chromosome, extra = "--indels"
    conda: "envs/whatshap.yaml"
    shell:
        """
        (whatshap phase {params.extra} \
            --chromosome {wildcards.chromosome} \
            --output {output} \
            --reference {input.ref} \
            {input.vcf} \
            {input.phaseinput}) > {log} 2>&1
        """


rule whatshap_bcftools_concat_round2:
    input:
        calls = lambda wc: expand(f"{output_dir}/{wc.condition}/whatshap/{{chromosome}}.deepvariant.phased.vcf.gz", chromosome=all_chroms),
        indices = lambda wc: expand(f"{output_dir}/{wc.condition}/whatshap/{{chromosome}}.deepvariant.phased.vcf.gz.tbi", chromosome=all_chroms)
    output: f"{output_dir}/{{condition}}/whatshap/deepvariant.phased.vcf.gz"
    log: f"{output_dir}/{{condition}}/logs/whatshap_bcftools_concat_round2.log"
    params: "-a -Oz"
    conda: "envs/bcftools.yaml"
    shell: "bcftools concat {params} -o {output} {input.calls} > {log} 2>&1"


rule whatshap_stats:
    input:
        vcf = rules.whatshap_bcftools_concat_round2.output,
        tbi = f"{output_dir}/{{condition}}/whatshap/deepvariant.phased.vcf.gz.tbi",
        chr_lengths = chromosome_lengths
    output:
        gtf = f"{output_dir}/{{condition}}/whatshap/deepvariant.phased.gtf",
        tsv = f"{output_dir}/{{condition}}/whatshap/deepvariant.phased.tsv",
        blocklist = f"{output_dir}/{{condition}}/whatshap/deepvariant.phased.blocklist"
    log: f"{output_dir}/{{condition}}/logs/whatshap_stats.log"
    conda: "envs/whatshap.yaml"
    shell:
        """
        (whatshap stats \
            --gtf {output.gtf} \
            --tsv {output.tsv} \
            --block-list {output.blocklist} \
            --chr-lengths {input.chr_lengths} \
            {input.vcf}) > {log} 2>&1
        """


rule whatshap_haplotag_round2:
    input:
        ref = reference,
        vcf = rules.whatshap_bcftools_concat_round2.output,
        tbi = f"{output_dir}/{{condition}}/whatshap/deepvariant.phased.vcf.gz.tbi",
        bam = rules.pbmm2_align.output.bam
    output: f"{output_dir}/{{condition}}/whatshap/{{movie}}.deepvariant.haplotagged.bam"
    log: f"{output_dir}/{{condition}}/logs/whatshap_haplotag_round2.{{movie}}.log"
    params: "--tag-supplementary"
    conda: "envs/whatshap.yaml"
    shell:
        """
        (whatshap haplotag {params} \
            --output {output} \
            --reference {input.ref} \
            {input.vcf} {input.bam}) > {log} 2>&1
        """


rule merge_haplotagged_bams:
    input: lambda wc: expand(f"{output_dir}/{wc.condition}/whatshap/{{movie}}.deepvariant.haplotagged.bam", movie=ubam_dict[wc.condition].keys())
    output: f"{output_dir}/{{condition}}/whatshap/deepvariant.haplotagged.bam"
    log: f"{output_dir}/{{condition}}/logs/merge_haplotagged_bams.log"
    threads: 8
    conda: "envs/samtools.yaml"
    shell: "(samtools merge -@ 7 {output} {input}) > {log} 2>&1"


rule mosdepth:
    input:
        bam = rules.merge_haplotagged_bams.output,
        bai = f"{output_dir}/{{condition}}/whatshap/deepvariant.haplotagged.bam.bai"
    output:
        summary = f"{output_dir}/{{condition}}/mosdepth/coverage.mosdepth.summary.txt",
        regions = f"{output_dir}/{{condition}}/mosdepth/coverage.regions.bed.gz",
        thresholds = f"{output_dir}/{{condition}}/mosdepth/coverage.thresholds.bed.gz"
    log: f"{output_dir}/{{condition}}/logs/mosdepth.log"
    params:
        by = "500",
        prefix = f"{output_dir}/{{condition}}/mosdepth/coverage",
        thresholds = "1,2,3,4,5,6,7,8,9,10",
        extra = "--no-per-base --use-median"
    threads: 4
    conda: "envs/mosdepth.yaml"
    shell:
        """
        (mosdepth \
            --threads {threads} --by {params.by} --thresholds {params.thresholds} \
            {params.extra} {params.prefix} {input.bam}) > {log} 2>&1
        """


rule calculate_m2_ratio:
    input: rules.mosdepth.output.summary
    output: f"{output_dir}/{{condition}}/mosdepth/M2_ratio.txt"
    log: f"{output_dir}/{{condition}}/logs/calculate_M2_ratio.log"
    conda: "envs/pandas.yaml"
    shell: f"(python3 {workflow_dir}/scripts/calculate_M2_ratio.py {{input}} > {{output}}) > {{log}} 2>&1"


rule calculate_gc_coverage:
    input:
        mosdepth_regions = rules.mosdepth.output.regions,
        ref = reference,
    output: f"{output_dir}/{{condition}}/mosdepth/gc_coverage.summary.txt"
    log: f"{output_dir}/{{condition}}/logs/calculate_gc_coverage.log"
    conda: "envs/gc_coverage.yaml"
    shell:
        """
        (bedtools nuc -fi {input.ref} -bed {input.mosdepth_regions} \
            | awk '($11==0) {{ print 0.05*int($6/0.05) "\t" $4; }}' \
            | sort -k1,1g \
            | datamash -g1 q1 2 median 2 q3 2 count 2 \
            | tr '\t' ',' \
            | awk 'BEGIN {{ print "#gc_perc,q1,median,q3,count"; }} {{ print $0; }}' > {output}) > {log} 2>&1
        """


rule calculate_coverage_thresholds:
    input: rules.mosdepth.output.thresholds
    output: f"{output_dir}/{{condition}}/mosdepth/coverage.thresholds.summary.txt"
    log: f"{output_dir}/{{condition}}/logs/calculate_coverage_thresholds.log"
    conda: "envs/gc_coverage.yaml"
    shell: 
        """
        (zcat {input} | awk 'NR==1' | cut -f 5-14 > {output};
            zcat {input} \
                | datamash -H mean 5-14 \
                | awk -v c=500 'NR!=1 {{ for (i = 1; i <= NF; ++i) $i /= c; print }}' >> {output}) > {log} 2>&1
        """


rule whatshap_compare:
    input:
        bench = lambda wc: whatshap_vcf_dict[wc.condition][wc.version],
        vcf = rules.whatshap_bcftools_concat_round2.output,
    output: f"{output_dir}/{{condition}}/whatshap/{{version}}.phase.eval.tsv"
    log: f"{output_dir}/{{condition}}/logs/whatshap_compare.{{version}}.log"
    conda: "envs/whatshap.yaml"
    shell: "(whatshap compare --names benchmark,whatshap --ignore-sample-name --tsv-pairwise {output} {input.bench} {input.vcf}) > {log} 2>&1"
