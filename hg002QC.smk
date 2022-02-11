# takes dict of dicts to permit multiple movies per condition in format ubam_dict[condition][movie]
ubam_dict = { 'condition1' : {  'm64012_190920_173625.ccs': '/pbi/dept/appslab/projects/old/2020/wr_pbRUGD_testing/test_datasets/HG002/m64012_190920_173625.ccs.bam',
                                'm64012_190921_234837.ccs': '/pbi/dept/appslab/projects/old/2020/wr_pbRUGD_testing/test_datasets/HG002/m64012_190921_234837.ccs.bam'}}

reference = '/pbi/dept/appslab/projects/old/2021/TAG-4941/reference/human_GRCh38_no_alt_analysis_set.fasta'
reference_index = '/pbi/dept/appslab/projects/old/2021/TAG-4941/reference/human_GRCh38_no_alt_analysis_set.fasta.fai'
hg002_sr_yak = '/pbi/dept/appslab/projects/old/2021/TAG-4941/resources/yak/HG002.sr.yak'
tandem_repeats_bed = '/pbi/dept/appslab/projects/old/2021/TAG-4941/reference/human_GRCh38_no_alt_analysis_set.trf.bed'
truvari_vcf = '/pbi/dept/appslab/projects/old/2021/TAG-4941/resources/truvari/GRCh38.HG002.latest.sv.vcf.gz'
truvari_bed = '/pbi/dept/appslab/projects/old/2021/TAG-4941/resources/truvari/GRCh38.HG002.latest.sv.bed'
all_chroms = ['chr'+ s for s in list(map(str,*[range(1,23)]))] + ['chrX','chrY','chrM']
DEEPVARIANT_VERSION = '1.2.0-gpu'  # GPU
N_SHARDS = 256

shards = [f"{x:05}" for x in range(N_SHARDS)]

localrules: all, calculate_m2_ratio, calculate_gc_coverage, calculate_coverage_thresholds,bcftools_concat_pbsv_vcf
localrules: split_deepvariant_vcf_round1, split_deepvariant_vcf_round2
localrules: whatshap_bcftools_concat_round1, whatshap_bcftools_concat_round2

targets = []
# targets.extend([f"conditions/{condition}/{filename}"
#                     for condition in ubam_dict.keys()
#                     for filename in ["smrtcell_stats/all_movies.read_length_and_quality.tsv",
#                                     "mosdepth/mosdepth.M2_ratio.txt",
#                                     "mosdepth/gc_coverage.summary.txt",
#                                     "mosdepth/coverage.thresholds.summary.txt",
#                                     "hifiasm/asm.p_ctg.fasta.stats.txt",
#                                     "hifiasm/asm.a_ctg.fasta.stats.txt",
#                                     "hifiasm/asm.p_ctg.qv.txt",
#                                     "hifiasm/asm.a_ctg.qv.txt",
#                                     "truvari/truvari.summary.txt"]])

targets.extend([f"conditions/{condition}/{filename}"
                    for condition in ubam_dict.keys()
                    for filename in ["smrtcell_stats/all_movies.read_length_and_quality.tsv",
                                    "hifiasm/asm.p_ctg.fasta.stats.txt",
                                    "hifiasm/asm.a_ctg.fasta.stats.txt",
                                    "hifiasm/asm.p_ctg.qv.txt",
                                    "hifiasm/asm.a_ctg.qv.txt",
                                    "truvari/truvari.summary.txt",
                                    "deepvariant/intermediate/deepvariant.vcf.gz"]])


rule all:
    input: targets

# outputs 3 columns (read name, read length, read quality)
# PLOT - boxplot of read length and quality
# this is per smrtcell/movie but could be combined by condition
rule smrtcell_stats:
    input: lambda wc: ubam_dict[wc.condition][wc.movie]
    output: "conditions/{condition}/smrtcell_stats/{movie}.read_length_and_quality.tsv"
    log: "conditions/{condition}/logs/smrtcell_stats.{movie}.log"
    conda: "envs/smrtcell_stats.yaml"
    shell: "(python3 workflow/scripts/extract_read_length_and_qual.py {input} > {output}) > {log} 2>&1"


rule combine_smrtcell_stats:
    input: lambda wc: expand(f"conditions/{wc.condition}/smrtcell_stats/{{movie}}.read_length_and_quality.tsv", movie=ubam_dict[wc.condition].keys())
    output: "conditions/{condition}/smrtcell_stats/all_movies.read_length_and_quality.tsv"
    log: "conditions/{condition}/logs/combine_smrtcell_stats.log"
    shell: "(cat {input} > {output}) > {log} 2>&1"    


rule pbmm2_align:
    input:
        ref = reference,
        ref_index = reference_index,
        query = lambda wc: ubam_dict[wc.condition][wc.movie]
    output:
        bam = "conditions/{condition}/aligned/{movie}.aligned.bam",
        bai = "conditions/{condition}/aligned/{movie}.aligned.bam.bai"
    log: "conditions/{condition}/logs/pbmm2_align.{movie}.log"
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


# PLOT: mean aligned coverage in "coverage.mosdepth.summary.txt" - 4th column of final row, can grep 'total_region'
# rule mosdepth:
#     input:
#         bam = "conditions/{condition}/WHATSHAP_BAM.bam",
#         bai = "conditions/{condition}/WHATSHAP_BAM.bam.bai"
#     output:
#         "conditions/{condition}/mosdepth/coverage.mosdepth.summary.txt",
#         "conditions/{condition}/mosdepth/coverage.regions.bed.gz",
#         "conditions/{condition}/mosdepth/coverage.thresholds.bed.gz"
#     log: "conditions/{condition}/logs/mosdepth.log"
#     params:
#         by = "500",
#         prefix = "conditions/{condition}/mosdepth/coverage",
#         thresholds = "1,2,3,4,5,6,7,8,9,10",
#         extra = "--no-per-base --use-median"
#     threads: 4
#     conda: "envs/mosdepth.yaml"
#     shell:
#         """
#         (mosdepth \
#             --threads {threads} --by {params.by} --thresholds {params.thresholds} \
#             {params.extra} {params.prefix} {input.bam}) > {log} 2>&1
#         """


# # outputs single value: ratio of chr2 coverage to chrM coverage
# # PLOT: bar chart of m2 ratio
# rule calculate_m2_ratio:
#     input: "conditions/{condition}/mosdepth/coverage.mosdepth.summary.txt"
#     output: "conditions/{condition}/mosdepth/mosdepth.M2_ratio.txt"
#     log: "conditions/{condition}/logs/calculate_M2_ratio.log"
#     conda: "envs/pandas.yaml"
#     shell: "(python3 workflow/scripts/calculate_M2_ratio.py {input} > {output}) > {log} 2>&1"


# # outputs 5 columns: gc percentage bin, q1 , median , q3 , count
# # q1, median, q3 columns are statistics for coverage at different gc percentages (e.g. median cover at 30% GC)
# # "count" refers to # of 500 bp windows that fall in that bin
# # PLOT: can pick a couple of key GC coverage bins and make box plots out of them
# rule calculate_gc_coverage:
#     input:
#         mosdepth_regions = "conditions/{condition}/mosdepth/coverage.regions.bed.gz",
#         ref = reference,
#     output: "conditions/{condition}/mosdepth/gc_coverage.summary.txt"
#     log: "conditions/{condition}/logs/calculate_gc_coverage.log"
#     conda: "envs/gc_coverage.yaml"
#     shell:
#         """
#         (bedtools nuc -fi {input.ref} -bed {input.mosdepth_regions} \
#             | awk '($11==0) {{ print 0.05*int($6/0.05) "\t" $4; }}' \
#             | sort -k1,1g \
#             | datamash -g1 q1 2 median 2 q3 2 count 2 \
#             | tr '\t' ',' \
#             | awk 'BEGIN {{ print "#gc_perc,q1,median,q3,count"; }} {{ print $0; }}' > {output}) > {log} 2>&1
#         """

# # outputs 10 columns corresponding to % of genome sequenced to minimum coverage depths (1X - 10X)
# # PLOT: maybe a line chart comparing the different coverage thresholds among conditions
# rule calculate_coverage_thresholds:
#     input: "conditions/{condition}/mosdepth/coverage.thresholds.bed.gz"
#     output: "conditions/{condition}/mosdepth/coverage.thresholds.summary.txt"
#     log: "conditions/{condition}/logs/calculate_coverage_thresholds.log"
#     conda: "envs/gc_coverage.yaml"
#     shell: 
#         """
#         (zcat {input} | awk 'NR==1' | cut -f 5-14 > {output};
#             zcat {input} \
#                 | datamash -H mean 5-14 \
#                 | awk -v c=500 'NR!=1 {{ for (i = 1; i <= NF; ++i) $i /= c; print }}' >> {output}) > {log} 2>&1
#         """


rule samtools_fasta:
    input: lambda wc: ubam_dict[wc.condition][wc.movie]
    output: "conditions/{condition}/fasta/{movie}.fasta"
    log: "conditions/{condition}/logs/samtools_fasta.{movie}.log"
    threads: 4
    conda: "envs/samtools.yaml"
    shell: "(samtools fasta -@ 3 {input} > {output}) > {log} 2>&1"


rule hifiasm_assemble:
    input: lambda wc: expand(f"conditions/{wc.condition}/fasta/{{movie}}.fasta", movie=ubam_dict[wc.condition].keys())
    output:
        "conditions/{condition}/hifiasm/asm.p_ctg.gfa",
        "conditions/{condition}/hifiasm/asm.a_ctg.gfa",
    log: "conditions/{condition}/logs/hifiasm.log"
    conda: "envs/hifiasm.yaml"
    params: prefix = "conditions/{condition}/hifiasm/asm"
    threads: 48
    shell: "(hifiasm -o {params.prefix} --primary -t {threads} {input}) > {log} 2>&1"


rule gfa2fa:
    input: "conditions/{condition}/hifiasm/asm.{infix}.gfa"
    output: "conditions/{condition}/hifiasm/asm.{infix}.fasta"
    log: "conditions/{condition}/logs/gfa2fa.{infix}.log"
    conda: "envs/gfatools.yaml"
    shell: "(gfatools gfa2fa {input} > {output}) 2> {log}"


rule bgzip_fasta:
    input: "conditions/{condition}/hifiasm/asm.{infix}.fasta"
    output: "conditions/{condition}/hifiasm/asm.{infix}.fasta.gz"
    log: "conditions/{condition}/logs/bgzip_fasta.{infix}.log"
    threads: 4
    conda: "envs/htslib.yaml"
    shell: "(bgzip --threads {threads} {input}) > {log} 2>&1"


# PLOT fasta stats SZ and and all NL lines as cumulative line plot, aun (heng li argued it's better than NG(X))
rule asm_stats:
    input: 
        fasta = "conditions/{condition}/hifiasm/asm.{infix}.fasta.gz",
        ref_index = reference_index
    output: "conditions/{condition}/hifiasm/asm.{infix}.fasta.stats.txt"
    log: "conditions/{condition}/logs/asm_stats.{infix}.log"
    conda: "envs/k8.yaml"
    shell: "(k8 workflow/scripts/calN50.js -f {input.ref_index} {input.fasta} > {output}) > {log} 2>&1"

# PLOT adjusted_quality_value `asm.a_ctg.qv.txt` for each haplotype
rule yak_qv:
    input:
        asm = "conditions/{condition}/hifiasm/asm.{infix}.fasta.gz",
        sr = hg002_sr_yak
    output: "conditions/{condition}/hifiasm/asm.{infix}.qv.txt"
    log: "conditions/{condition}/logs/yak_qv.{infix}.log"
    threads: 24
    conda: "envs/yak.yaml"
    shell: "(yak qv -t{threads} -p -K3.2g -l100k {input.sr} <(zcat {input.asm}) > {output}) > {log} 2>&1"

localrules: bgzip_vcf, tabix_vcf


rule bgzip_vcf:
    input: "conditions/{condition}/{prefix}.vcf"
    output: "conditions/{condition}/{prefix}.vcf.gz"
    log: "conditions/{condition}/logs/{prefix}.bgzip_vcf.log"
    threads: 2
    conda: "envs/htslib.yaml"
    shell: "(bgzip --threads {threads} {input}) > {log} 2>&1"


rule tabix_vcf:
    input: "conditions/{condition}/{prefix}.vcf.gz"
    output: "conditions/{condition}/{prefix}.vcf.gz.tbi"
    log: "conditions/{condition}/logs/{prefix}.tabix_vcf.log"
    params: "-p vcf"
    conda: "envs/htslib.yaml"
    shell: "tabix {params} {input} > {log} 2>&1"


# rule samtools_index_bam:
#     input: "conditions/{condition}/{prefix}.bam"
#     output: "conditions/{condition}/{prefix}.bam.bai"
#     log: "conditions/{condition}/logs/{prefix}.samtools_index_bam.log"
#     threads: 4
#     conda: "envs/samtools.yaml"
#     shell: "(samtools index -@ 3 {input}) > {log} 2>&1"

rule pbsv_discover:
    input:
        bam = "conditions/{condition}/aligned/{movie}.aligned.bam",
        bai = "conditions/{condition}/aligned/{movie}.aligned.bam.bai",
        tr_bed = tandem_repeats_bed
    output: "conditions/{condition}/pbsv/intermediate/{movie}.{region}.svsig.gz"
    log: "conditions/{condition}/logs/pbsv_discover.{movie}.{region}.log"
    params:
        extra = "--hifi",
        region = lambda wc: wc.region,
        loglevel = "INFO"
    conda: "envs/pbsv.yaml"
    shell:
        """
        (pbsv discover {params.extra} \
            --log-level {params.loglevel} \
            --region {wildcards.region} \
            --tandem-repeats {input.tr_bed} \
            {input.bam} {output}) > {log} 2>&1
        """


rule pbsv_call:
    input:
        svsigs = lambda wc: expand([f"conditions/{wc.condition}/pbsv/intermediate/{{movie}}.{wc.region}.svsig.gz"], movie=ubam_dict[wc.condition].keys()),
        ref = reference
    output: "conditions/{condition}/pbsv/intermediate/{region}.pbsv.vcf"
    log: "conditions/{condition}/logs/pbsv_call.{region}.log"
    params:
        region = lambda wc: wc.region,
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
        calls = lambda wc: expand(f"conditions/{wc.condition}/pbsv/intermediate/{{region}}.pbsv.vcf.gz", region=all_chroms),
        indices = lambda wc: expand(f"conditions/{wc.condition}/pbsv/intermediate/{{region}}.pbsv.vcf.gz.tbi", region=all_chroms)
    output: "conditions/{condition}/all_chroms.pbsv.vcf"
    log: "conditions/{condition}/logs/bcftools_concat_pbsv_vcf.log"
    conda: "envs/bcftools.yaml"
    shell: "(bcftools concat -a -o {output} {input.calls}) > {log} 2>&1"


rule truvari_benchmark:
    input:
        ref = reference,
        query_vcf = "conditions/{condition}/all_chroms.pbsv.vcf.gz",
        query_tbi = "conditions/{condition}/all_chroms.pbsv.vcf.gz.tbi",
        bench_vcf = truvari_vcf,
        bench_bed = truvari_bed
    output: 
        "conditions/{condition}/truvari/truvari.tp-call.vcf",
        "conditions/{condition}/truvari/truvari.tp-base.vcf",
        "conditions/{condition}/truvari/truvari.fn.vcf",
        "conditions/{condition}/truvari/truvari.fp.vcf",
        "conditions/{condition}/truvari/truvari.base-filter.vcf",
        "conditions/{condition}/truvari/truvari.call-filter.vcf",
        "conditions/{condition}/truvari/truvari.summary.txt",
        "conditions/{condition}/truvari/truvari.log.txt",
        "conditions/{condition}/truvari/truvari.giab_report.txt"
    log: "conditions/{condition}/logs/truvari_benchmark.log"
    params: prefix = "conditions/{condition}/truvari/truvari"
    shell:
        """
        (truvari \
            -f {input.ref} \
            -b {input.bench_vcf} \
            --includebed {input.bench_bed} \
            -o {params.prefix} \
            --passonly --giabreport \
            -r 1000 -p 0.00 \
            -c {input.query_vcf}) > {log} 2>&1
        """



rule deepvariant_make_examples_round1:
    input:
        bams = lambda wc: [f"conditions/{wc.condition}/aligned/{movie}.aligned.bam" for movie in ubam_dict[wc.condition].keys()],
        bais = lambda wc: [f"conditions/{wc.condition}/aligned/{movie}.aligned.bam.bai" for movie in ubam_dict[wc.condition].keys()],
        ref = reference
    output:
        tfrecord = f"conditions/{{condition}}/deepvariant/intermediate/examples/examples.tfrecord-{{shard}}-of-{N_SHARDS:05}.gz"
    log: f"conditions/{{condition}}/logs/deepvariant_make_examples_round1.{{shard}}-of-{N_SHARDS:05}.log"
    container: f"docker://google/deepvariant:{DEEPVARIANT_VERSION}"
    params:
        vsc_min_fraction_indels = "0.12",
        pileup_image_width = 199,
        shard = lambda wc: wc.shard,
        reads = lambda wc: ','.join(["conditions/{wc.condition}/aligned/{movie}.aligned.bam" for movie in ubam_dict[wc.condition].keys()])
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
            --examples conditions/{{wildcards.condition}}/deepvariant/intermediate/examples/examples.tfrecord@{N_SHARDS}.gz \
            --task {{wildcards.shard}}) > {{log}} 2>&1
        """


rule deepvariant_call_variants_gpu_round1:
    input: lambda wc: expand(f"conditions/{wc.condition}/deepvariant/intermediate/examples/examples.tfrecord-{{shard}}-of-{N_SHARDS:05}.gz", shard=shards)
    output: "conditions/{condition}/deepvariant/intermediate/call_variants_output.tfrecord.gz"
    log: "conditions/{condition}/logs/deepvariants_call_variants_gpu_round1.log"
    container: f"docker://google/deepvariant:{DEEPVARIANT_VERSION}"
    params: model = "/opt/models/pacbio/model.ckpt"
    threads: 8
    shell:
        f"""
        (/opt/deepvariant/bin/call_variants \
            --outfile {{output}} \
            --examples conditions/{{wildcards.condition}}/deepvariant/intermediate/examples/examples.tfrecord@{N_SHARDS}.gz \
            --checkpoint {{params.model}}) > {{log}} 2>&1
        """


rule deepvariant_postprocess_variants_round1:
    input:
        tfrecord = "conditions/{condition}/deepvariant/intermediate/call_variants_output.tfrecord.gz",
        ref = reference
    output:
        vcf = "conditions/{condition}/deepvariant/intermediate/deepvariant.vcf.gz",
        vcf_index = "conditions/{condition}/deepvariant/intermediate/deepvariant.vcf.gz.tbi",
        report = "conditions/{condition}/deepvariant/intermediate/deepvariant.visual_report.html"
    log: "conditions/{condition}/logs/deepvariant/intermediate/postprocess_variants/deepvariant_postprocess_variants_round1.log"
    container: f"docker://google/deepvariant:{DEEPVARIANT_VERSION}"
    threads: 4
    shell:
        """
        (/opt/deepvariant/bin/postprocess_variants \
            --ref {input.ref} \
            --infile {input.tfrecord} \
            --outfile {output.vcf}) > {log} 2>&1
        """


rule split_deepvariant_vcf_round1:
    input: "conditions/{condition}/deepvariant/intermediate/deepvariant.vcf.gz",
    output: "conditions/{condition}/whatshap/intermediate/{chromosome}.deepvariant.vcf"
    log: "conditions/{condition}/logs/tabix/query/{chromosome}.deepvariant/intermediate.vcf.log"
    params: extra = '-h'
    conda: "envs/htslib.yaml"
    shell: "tabix {params.extra} {input} {wildcards.chromosome} > {output} 2> {log}"


rule whatshap_phase_round1:
    input:
        ref = reference,
        vcf = "conditions/{condition}/whatshap/intermediate/{chromosome}.deepvariant.vcf.gz",
        tbi = "conditions/{condition}/whatshap/intermediate/{chromosome}.deepvariant.vcf.gz.tbi",
        phaseinput = lambda wc: ["conditions/{wc.condition}/aligned/{movie}.aligned.bam" for movie in ubam_dict[wc.condition].keys()],
        phaseinputindex = lambda wc: ["conditions/{wc.condition}/aligned/{movie}.aligned.bam.bai" for movie in ubam_dict[wc.condition].keys()],
    output: "conditions/{condition}/whatshap/intermediate/{chromosome}.deepvariant.phased.vcf.gz"
    log: "conditions/{condition}/logs/whatshap_intermediate.{chromosome}.log"
    params: chromosome = lambda wc: wc.chromosome
    conda: "envs/whatshap.yaml"
    shell:
        """
        (whatshap phase \
            --chromosome {wildcards.chromosome} \
            --output {output} \
            --reference {input.reference} \
            {input.vcf} {input.phaseinput}) > {log} 2>&1
        """


# rule whatshap_bcftools_concat_round1:
#     input:
#         calls = expand(f"conditions/{condition}/whatshap_intermediate/{sample}.{ref}.regions/{sample}.{ref}.{{region}}.deepvariant.phased.vcf.gz", region=all_chroms),
#         indices = expand(f"conditions/{condition}/whatshap_intermediate/{sample}.{ref}.regions/{sample}.{ref}.{{region}}.deepvariant.phased.vcf.gz.tbi", region=all_chroms)
#     output: f"conditions/{condition}/whatshap_intermediate/{sample}.{ref}.deepvariant.phased.vcf.gz"
#     log: f"conditions/{condition}/logs/bcftools/concat/{sample}.{ref}.whatshap_intermediate.log"
#     benchmark: f"conditions/{condition}/benchmarks/bcftools/concat/{sample}.{ref}.whatshap_intermediate.tsv"
#     params: "-a -Oz"
#     conda: "envs/bcftools.yaml"
#     message: "Executing {rule}: Concatenating WhatsHap phased VCFs: {input.calls}"
#     shell: "bcftools concat {params} -o {output} {input.calls} > {log} 2>&1"


# rule whatshap_haplotag_round1:
#     input:
#         reference = config['ref']['fasta'],
#         vcf = f"conditions/{condition}/whatshap_intermediate/{sample}.{ref}.deepvariant.phased.vcf.gz",
#         tbi = f"conditions/{condition}/whatshap_intermediate/{sample}.{ref}.deepvariant.phased.vcf.gz.tbi",
#         bam = lambda wildcards: abam_dict[wildcards.movie]
#     output: temp(f"conditions/{condition}/whatshap_intermediate/{sample}.{ref}.{{movie}}.deepvariant.haplotagged.bam")
#     log: f"conditions/{condition}/logs/whatshap/haplotag/{sample}.{ref}.{{movie}}.whatshap_intermediate.log"
#     benchmark: f"conditions/{condition}/benchmarks/whatshap/haplotag/{sample}.{ref}.{{movie}}.whatshap_intermediate.tsv"
#     params: "--tag-supplementary"
#     conda: "envs/whatshap.yaml"
#     message: "Executing {rule}: Haplotagging {input.bam} using phase information from {input.vcf}."
#     shell:
#         """
#         (whatshap haplotag {params} \
#             --output {output} \
#             --reference {input.reference} \
#             {input.vcf} {input.bam}) > {log} 2>&1
#         """











# haplotagged_abams = [f"conditions/{condition}/whatshap_intermediate/{sample}.{ref}.{movie}.deepvariant.haplotagged.bam" for movie in movies]


# rule deepvariant_make_examples_round2:
#     input:
#         bams = haplotagged_abams,
#         bais = [f"{x}.bai" for x in haplotagged_abams],
#         reference = config['ref']['fasta']
#     output:
#         tfrecord = temp(f"conditions/{condition}/deepvariant/examples/examples.tfrecord-{{shard}}-of-{N_SHARDS:05}.gz"),
#         nonvariant_site_tfrecord = temp(f"conditions/{condition}/deepvariant/examples/gvcf.tfrecord-{{shard}}-of-{N_SHARDS:05}.gz")
#     log: f"conditions/{condition}/logs/deepvariant/make_examples/{sample}.{ref}.{{shard}}-of-{N_SHARDS:05}.log"
#     benchmark: f"conditions/{condition}/benchmarks/deepvariant/make_examples/{sample}.{ref}.{{shard}}-of-{N_SHARDS:05}.tsv"
#     container: f"docker://google/deepvariant:{DEEPVARIANT_VERSION}"
#     params:
#         vsc_min_fraction_indels = "0.12",
#         pileup_image_width = 199,
#         shard = lambda wildcards: wildcards.shard,
#         reads = ','.join(haplotagged_abams)
#     message: "Executing {rule}: DeepVariant make_examples {wildcards.shard} for {input.bams}."
#     shell:
#         f"""
#         (/opt/deepvariant/bin/make_examples \
#             --norealign_reads \
#             --vsc_min_fraction_indels {{params.vsc_min_fraction_indels}} \
#             --pileup_image_width {{params.pileup_image_width}} \
#             --alt_aligned_pileup=diff_channels \
#             --add_hp_channel \
#             --sort_by_haplotypes \
#             --parse_sam_aux_fields \
#             --mode calling \
#             --ref {{input.reference}} \
#             --reads {{params.reads}} \
#             --examples conditions/{condition}/deepvariant/examples/examples.tfrecord@{N_SHARDS}.gz \
#             --gvcf conditions/{condition}/deepvariant/examples/gvcf.tfrecord@{N_SHARDS}.gz \
#             --task {{wildcards.shard}}) > {{log}} 2>&1
#         """


# rule deepvariant_call_variants_gpu_round2:
#     input: expand(f"conditions/{condition}/deepvariant/examples/examples.tfrecord-{{shard}}-of-{N_SHARDS:05}.gz", shard=shards)
#     output: temp(f"conditions/{condition}/deepvariant/{sample}.{ref}.call_variants_output.tfrecord.gz")
#     log: f"conditions/{condition}/logs/deepvariant/call_variants/{sample}.{ref}.log"
#     benchmark: f"conditions/{condition}/benchmarks/deepvariant/call_variants/{sample}.{ref}.tsv"
#     container: f"docker://google/deepvariant:{DEEPVARIANT_VERSION}"
#     params: model = "/opt/models/pacbio/model.ckpt"
#     message: "Executing {rule}: DeepVariant call_variants for {input}."
#     threads: 8
#     shell:
#         f"""
#         (echo "CUDA_VISIBLE_DEVICES=" $CUDA_VISIBLE_DEVICES; \
#         /opt/deepvariant/bin/call_variants \
#             --outfile {{output}} \
#             --examples conditions/{condition}/deepvariant/examples/examples.tfrecord@{N_SHARDS}.gz \
#             --checkpoint {{params.model}}) > {{log}} 2>&1
#         """


# rule deepvariant_postprocess_variants_round2:
#     input:
#         tfrecord = f"conditions/{condition}/deepvariant/{sample}.{ref}.call_variants_output.tfrecord.gz",
#         nonvariant_site_tfrecord = expand(f"conditions/{condition}/deepvariant/examples/gvcf.tfrecord-{{shard:05}}-of-{N_SHARDS:05}.gz",
#                                           shard=range(N_SHARDS)),
#         reference = config['ref']['fasta']
#     output:
#         vcf = f"conditions/{condition}/deepvariant/{sample}.{ref}.deepvariant.vcf.gz",
#         vcf_index = f"conditions/{condition}/deepvariant/{sample}.{ref}.deepvariant.vcf.gz.tbi",
#         gvcf = f"conditions/{condition}/deepvariant/{sample}.{ref}.deepvariant.g.vcf.gz",
#         gvcf_index = f"conditions/{condition}/deepvariant/{sample}.{ref}.deepvariant.g.vcf.gz.tbi",
#         report = f"conditions/{condition}/deepvariant/{sample}.{ref}.deepvariant.visual_report.html"
#     log: f"conditions/{condition}/logs/deepvariant/postprocess_variants/{sample}.{ref}.log"
#     benchmark: f"conditions/{condition}/benchmarks/deepvariant/postprocess_variants/{sample}.{ref}.tsv"
#     container: f"docker://google/deepvariant:{DEEPVARIANT_VERSION}"
#     message: "Executing {rule}: DeepVariant postprocess_variants for {input.tfrecord}."
#     threads: 4
#     shell:
#         f"""
#         (/opt/deepvariant/bin/postprocess_variants \
#             --ref {{input.reference}} \
#             --infile {{input.tfrecord}} \
#             --outfile {{output.vcf}} \
#             --nonvariant_site_tfrecord_path conditions/{condition}/deepvariant/examples/gvcf.tfrecord@{N_SHARDS}.gz \
#             --gvcf_outfile {{output.gvcf}}) > {{log}} 2>&1
#         """


# rule deepvariant_bcftools_stats:
#     input: f"conditions/{condition}/deepvariant/{sample}.{ref}.deepvariant.vcf.gz"
#     output: f"conditions/{condition}/deepvariant/{sample}.{ref}.deepvariant.vcf.stats.txt"
#     log: f"conditions/{condition}/logs/bcftools/stats/{sample}.{ref}.deepvariant.vcf.log"
#     benchmark: f"conditions/{condition}/benchmarks/bcftools/stats/{sample}.{ref}.deepvariant.vcf.tsv"
#     params: f"--fasta-ref {config['ref']['fasta']} --apply-filters PASS -s {sample}"
#     threads: 4
#     conda: "envs/bcftools.yaml"
#     message: "Executing {rule}: Calculating VCF statistics for {input}."
#     shell: "(bcftools stats --threads 3 {params} {input} > {output}) > {log} 2>&1"


# rule deepvariant_bcftools_roh:
#     input: f"conditions/{condition}/deepvariant/{sample}.{ref}.deepvariant.vcf.gz"
#     output: f"conditions/{condition}/deepvariant/{sample}.{ref}.deepvariant.roh.bed"
#     log: f"conditions/{condition}/logs/bcftools/stats/{sample}.{ref}.deepvariant.vcf.log"
#     benchmark: f"conditions/{condition}/benchmarks/bcftools/stats/{sample}.{ref}.deepvariant.vcf.tsv"
#     params: default_allele_frequency = 0.4
#     conda: "envs/bcftools.yaml"
#     message: "Executing {rule}: Calculating runs of homozygosity for {input}."
#     shell:
#         """
#         (echo -e "#chr\tstart\tend\tqual" > {output}
#         bcftools roh --AF-dflt {params.default_allele_frequency} {input} \
#         | awk -v OFS='\t' '$1=="RG" {{ print $3, $4, $5, $8 }}' \
#         >> {output}) > {log} 2>&1
#         """


# rule split_deepvariant_vcf_round2:
#     input: f"conditions/{condition}/deepvariant/{sample}.{ref}.deepvariant.vcf.gz",
#     output: temp(f"conditions/{condition}/whatshap/{sample}.{ref}.regions/{sample}.{ref}.{{region}}.deepvariant.vcf")
#     log: f"conditions/{condition}/logs/tabix/query/{sample}.{ref}.{{region}}.deepvariant.vcf.log"
#     benchmark: f"conditions/{condition}/benchmarks/tabix/query/{sample}.{ref}.{{region}}.deepvariant.vcf.tsv"
#     params: region = lambda wildcards: wildcards.region, extra = '-h'
#     conda: "envs/htslib.yaml"
#     message: "Executing {rule}: Extracting {wildcards.region} variants from {input}."
#     shell: "tabix {params.extra} {input} {params.region} > {output} 2> {log}"


# rule whatshap_phase_round2:
#     input:
#         reference = config['ref']['fasta'],
#         vcf = f"conditions/{condition}/whatshap/{sample}.{ref}.regions/{sample}.{ref}.{{chromosome}}.deepvariant.vcf.gz",
#         tbi = f"conditions/{condition}/whatshap/{sample}.{ref}.regions/{sample}.{ref}.{{chromosome}}.deepvariant.vcf.gz.tbi",
#         phaseinput = abams,
#         phaseinputindex = [f"{x}.bai" for x in abams]
#     output: temp(f"conditions/{condition}/whatshap/{sample}.{ref}.regions/{sample}.{ref}.{{chromosome}}.deepvariant.phased.vcf.gz")
#     log: f"conditions/{condition}/logs/whatshap/phase/{sample}.{ref}.{{chromosome}}.log"
#     benchmark: f"conditions/{condition}/benchmarks/whatshap/phase/{sample}.{ref}.{{chromosome}}.tsv"
#     params: chromosome = lambda wildcards: wildcards.chromosome, extra = "--indels"
#     conda: "envs/whatshap.yaml"
#     message: "Executing {rule}: Phasing {input.vcf} using {input.phaseinput} for chromosome {wildcards.chromosome}."
#     shell:
#         """
#         (whatshap phase {params.extra} \
#             --chromosome {wildcards.chromosome} \
#             --output {output} \
#             --reference {input.reference} \
#             {input.vcf} \
#             {input.phaseinput}) > {log} 2>&1
#         """


# rule whatshap_bcftools_concat_round2:
#     input:
#         calls = expand(f"conditions/{condition}/whatshap/{sample}.{ref}.regions/{sample}.{ref}.{{region}}.deepvariant.phased.vcf.gz", region=all_chroms),
#         indices = expand(f"conditions/{condition}/whatshap/{sample}.{ref}.regions/{sample}.{ref}.{{region}}.deepvariant.phased.vcf.gz.tbi", region=all_chroms)
#     output: f"conditions/{condition}/whatshap/{sample}.{ref}.deepvariant.phased.vcf.gz"
#     log: f"conditions/{condition}/logs/bcftools/concat/{sample}.{ref}.whatshap.log"
#     benchmark: f"conditions/{condition}/benchmarks/bcftools/concat/{sample}.{ref}.whatshap.tsv"
#     params: "-a -Oz"
#     conda: "envs/bcftools.yaml"
#     message: "Executing {rule}: Concatenating WhatsHap phased VCFs: {input.calls}"
#     shell: "bcftools concat {params} -o {output} {input.calls} > {log} 2>&1"


# rule whatshap_stats:
#     input:
#         vcf = f"conditions/{condition}/whatshap/{sample}.{ref}.deepvariant.phased.vcf.gz",
#         tbi = f"conditions/{condition}/whatshap/{sample}.{ref}.deepvariant.phased.vcf.gz.tbi",
#         chr_lengths = config['ref']['chr_lengths']
#     output:
#         gtf = f"conditions/{condition}/whatshap/{sample}.{ref}.deepvariant.phased.gtf",
#         tsv = f"conditions/{condition}/whatshap/{sample}.{ref}.deepvariant.phased.tsv",
#         blocklist = f"conditions/{condition}/whatshap/{sample}.{ref}.deepvariant.phased.blocklist"
#     log: f"conditions/{condition}/logs/whatshap/stats/{sample}.{ref}.log"
#     benchmark: f"conditions/{condition}/benchmarks/whatshap/stats/{sample}.{ref}.tsv"
#     conda: "envs/whatshap.yaml"
#     message: "Executing {rule}: Calculating phasing stats for {input.vcf}."
#     shell:
#         """
#         (whatshap stats \
#             --gtf {output.gtf} \
#             --tsv {output.tsv} \
#             --block-list {output.blocklist} \
#             --chr-lengths {input.chr_lengths} \
#             {input.vcf}) > {log} 2>&1
#         """


# rule whatshap_haplotag_round2:
#     input:
#         reference = config['ref']['fasta'],
#         vcf = f"conditions/{condition}/whatshap/{sample}.{ref}.deepvariant.phased.vcf.gz",
#         tbi = f"conditions/{condition}/whatshap/{sample}.{ref}.deepvariant.phased.vcf.gz.tbi",
#         bam = lambda wildcards: abam_dict[wildcards.movie]
#     output: temp(f"conditions/{condition}/whatshap/{sample}.{ref}.{{movie}}.deepvariant.haplotagged.bam")
#     log: f"conditions/{condition}/logs/whatshap/haplotag/{sample}.{ref}.{{movie}}.log"
#     benchmark: f"conditions/{condition}/benchmarks/whatshap/haplotag/{sample}.{ref}.{{movie}}.tsv"
#     params: "--tag-supplementary"
#     conda: "envs/whatshap.yaml"
#     message: "Executing {rule}: Haplotagging {input.bam} using phase information from {input.vcf}."
#     shell:
#         """
#         (whatshap haplotag {params} \
#             --output {output} \
#             --reference {input.reference} \
#             {input.vcf} {input.bam}) > {log} 2>&1
#         """


# rule merge_haplotagged_bams:
#     input: expand(f"conditions/{condition}/whatshap/{sample}.{ref}.{{movie}}.deepvariant.haplotagged.bam", movie=movies)
#     output: f"conditions/{condition}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam"
#     log: f"conditions/{condition}/logs/samtools/merge/{sample}.{ref}.haplotag.log"
#     benchmark: f"conditions/{condition}/benchmarks/samtools/merge/{sample}.{ref}.haplotag.tsv"
#     threads: 8
#     conda: "envs/samtools.yaml"
#     message: "Executing {rule}: Merging {input}."
#     shell: "(samtools merge -@ 7 {output} {input}) > {log} 2>&1"


#TANDEM GENOTYPES

# rule download_tg_list:
#     output: config['ref']['tg_list']
#     log: "logs/download_tg_list.log"
#     params: url = config['ref']['tg_list_url']
#     message: "Executing {rule}: Downloading a list of loci with disease-associated repeats to {output}."
#     shell: "(wget -qO - {params.url} > {output}) > {log} 2>&1"


# rule generate_tg_bed:
#     input:
#         tg_list = config['ref']['tg_list'],
#         fai = config['ref']['index']
#     output: config['ref']['tg_bed']
#     log: "logs/generate_tg_bed.log"
#     conda: "envs/bedtools.yaml"
#     params: slop = 1000
#     message: "Executing {rule}: Adding {params.slop}bp slop to {input.tg_list} to generate {output}."
#     shell: 
#         """
#         (grep -v '^#' {input.tg_list} | sort -k1,1V -k2,2g \
#         | bedtools slop -b {params.slop} -g {input.fai} -i - \
#         > {output}) > {log} 2>&1
#         """


# rule generate_last_index:
#     input: config['ref']['fasta']
#     output:
#         [f"{config['ref']['last_index']}.{suffix}"
#          for suffix in ['bck', 'des', 'prj', 'sds', 'ssp', 'suf', 'tis']]
#     log: "logs/generate_last_index.log"
#     conda: "envs/last.yaml"
#     params: "-uRY32 -R01"
#     threads: 24
#     message: "Executing {rule}: Generating last index of {input} using params: {params}"
#     shell: f"(lastdb -P{{threads}} {{params}} {config['ref']['last_index']} {{input}}) > {{log}} 2>&1"


# rule last_align:
#     input:
#         [f"{config['ref']['last_index']}.{suffix}"
#          for suffix in ['bck', 'des', 'prj', 'sds', 'ssp', 'suf', 'tis']],
#         bam = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam",
#         bai = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam.bai",
#         bed = config['ref']['tg_bed'],
#         score_matrix = config['score_matrix']
#     output: temp(f"samples/{sample}/tandem-genotypes/{sample}.maf.gz")
#     log: f"samples/{sample}/logs/last/align/{sample}.log"
#     benchmark: f"samples/{sample}/benchmarks/last/align/{sample}.tsv"
#     conda: "envs/last.yaml"
#     params: 
#         last_index = config['ref']['last_index'],
#         extra = "-C2"
#     threads: 24
#     message: "Executing {rule}: Aligning {input.bed} regions of {input.bam} to {params.last_index} using lastal with {input.score_matrix} score matrix."
#     shell:
#         """
#         (samtools view -@3 -bL {input.bed} {input.bam} | samtools fasta \
#          | lastal -P20 -p {input.score_matrix} {params.extra} {params.last_index} - \
#          | last-split | bgzip > {output}) > {log} 2>&1
#         """


# rule tandem_genotypes:
#     input:
#         maf = f"samples/{sample}/tandem-genotypes/{sample}.maf.gz",
#         repeats = config['ref']['tg_list']
#     output: f"samples/{sample}/tandem-genotypes/{sample}.tandem-genotypes.txt"
#     log: f"samples/{sample}/logs/tandem-genotypes/{sample}.log"
#     benchmark: f"samples/{sample}/benchmarks/tandem-genotypes/{sample}.tsv"
#     conda: "envs/tandem-genotypes.yaml"
#     message: "Executing {rule}: Genotyping tandem repeats from {input.repeats} regions in {input.maf}."
#     shell: "(tandem-genotypes {input.repeats} {input.maf} > {output}) > {log} 2>&1"


# rule tandem_genotypes_absolute_count:
#     input: f"samples/{sample}/tandem-genotypes/{sample}.tandem-genotypes.txt"
#     output: f"samples/{sample}/tandem-genotypes/{sample}.tandem-genotypes.absolute.txt"
#     log: f"samples/{sample}/logs/tandem-genotypes/{sample}.absolute.log"
#     benchmark: f"samples/{sample}/benchmarks/tandem-genotypes/{sample}.absolute.tsv"
#     message: "Executing {rule}: Adjusting repeat count with reference counts for {input}."
#     shell:
#         """
#         (awk -v OFS='\t' \
#             '$0 ~ /^#/ {{print $0 " modified by adding reference repeat count"}}
#             $0 !~ /^#/ {{
#                 ref_count=int(($3-$2)/length($4));
#                 num_fwd=split($7, fwd, ",");
#                 num_rev=split($8, rev, ",");
#                 new_fwd=result=fwd[1] + ref_count;
#                 for (i=2; i<=num_fwd; i++)
#                     new_fwd = new_fwd "," fwd[i] + ref_count;
#                 new_rev=rev[1] + ref_count;
#                 for (i=2; i<=num_rev; i++)
#                     new_rev = new_rev "," rev[i] + ref_count;
#                 print $1, $2, $3, $4, $5, $6, new_fwd, new_rev;
#             }}' {input} > {output} \
#         ) > {log} 2>&1
#         """


# rule tandem_genotypes_plot:
#     input: f"samples/{sample}/tandem-genotypes/{sample}.tandem-genotypes.txt"
#     output: f"samples/{sample}/tandem-genotypes/{sample}.tandem-genotypes.pdf"
#     log: f"samples/{sample}/logs/tandem-genotypes/plot/{sample}.log"
#     benchmark: f"samples/{sample}/benchmarks/tandem-genotypes/plot/{sample}.tsv"
#     conda: "envs/tandem-genotypes.yaml"
#     params: top_N_plots = 100
#     message: "Executing {rule}: Plotting tandem repeats from {input}."
#     shell: "(tandem-genotypes-plot -n {params.top_N_plots} {input} {output}) > {log} 2>&1"


# rule tandem_repeat_coverage_dropouts:
#     input:
#         bam = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam",
#         bai = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam.bai",
#         bed = config['ref']['tg_bed']
#     output: f"samples/{sample}/tandem-genotypes/{sample}.tandem-genotypes.dropouts.txt"
#     log: f"samples/{sample}/logs/tandem-genotypes/{sample}.dropouts.log"
#     benchmark: f"samples/{sample}/benchmarks/tandem-genotypes/{sample}.dropouts.tsv"
#     conda: "envs/tandem-genotypes.yaml"
#     message: "Executing {rule}: Identify coverage dropouts in {input.bed} regions in {input.bam}."
#     shell: "(python3 workflow/scripts/check_tandem_repeat_coverage.py {input.bed} {input.bam} > {output}) > {log} 2>&1"