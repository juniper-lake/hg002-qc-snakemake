# test data
# takes dict of dicts to permit multiple movies per condition in format ubam_dict[condition][movie]
ubam_dict = { 'condition1' : {  'movie1': 'test_data/m64012_190920_173625.ccs.bam',
                                'movie2': 'test_data/m64012_190921_234837.ccs.bam'}}

# resources
reference = 'resources/human_GRCh38_no_alt_analysis_set.fasta'
reference_index = 'resources/human_GRCh38_no_alt_analysis_set.fasta.fai'
chromosome_lengths = 'resources/human_GRCh38_no_alt_analysis_set.chr_lengths.txt'
hg002_sr_yak = 'resources/HG002.sr.yak'
tandem_repeats_bed = 'resources/human_GRCh38_no_alt_analysis_set.trf.bed'
truvari_vcf = 'resources/GRCh38.HG002.latest.sv.vcf.gz'
truvari_bed = 'resources/GRCh38.HG002.latest.sv.bed'
happy = { 'all' : { 'bed' : 'resources/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed.gz',
                    'vcf' : 'resources/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz'},
           'cmrg' : { 'bed' : 'resources/HG002_GRCh38_CMRG_smallvar_v1.00.bed.gz',
                      'vcf' : 'resources/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz'}}
stratifications = 'resources/stratifications/v2.0-GRCh38-stratifications.tsv'
sdf = 'resources/human_GRCh38_no_alt_analysis_set.sdf'
all_chroms = ['chr'+ s for s in list(map(str,*[range(1,23)]))] + ['chrX','chrY','chrM']
deepvariant_version = '1.2.0-gpu'  # GPU
n_shards = 256
shards = [f"{x:05}" for x in range(n_shards)]


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


targets = [f"conditions/{condition}/{filename}"
                    for condition in ubam_dict.keys()
                    for filename in ["smrtcell_stats/all_movies.read_length_and_quality.tsv",
                                    "hifiasm/asm.p_ctg.fasta.stats.txt",
                                    "hifiasm/asm.a_ctg.fasta.stats.txt",
                                    "hifiasm/asm.p_ctg.qv.txt",
                                    "hifiasm/asm.a_ctg.qv.txt",
                                    "pbsv/all_chroms.pbsv.vcf.gz",
                                    "truvari/summary.txt",
                                    "deepvariant/deepvariant.vcf.stats.txt",
                                    "whatshap/deepvariant.phased.tsv",
                                    "whatshap/deepvariant.phased.blocklist",
                                    "happy/all.summary.csv",
                                    "happy/all.extended.csv",
                                    "happy/cmrg.summary.csv",
                                    "happy/cmrg.extended.csv",
                                    "mosdepth/coverage.mosdepth.summary.txt",
                                    "mosdepth/mosdepth.M2_ratio.txt",
                                    "mosdepth/gc_coverage.summary.txt",
                                    "mosdepth/coverage.thresholds.summary.txt"]]

# targets = [f"conditions/{condition}/{filename}"
#                     for condition in ubam_dict.keys()
#                     for filename in ["truvari/truvari.summary.txt"]]


rule all:
    input: targets


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


rule asm_stats:
    input: 
        fasta = "conditions/{condition}/hifiasm/asm.{infix}.fasta.gz",
        ref_index = reference_index
    output: "conditions/{condition}/hifiasm/asm.{infix}.fasta.stats.txt"
    log: "conditions/{condition}/logs/asm_stats.{infix}.log"
    conda: "envs/k8.yaml"
    shell: "(k8 workflow/scripts/calN50.js -f {input.ref_index} {input.fasta} > {output}) > {log} 2>&1"


rule yak_qv:
    input:
        asm = "conditions/{condition}/hifiasm/asm.{infix}.fasta.gz",
        sr = hg002_sr_yak
    output: "conditions/{condition}/hifiasm/asm.{infix}.qv.txt"
    log: "conditions/{condition}/logs/yak_qv.{infix}.log"
    threads: 24
    conda: "envs/yak.yaml"
    shell: "(yak qv -t{threads} -p -K3.2g -l100k {input.sr} <(zcat {input.asm}) > {output}) > {log} 2>&1"



rule bgzip_vcf:
    input: "conditions/{condition}/{prefix}.vcf"
    output: "conditions/{condition}/{prefix}.vcf.gz"
    log: "conditions/{condition}/logs/bgzip_vcf/{prefix}.log"
    threads: 2
    conda: "envs/htslib.yaml"
    shell: "(bgzip --threads {threads} {input}) > {log} 2>&1"


rule tabix_vcf:
    input: "conditions/{condition}/{prefix}.vcf.gz"
    output: "conditions/{condition}/{prefix}.vcf.gz.tbi"
    log: "conditions/{condition}/logs/tabix_vcf/{prefix}.log"
    params: "-p vcf"
    conda: "envs/htslib.yaml"
    shell: "tabix {params} {input} > {log} 2>&1"


rule pbsv_discover:
    input:
        bam = rules.pbmm2_align.output.bam,
        bai = rules.pbmm2_align.output.bai,
        tr_bed = tandem_repeats_bed
    output: "conditions/{condition}/pbsv/svsig/{movie}.{chromosome}.svsig.gz"
    log: "conditions/{condition}/logs/pbsv_discover.{movie}.{chromosome}.log"
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
        svsigs = lambda wc: [f"conditions/{wc.condition}/pbsv/svsig/{movie}.{wc.chromosome}.svsig.gz" for movie in ubam_dict[wc.condition].keys()],
        ref = reference
    output: "conditions/{condition}/pbsv/chrom_vcfs/{chromosome}.pbsv.vcf"
    log: "conditions/{condition}/logs/pbsv_call.{chromosome}.log"
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
        calls = lambda wc: expand(f"conditions/{wc.condition}/pbsv/chrom_vcfs/{{chromosome}}.pbsv.vcf.gz", chromosome=all_chroms),
        indices = lambda wc: expand(f"conditions/{wc.condition}/pbsv/chrom_vcfs/{{chromosome}}.pbsv.vcf.gz.tbi", chromosome=all_chroms)
    output: "conditions/{condition}/pbsv/all_chroms.pbsv.vcf"
    log: "conditions/{condition}/logs/bcftools_concat_pbsv_vcf.log"
    conda: "envs/bcftools.yaml"
    shell: "(bcftools concat -a -o {output} {input.calls}) > {log} 2>&1"


rule truvari_benchmark:
    input:
        ref = reference,
        query_vcf = "conditions/{condition}/pbsv/all_chroms.pbsv.vcf.gz",
        query_tbi = "conditions/{condition}/pbsv/all_chroms.pbsv.vcf.gz.tbi",
        bench_vcf = truvari_vcf,
        bench_bed = truvari_bed
    output: 
        "conditions/{condition}/truvari/tp-call.vcf",
        "conditions/{condition}/truvari/tp-base.vcf",
        "conditions/{condition}/truvari/fn.vcf",
        "conditions/{condition}/truvari/fp.vcf",
        "conditions/{condition}/truvari/base-filter.vcf",
        "conditions/{condition}/truvari/call-filter.vcf",
        "conditions/{condition}/truvari/summary.txt",
        "conditions/{condition}/truvari/log.txt",
        "conditions/{condition}/truvari/giab_report.txt"
    log: "conditions/{condition}/logs/truvari_benchmark.log"
    params: prefix = "conditions/{condition}/truvari"
    conda: "envs/truvari.yaml"
    shell:
        """
        (rmdir {params.prefix} && \
        truvari \
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
        bams = lambda wc: [f"conditions/{{condition}}/aligned/{movie}.aligned.bam" for movie in ubam_dict[wc.condition].keys()],
        bais = lambda wc: [f"conditions/{{condition}}/aligned/{movie}.aligned.bam.bai" for movie in ubam_dict[wc.condition].keys()],
        ref = reference
    output:
        tfrecord = f"conditions/{{condition}}/deepvariant_intermediate/examples/examples.tfrecord-{{shard}}-of-{n_shards:05}.gz"
    log: f"conditions/{{condition}}/logs/deepvariant_make_examples_round1.{{shard}}-of-{n_shards:05}.log"
    container: f"docker://google/deepvariant:{deepvariant_version}"
    params:
        vsc_min_fraction_indels = "0.12",
        pileup_image_width = 199,
        shard = lambda wc: wc.shard,
        reads = lambda wc: ','.join([f"conditions/{wc.condition}/aligned/{movie}.aligned.bam" for movie in ubam_dict[wc.condition].keys()])
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
            --examples conditions/{{wildcards.condition}}/deepvariant_intermediate/examples/examples.tfrecord@{n_shards}.gz \
            --task {{wildcards.shard}}) > {{log}} 2>&1
        """


rule deepvariant_call_variants_gpu_round1:
    input: lambda wc: expand(f"conditions/{wc.condition}/deepvariant_intermediate/examples/examples.tfrecord-{{shard}}-of-{n_shards:05}.gz", shard=shards)
    output: "conditions/{condition}/deepvariant_intermediate/call_variants_output.tfrecord.gz"
    log: "conditions/{condition}/logs/deepvariants_call_variants_gpu_round1.log"
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
            --examples conditions/{{wildcards.condition}}/deepvariant_intermediate/examples/examples.tfrecord@{n_shards}.gz \
            --checkpoint {{params.model}}) > {{log}} 2>&1
        """


rule deepvariant_postprocess_variants_round1:
    input:
        tfrecord = rules.deepvariant_call_variants_gpu_round1.output,
        ref = reference
    output:
        vcf = "conditions/{condition}/deepvariant_intermediate/deepvariant.vcf.gz",
        vcf_index = "conditions/{condition}/deepvariant_intermediate/deepvariant.vcf.gz.tbi",
        report = "conditions/{condition}/deepvariant_intermediate/deepvariant.visual_report.html"
    log: "conditions/{condition}/logs/deepvariant_postprocess_variants_round1.log"
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
    output: "conditions/{condition}/whatshap_intermediate/{chromosome}.deepvariant.vcf"
    log: "conditions/{condition}/logs/split_deepvariant_vcf_round1.{chromosome}.log"
    params: extra = '-h'
    conda: "envs/htslib.yaml"
    shell: "tabix {params.extra} {input} {wildcards.chromosome} > {output} 2> {log}"


rule whatshap_phase_round1:
    input:
        ref = reference,
        vcf = "conditions/{condition}/whatshap_intermediate/{chromosome}.deepvariant.vcf.gz",
        tbi = "conditions/{condition}/whatshap_intermediate/{chromosome}.deepvariant.vcf.gz.tbi",
        phaseinput = lambda wc: [f"conditions/{{condition}}/aligned/{movie}.aligned.bam" for movie in ubam_dict[wc.condition].keys()],
        phaseinputindex = lambda wc: [f"conditions/{{condition}}/aligned/{movie}.aligned.bam.bai" for movie in ubam_dict[wc.condition].keys()],
    output: "conditions/{condition}/whatshap_intermediate/{chromosome}.deepvariant.phased.vcf.gz"
    log: "conditions/{condition}/logs/whatshap_phase_round1.{chromosome}.log"
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
        calls = lambda wc: expand(f"conditions/{wc.condition}/whatshap_intermediate/{{chromosome}}.deepvariant.phased.vcf.gz", chromosome=all_chroms),
        indices = lambda wc: expand(f"conditions/{wc.condition}/whatshap_intermediate/{{chromosome}}.deepvariant.phased.vcf.gz.tbi", chromosome=all_chroms)
    output: "conditions/{condition}/whatshap_intermediate/deepvariant.phased.vcf.gz"
    log: "conditions/{condition}/logs/whatshap_bcftools_concat_round1.log"
    params: "-a -Oz"
    conda: "envs/bcftools.yaml"
    shell: "bcftools concat {params} -o {output} {input.calls} > {log} 2>&1"


rule whatshap_haplotag_round1:
    input:
        ref = reference,
        vcf = rules.whatshap_bcftools_concat_round1.output,
        tbi = "conditions/{condition}/whatshap_intermediate/deepvariant.phased.vcf.gz.tbi",
        bam = rules.pbmm2_align.output.bam
    output: "conditions/{condition}/whatshap_intermediate/{movie}.deepvariant.haplotagged.bam"
    log: "conditions/{condition}/logs/whatshap_haplotag_round1.{movie}.log"
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
    input: "conditions/{condition}/{whatshap_folder}/{prefix}.bam"
    output: "conditions/{condition}/{whatshap_folder}/{prefix}.bam.bai"
    log: "conditions/{condition}/logs/samtools_index_bam/{whatshap_folder}/{prefix}.log"
    threads: 4
    conda: "envs/samtools.yaml"
    shell: "(samtools index -@ 3 {input}) > {log} 2>&1"


rule deepvariant_make_examples_round2:
    input:
        bams = lambda wc: [f"conditions/{{condition}}/whatshap_intermediate/{movie}.deepvariant.haplotagged.bam" for movie in ubam_dict[wc.condition].keys()],
        bais = lambda wc: [f"conditions/{{condition}}/whatshap_intermediate/{movie}.deepvariant.haplotagged.bam.bai" for movie in ubam_dict[wc.condition].keys()],
        ref = reference
    output:
        tfrecord = f"conditions/{{condition}}/deepvariant/examples/examples.tfrecord-{{shard}}-of-{n_shards:05}.gz",
        nonvariant_site_tfrecord = f"conditions/{{condition}}/deepvariant/examples/gvcf.tfrecord-{{shard}}-of-{n_shards:05}.gz"
    log: f"conditions/{{condition}}/logs/deepvariant_make_examples_round2.{{shard}}-of-{n_shards:05}.log"
    container: f"docker://google/deepvariant:{deepvariant_version}"
    params:
        vsc_min_fraction_indels = "0.12",
        pileup_image_width = 199,
        shard = lambda wc: wc.shard,
        reads = lambda wc: ','.join([f"conditions/{wc.condition}/whatshap_intermediate/{movie}.deepvariant.haplotagged.bam" for movie in ubam_dict[wc.condition].keys()])
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
            --examples conditions/{{wildcards.condition}}/deepvariant/examples/examples.tfrecord@{n_shards}.gz \
            --gvcf conditions/{{wildcards.condition}}/deepvariant/examples/gvcf.tfrecord@{n_shards}.gz \
            --task {{wildcards.shard}}) > {{log}} 2>&1
        """


rule deepvariant_call_variants_gpu_round2:
    input: lambda wc: expand(f"conditions/{wc.condition}/deepvariant/examples/examples.tfrecord-{{shard}}-of-{n_shards:05}.gz", shard=shards)
    output: "conditions/{condition}/deepvariant/call_variants_output.tfrecord.gz"
    log: "conditions/{condition}/logs/deepvariant_call_variants_gpu_round2.log"
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
            --examples conditions/{{wildcards.condition}}/deepvariant/examples/examples.tfrecord@{n_shards}.gz \
            --checkpoint {{params.model}}) > {{log}} 2>&1
        """


rule deepvariant_postprocess_variants_round2:
    input:
        tfrecord = "conditions/{condition}/deepvariant/call_variants_output.tfrecord.gz",
        nonvariant_site_tfrecord = lambda wc: expand(f"conditions/{wc.condition}/deepvariant/examples/gvcf.tfrecord-{{shard:05}}-of-{n_shards:05}.gz",
                                          shard=range(n_shards)),
        ref = reference
    output:
        vcf = "conditions/{condition}/deepvariant/deepvariant.vcf.gz",
        vcf_index = "conditions/{condition}/deepvariant/deepvariant.vcf.gz.tbi",
        gvcf = "conditions/{condition}/deepvariant/deepvariant.g.vcf.gz",
        gvcf_index = "conditions/{condition}/deepvariant/deepvariant.g.vcf.gz.tbi",
        report = "conditions/{condition}/deepvariant/deepvariant.visual_report.html"
    log: "conditions/{condition}/logs/deepvariant_postprocess_variants_round2.log"
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
            --nonvariant_site_tfrecord_path conditions/{{wildcards.condition}}/deepvariant/examples/gvcf.tfrecord@{n_shards}.gz \
            --gvcf_outfile {{output.gvcf}}) > {{log}} 2>&1
        """


rule happy_benchmark_deepvariant:
    input:
        ref = reference,
        query_vcf = rules.deepvariant_postprocess_variants_round2.output.vcf,
        query_tbi = rules.deepvariant_postprocess_variants_round2.output.vcf_index,
        bench_vcf = lambda wc: happy[wc.ver]['vcf'],
        bench_bed = lambda wc: happy[wc.ver]['bed'],
        strats = stratifications,
        sdf = sdf
    output:
        "conditions/{condition}/happy/{ver}.extended.csv",
        "conditions/{condition}/happy/{ver}.metrics.json.gz",
        "conditions/{condition}/happy/{ver}.roc.all.csv.gz",
        "conditions/{condition}/happy/{ver}.roc.Locations.INDEL.csv.gz",
        "conditions/{condition}/happy/{ver}.roc.Locations.INDEL.PASS.csv.gz",
        "conditions/{condition}/happy/{ver}.roc.Locations.SNP.csv.gz",
        "conditions/{condition}/happy/{ver}.roc.Locations.SNP.PASS.csv.gz",
        "conditions/{condition}/happy/{ver}.runinfo.json",
        "conditions/{condition}/happy/{ver}.summary.csv",
        "conditions/{condition}/happy/{ver}.vcf.gz",
        "conditions/{condition}/happy/{ver}.vcf.gz.tbi"
    log: "conditions/{condition}/logs/happy_benchmark_deepvariant.{ver}.log"
    container: "docker://pkrusche/hap.py:latest"
    params:
        prefix = "conditions/{condition}/happy/{ver}"
    threads: 12
    shell:
        """
        (/opt/hap.py/bin/hap.py \
            --threads {threads} \
            -r {input.ref} -f {input.bench_bed} \
            -o {params.prefix} \
            --engine=vcfeval --engine-vcfeval-template {input.sdf} \
            --stratification {input.strats} \
            {input.bench_vcf} {input.query_vcf}) > {log} 2>&1
        """


rule deepvariant_bcftools_stats:
    input: rules.deepvariant_postprocess_variants_round2.output.vcf
    output: "conditions/{condition}/deepvariant/deepvariant.vcf.stats.txt"
    log: "conditions/{condition}/logs/deepvariant_bcftools_stats.log"
    params: f"--fasta-ref {reference} --apply-filters PASS -s {{condition}}"
    threads: 4
    conda: "envs/bcftools.yaml"
    shell: "(bcftools stats --threads 3 {params} {input} > {output}) > {log} 2>&1"


rule split_deepvariant_vcf_round2:
    input: rules.deepvariant_postprocess_variants_round2.output.vcf
    output: "conditions/{condition}/whatshap/{chromosome}.deepvariant.vcf"
    log: "conditions/{condition}/logs/split_deepvariant_vcf_round2.{chromosome}.log"
    params: region = lambda wc: wc.chromosome, extra = '-h'
    conda: "envs/htslib.yaml"
    shell: "tabix {params.extra} {input} {params.region} > {output} 2> {log}"


rule whatshap_phase_round2:
    input:
        ref = reference,
        vcf = "conditions/{condition}/whatshap/{chromosome}.deepvariant.vcf.gz",
        tbi = "conditions/{condition}/whatshap/{chromosome}.deepvariant.vcf.gz.tbi",
        phaseinput = lambda wc: [f"conditions/{wc.condition}/aligned/{movie}.aligned.bam" for movie in ubam_dict[wc.condition].keys()],
        phaseinputindex = lambda wc: [f"conditions/{wc.condition}/aligned/{movie}.aligned.bam.bai" for movie in ubam_dict[wc.condition].keys()]
    output: "conditions/{condition}/whatshap/{chromosome}.deepvariant.phased.vcf.gz"
    log: "conditions/{condition}/logs/whatshap_phase_round2.{chromosome}.log"
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
        calls = lambda wc: expand(f"conditions/{wc.condition}/whatshap/{{chromosome}}.deepvariant.phased.vcf.gz", chromosome=all_chroms),
        indices = lambda wc: expand(f"conditions/{wc.condition}/whatshap/{{chromosome}}.deepvariant.phased.vcf.gz.tbi", chromosome=all_chroms)
    output: "conditions/{condition}/whatshap/deepvariant.phased.vcf.gz"
    log: "conditions/{condition}/logs/whatshap_bcftools_concat_round2.log"
    params: "-a -Oz"
    conda: "envs/bcftools.yaml"
    shell: "bcftools concat {params} -o {output} {input.calls} > {log} 2>&1"


rule whatshap_stats:
    input:
        vcf = rules.whatshap_bcftools_concat_round2.output,
        tbi = "conditions/{condition}/whatshap/deepvariant.phased.vcf.gz.tbi",
        chr_lengths = chromosome_lengths
    output:
        gtf = "conditions/{condition}/whatshap/deepvariant.phased.gtf",
        tsv = "conditions/{condition}/whatshap/deepvariant.phased.tsv",
        blocklist = "conditions/{condition}/whatshap/deepvariant.phased.blocklist"
    log: "conditions/{condition}/logs/whatshap_stats.log"
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
        tbi = "conditions/{condition}/whatshap/deepvariant.phased.vcf.gz.tbi",
        bam = rules.pbmm2_align.output.bam
    output: "conditions/{condition}/whatshap/{movie}.deepvariant.haplotagged.bam"
    log: "conditions/{condition}/logs/whatshap_haplotag_round2.{movie}.log"
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
    input: lambda wc: expand(f"conditions/{wc.condition}/whatshap/{{movie}}.deepvariant.haplotagged.bam", movie=ubam_dict[wc.condition].keys())
    output: "conditions/{condition}/whatshap/deepvariant.haplotagged.bam"
    log: "conditions/{condition}/logs/merge_haplotagged_bams.log"
    threads: 8
    conda: "envs/samtools.yaml"
    shell: "(samtools merge -@ 7 {output} {input}) > {log} 2>&1"


rule mosdepth:
    input:
        bam = rules.merge_haplotagged_bams.output,
        bai = "conditions/{condition}/whatshap/deepvariant.haplotagged.bam.bai"
    output:
        summary = "conditions/{condition}/mosdepth/coverage.mosdepth.summary.txt",
        regions = "conditions/{condition}/mosdepth/coverage.regions.bed.gz",
        thresholds = "conditions/{condition}/mosdepth/coverage.thresholds.bed.gz"
    log: "conditions/{condition}/logs/mosdepth.log"
    params:
        by = "500",
        prefix = "conditions/{condition}/mosdepth/coverage",
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
    output: "conditions/{condition}/mosdepth/mosdepth.M2_ratio.txt"
    log: "conditions/{condition}/logs/calculate_M2_ratio.log"
    conda: "envs/pandas.yaml"
    shell: "(python3 workflow/scripts/calculate_M2_ratio.py {input} > {output}) > {log} 2>&1"


rule calculate_gc_coverage:
    input:
        mosdepth_regions = rules.mosdepth.output.regions,
        ref = reference,
    output: "conditions/{condition}/mosdepth/gc_coverage.summary.txt"
    log: "conditions/{condition}/logs/calculate_gc_coverage.log"
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
    output: "conditions/{condition}/mosdepth/coverage.thresholds.summary.txt"
    log: "conditions/{condition}/logs/calculate_coverage_thresholds.log"
    conda: "envs/gc_coverage.yaml"
    shell: 
        """
        (zcat {input} | awk 'NR==1' | cut -f 5-14 > {output};
            zcat {input} \
                | datamash -H mean 5-14 \
                | awk -v c=500 'NR!=1 {{ for (i = 1; i <= NF; ++i) $i /= c; print }}' >> {output}) > {log} 2>&1
        """
