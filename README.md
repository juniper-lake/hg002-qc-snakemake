# WGS Benchmark Snakemake

## To Run

**Warning: Several steps of this workflow require minimum coverage. It's recommended that this workflow not be run when yield in base pairs is insufficient to produce at least 15X coverage (i.e. yield/3099922541 >= 15x).**

```text
# clone repo
git clone --recursive https://github.com/juniper-lake/hg002-qc-snakemake.git workflow

# make necessary directories
mkdir cluster_logs resources

# put necessary resources into resources folder (check config file for list of resources)

# update config.yaml as necessary
nano workflow/config.yaml

# create sample sheet TSV
nano sample_sheet.tsv

# create conda environment
conda env create --file workflow/environment.yaml --prefix ./WGSbenchmark_env

# activate conda environment
conda activate ./WGSbenchmark_env

# dry run
sbatch workflow/run_WGSbenchmark.sh sample_sheet.tsv -n

# submit job
sbatch workflow/run_WGSbenchmark.sh sample_sheet.tsv
```

## Plots (NEEDS TO BE UPDATED)

A list of important stats from target files that would be good for plotting.

```text
targets = [f"conditions/{condition}/{filename}"
                    for condition in ubam_dict.keys()
                    for filename in ["smrtcell_stats/all_movies.read_length_and_quality.tsv",
                                    "hifiasm/asm.p_ctg.fasta.stats.txt",
                                    "hifiasm/asm.a_ctg.fasta.stats.txt",
                                    "hifiasm/asm.p_ctg.qv.txt",
                                    "hifiasm/asm.a_ctg.qv.txt",
                                    "truvari/summary.txt",
                                    "pbsv/all_chroms.pbsv.vcf.gz",
                                    "deepvariant/deepvariant.vcf.stats.txt",
                                    "whatshap/deepvariant.phased.tsv",
                                    "happy/all.summary.csv",
                                    "happy/all.extended.csv",
                                    "happy/cmrg.summary.csv",
                                    "happy/cmrg.extended.csv",
                                    "mosdepth/coverage.mosdepth.summary.txt",
                                    "mosdepth/mosdepth.M2_ratio.txt",
                                    "mosdepth/gc_coverage.summary.txt",
                                    "mosdepth/coverage.thresholds.summary.txt",
                                    "whatshap/phase.eval.tsv"]]
```

- `smrtcell_stats/all_movies.read_length_and_quality.tsv`
  - outputs 3 columns (read name, read length, read quality)
  - boxplots of read length and quality
- `hifiasm/asm.p_ctg.fasta.stats.txt` (primary) + `hifiasm/asm.a_ctg.fasta.stats.txt` (alternate)
  - all stats below should be collected for both primary (p_ctg) and alternate (p_atg) assemblies
  - assembly size `awk '$1=="SZ" {print $2}' <filename>`
  - auN (area under the curve) `awk '$1=="AU" {print $2}' <filename>`
  - NGx - line plot of NG10 through NG90 `awk '$1=="NL" {print $2 $3}' <filename>` ($2 is x-axis, $3 y-axis) like this: [example plot](http://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity)
- `hifiasm/asm.p_ctg.qv.txt` + `hifiasm/asm.a_ctg.qv.txt`
  - adjusted assembly quality `awk '$1=="QV" {print $3}' <filename>` for primary and alternate assemblies
- `truvari/truvari.summary.txt`
  - structural variant recall `jq .recall <filename>`
  - structural variant precision `jq .precision <filename>`
  - structural variant f1 `jq .f1 <filename>`
  - number of calls `jq '."call cnt"' <filename>`
  - FP `jq .FP <filename>`
  - TP-call `jq .TP-call <filename>`
  - FN `jq .FN <filename>`
  - TP-base `jq .TP-base <filename>`
- `pbsv/all_chroms.pbsv.vcf.gz`
  - counts of each type of variant `bcftools query -i 'FILTER=="PASS"' -f '%INFO/SVTYPE\n' <filename> | awk '{A[$1]++}END{for(i in A)print i,A[i]}'`
  - can also do size distributions of indels `bcftools query -i 'FILTER=="PASS" && (INFO/SVTYPE=="INS" | INFO/SVTYPE=="DEL")' -f '%INFO/SVTYPE\t%INFO/SVLEN\n' <filename>`
- `deepvariant/deepvariant.vcf.stats.txt`
  - several values in lines starting with 'SN' `awk '$1=="SN"' <filename>`
    - number of SNPS
    - number INDELs
    - number of multi-allelic sites
    - number of multi-allelic SNP sites
  - ratio of transitions to transversions `awk '$1=="TSTV" {print$5}' <filename>`
  - can monitor substitution types `awk '$1=="ST"' <filename>`
  - SNP heterozygous : non-ref homozygous ratio `awk '$1=="PSC" {print $6/$5}' <filename>`
  - SNP transitions : transversions `awk '$1=="PSC" {print $7/$8}' <filename>`
  - Number of heterozygous insertions : number of homozgyous alt insertions `awk '$1=="PSI" {print $8/$10}' <filename>`
  - Number of heterozygous deletions : number of homozgyous alt deletions `awk '$1=="PSI" {print $9/$11}' <filename>`
  - Total INDEL heterozygous:homozygous ratio `awk '$1=="PSI" {print ($8+$9)/($10+$11)}' <filename>`8+9:10+11 indel het:hom)
- `whatshap/deepvariant.phased.tsv`
  - phase block N50 `awk '$2=="ALL" {print $22}' <filename>`
  - bp_per_block_sum (total number of phased bases) `awk '$2=="ALL" {print $18}' <filename>`
- `whatshap/deepvariant.phased.blocklist`
  - calculate phase block size (to - from) and reverse order them (`awk 'NR>1 {print $5-$4}' <filename> |sort -nr`), then plot as cumulative line graph like for assembly, N_0 to N90 [example plot](http://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity)
- `happy/all.summary.csv` + `happy/cmrg.summary.csv`
  - stats should be collected for `all` variants and `cmrg` challenging medically relevant genes
    - SNP recall `awk -F, '$1=="SNP" && $2=="PASS" {print $10}' <filename>`
    - SNP precision `awk -F, '$1=="SNP" && $2=="PASS" {print $11}' <filename>`
    - SNP F1 `awk -F, '$1=="SNP" && $2=="PASS" {print $13}' <filename>`
    - INDEL recall `awk -F, '$1=="INDEL" && $2=="PASS" {print $10}' <filename>`
    - INDEL precision `awk -F, '$1=="INDEL" && $2=="PASS" {print $11}' <filename>`
    - INDEL F1 `awk -F, '$1=="INDEL" && $2=="PASS" {print $13}' <filename>`
- `happy/all.extended.csv` + `happy/cmrg.extended.csv`
  - The following commands are just for one stratification "GRCh38_lowmappabilityall.bed.gz". Important stratifications include:
    - GRCh38_lowmappabilityall.bed.gz
    - GRCh38_nonunique_l250_m0_e0.bed.gz
    - GRCh38_segdups_gt10kb.bed.gz
    - GRCh38_SimpleRepeat_homopolymer_gt11_slop5.bed.gz
    - GRCh38_AllTandemRepeats_gt100bp_slop5.bed.gz
  - SNP GRCh38_lowmappabilityall recall `awk -F, '$1=="SNP" && $2=="*" && $3=="GRCh38_lowmappabilityall.bed.gz" && $4=="PASS" {print $8}' <filename>`
  - SNP GRCh38_lowmappabilityall precision `awk -F, '$1=="SNP" && $2=="*" && $3=="GRCh38_lowmappabilityall.bed.gz" && $4=="PASS" {print $9}' <filename>`
  - SNP GRCh38_lowmappabilityall F1 `awk -F, '$1=="SNP" && $2=="*" && $3=="GRCh38_lowmappabilityall.bed.gz" && $4=="PASS" {print $11}' <filename>` 
  - INDEL GRCh38_lowmappabilityall recall `awk -F, '$1=="INDEL" && $2=="*" && $3=="GRCh38_lowmappabilityall.bed.gz" && $4=="PASS" {print $8}' <filename>`
  - INDEL GRCh38_lowmappabilityall precision `awk -F, '$1=="INDEL" && $2=="*" && $3=="GRCh38_lowmappabilityall.bed.gz" && $4=="PASS" {print $9}' <filename>`
  - INDEL GRCh38_lowmappabilityall F1 `awk -F, '$1=="INDEL" && $2=="*" && $3=="GRCh38_lowmappabilityall.bed.gz" && $4=="PASS" {print $11}' <filename>`
- `mosdepth/coverage.mosdepth.summary.txt`
  - mean aligned coverage in "coverage.mosdepth.summary.txt" - 4th column of final row, can grep 'total_region'
- `mosdepth/mosdepth.M2_ratio.txt`
  - outputs single value: ratio of chr2 coverage to chrM coverage
  - bar chart of m2 ratio
- `mosdepth/gc_coverage.summary.txt`
  - outputs 5 columns: gc percentage bin, q1 , median , q3 , count
  - q1, median, q3 columns are statistics for coverage at different gc percentages (e.g. median cover at 30% GC)
  - "count" refers to # of 500 bp windows that fall in that bin
  - can pick a couple of key GC coverage bins and make box plots out of them
- `mosdepth/coverage.thresholds.summary.txt`
  - outputs 10 columns corresponding to % of genome sequenced to minimum coverage depths (1X - 10X)
  - maybe a line chart comparing the different coverage thresholds among conditions
- `whatshap/phase.eval.tsv`
  - each row is a different chromosome
  - important values for each chromosome might include columns:
    - 9  all_assessed_pairs (sample size for error rates)
    - 11  all_switch_rate
    - 13  all_switchflip_rate
    - 15  blockwise_hamming_rate
