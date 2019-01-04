###
#Run with the following:
#snakemake --jobs 100 \
#  --jobscript custom_jobscript.sh \
#  --cluster-config cluster.json \
#  --cluster "msub -V -l walltime={cluster.walltime},nodes=1:ppn={cluster.threads} -q {cluster.queue} -A {cluster.project} -j oe"

# BATS samples are found at: https://www.imicrobe.us/#/projects/276

###
import pandas as pd
configfile: "config.yaml"

	samples = pd.read_table(config["samples"],  comment='#').set_index("sample")

rule all:
    input:
		expand(["/gpfs/ts0/home/bt273/BIOS-SCOPE/metag/shared_store/biller/reads/{location}/{sample}.fwd.fq.gz",
	"/gpfs/ts0/home/bt273/BIOS-SCOPE/metag/shared_store/biller/reads/{location}/{sample}.rev.fq.gz"],
	sample=samples.index, location=samples.Location)

rule download_fwd:
	output:
		"/gpfs/ts0/home/bt273/BIOS-SCOPE/metag/shared_store/biller/reads/{location}/{sample}.fwd.fq.gz"
	params:
		url = lambda wildcards: samples.fwd[wildcards.sample],
	shell:
		"wget -O {input} {params.url}"


rule download_rev:
	output:
		"/gpfs/ts0/home/bt273/BIOS-SCOPE/metag/shared_store/biller/reads/{location}/{sample}.rev.fq.gz"
	params:
		url = lambda wildcards: samples.rev[wildcards.sample],
	shell:
		"wget -O {input} {params.url}"


rule download_contigs:
	output:
		"/gpfs/ts0/home/bt273/BIOS-SCOPE/metag/shared_store/biller/contigs/{location}/{sample}.contigs.fa.gz"
	params:
		url = lambda wildcards: samples.contigs[wildcards.sample]
	shell:
		"wget -O {input} {params.url}"

rule size_filter_contigs:
	input:
		rules.download_contigs.output
	output:
		"/gpfs/ts0/home/bt273/BIOS-SCOPE/metag/shared_store/biller/contigs/10k/{location}/{sample}.contigs.10k.fa.gz"
	params:
		min_len = 10000
	conda:
		"envs/biller_virome.yaml"
	script:
		"scripts/rename_and_size_select.py"

