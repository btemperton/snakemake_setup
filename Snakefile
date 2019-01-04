###
#Run with the following:
#snakemake --jobs 100 \
# --use-conda \
#  --cluster-config cluster.json \
#  --cluster "msub -V -l walltime={cluster.walltime},nodes=1:ppn={cluster.threads} -q {cluster.queue} -A {cluster.project} -j oe"

# BATS samples are found at: https://www.imicrobe.us/#/projects/276

###
import pandas as pd
configfile: "config.yaml"

	samples = pd.read_csv('samples_complete.csv',  comment='#', index_col=0)
	print(samples.columns)

rule all:
	input:
		expand(["/gpfs/ts0/home/bt273/BIOS-SCOPE/metag/shared_store/biller/reads/{sample}/{sample}.fwd.fq.gz",
		"/gpfs/ts0/home/bt273/BIOS-SCOPE/metag/shared_store/biller/reads/{sample}/{sample}.rev.fq.gz"], sample=samples.index)

rule download_fwd:
	output:
		"/gpfs/ts0/home/bt273/BIOS-SCOPE/metag/shared_store/biller/reads/{sample}/{sample}.fwd.fq.gz"
	params:
		url = lambda wildcards: samples.fwd[wildcards.sample]
	shell:
		"wget -O {output} {params.url}"

rule download_rev:
	output:
		"/gpfs/ts0/home/bt273/BIOS-SCOPE/metag/shared_store/biller/reads/{sample}/{sample}.rev.fq.gz"
	params:
		url = lambda wildcards: samples.rev[wildcards.sample]
	shell:
		"wget -O {output} {params.url}"
