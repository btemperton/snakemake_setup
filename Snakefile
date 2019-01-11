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

samples = pd.read_csv('samples.csv',  comment='#', index_col=0)

rule all:
	input:
		expand(["/gpfs/ts0/home/bt273/BIOS-SCOPE/metag/shared_store/biller/reads/{sample}/{sample}.fwd.fq.gz",
		"/gpfs/ts0/home/bt273/BIOS-SCOPE/metag/shared_store/biller/reads/{sample}/{sample}.rev.fq.gz",
		"/gpfs/ts0/home/bt273/BIOS-SCOPE/metag/shared_store/biller/reads/{sample}/{sample}.contigs.fa",
		"/gpfs/ts0/home/bt273/BIOS-SCOPE/metag/shared_store/biller/reads/{sample}/{sample}.contigs.10k.fa.gz"], sample=samples.index),
		"combined_10k.contigs.fa.gz",
		"combined.contigs.virsorter.10k.tgz"

rule download_fwd:
	output:
		"/gpfs/ts0/home/bt273/BIOS-SCOPE/metag/shared_store/biller/reads/{sample}/{sample}.fwd.fq.gz"
	params:
		url = lambda wildcards: samples.fwd_download[wildcards.sample]
	shell:
		"wget -O {output} {params.url}"

rule download_rev:
	output:
		"/gpfs/ts0/home/bt273/BIOS-SCOPE/metag/shared_store/biller/reads/{sample}/{sample}.rev.fq.gz"
	params:
		url = lambda wildcards: samples.rev_download[wildcards.sample]
	shell:
		"wget -O {output} {params.url}"

rule download_contigs:
	output:
		"/gpfs/ts0/home/bt273/BIOS-SCOPE/metag/shared_store/biller/reads/{sample}/{sample}.contigs.fa"
	params:
		url = lambda wildcards: samples.contig_download[wildcards.sample]
	shell:
		"wget -O {output} {params.url}"

rule rename_and_size_select:
	input:
		rules.download_contigs.output
	output:
		"/gpfs/ts0/home/bt273/BIOS-SCOPE/metag/shared_store/biller/reads/{sample}/{sample}.contigs.10k.fa.gz"
	params:
		min_len = 10000
	script:
		"scripts/rename_and_size_select.py"

rule combine_contigs:
	output:
		"combined_10k.contigs.fa.gz"
	shell:
		"""
			find /gpfs/ts0/home/bt273/BIOS-SCOPE/metag/shared_store/biller/reads/ -name \"*.contigs.10k.fa.gz\" > file_list;
			for i in $(cat file_list); do cat $i >> tmp; done;
			mv tmp {output};
			rm file_list
		"""

rule run_virsorter:
	input:
		rules.combine_contigs.output
	output:
		"combined.contigs.virsorter.10k.tgz"
	params:
		data_dir = "/gpfs/ts0/home/bt273/BIOS-SCOPE/tools/virsorter-data"
	threads:
		16
	shell:
		"""
		gunzip -c {input} > tmp;
		wrapper_phage_contigs_sorter_iPlant.pl -f tmp --dataset "Biller_VS" --db 2 --wdir tmp.virsorter --ncpu {threads} --data-dir {params.data_dir} --diamond;
		rm tmp;
		tar czvf {output} tmp.virsorter;
		rm -r tmp.virsorter
		"""