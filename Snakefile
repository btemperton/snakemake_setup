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
samples = pd.read_table(config["samples"]).set_index("sample")

rule all:
    input:
        expand(["samples/{sample}/{sample}.fwd.fq.bz2",
				"samples/{sample}/{sample}.rev.fq.bz2",
	"samples/{sample}/{sample}.contigs.10k.fa.gz",
	"samples/{sample}/{sample}.virsorter.tgz"],
				sample=samples.index)

rule download_fwd:
	output:
		"samples/{sample}/{sample}.fwd.fq.bz2"
	params:
		url = lambda wildcards: samples.dir[wildcards.sample] + "/" + samples.fwd[wildcards.sample],
		basename = lambda wildcards: samples.fwd[wildcards.sample]
	benchmark:
		"logs/{sample}/download_fwd.bmk"
	conda:
		"envs/biller_virome.yaml"
	shell:
		"""
		iget -K {params.url};
		mv {params.basename} {output}
		"""

rule download_rev:
	output:
		"samples/{sample}/{sample}.rev.fq.bz2"
	params:
		url = lambda wildcards: samples.dir[wildcards.sample] + "/" + samples.rev[wildcards.sample],
		basename = lambda wildcards: samples.rev[wildcards.sample]
	benchmark:
		"logs/{sample}/download_rev.bmk"
	conda:
		"envs/biller_virome.yaml"
	shell:
		"""
		iget -K {params.url};
		mv {params.basename} {output}
		"""
rule download_contigs:
	output:
		temp("samples/{sample}/{sample}.contigs.fa.gz")
	params:
		url = lambda wildcards: samples.dir[wildcards.sample] + "/" + samples.contigs[wildcards.sample],
		basename = lambda wildcards: samples.contigs[wildcards.sample]
	benchmark:
		"logs/{sample}/download_contigs.bmk"
	conda:
		"envs/biller_virome.yaml"
	shell:
		"""
		iget -K {params.url};
		mv {params.basename} {output}
		"""

rule size_filter_contigs:
	input:
		rules.download_contigs.output
	output:
		"samples/{sample}/{sample}.contigs.10k.fa.gz"
	params:
		min_len = 10000
	benchmark:
		"logs/{sample}/size_filter_contigs.bmk"
	conda:
		"envs/biller_virome.yaml"
	script:
		"scripts/rename_and_size_select.py"

rule run_virsorter:
	input:
		rules.size_filter_contigs.output
	output:
		"samples/{sample}/{sample}.virsorter.tgz"
	params:
		data_dir = "/gpfs/ts0/projects/Research_Project-172179/tools/VirSorter/virsorter-data"
	threads: 16
	log:
		"logs/{sample}/{sample}.virsorter.log"
	script:
			"""
		gunzip -c {input} > tmp.{wildcards.sample};
		wrapper_phage_contigs_sorter_iPlant.pl \
			-f tmp.{wildcards.sample} \
			--db 2 \
			--wdir {wildcards.sample}.virsorter \
			--ncpu {threads} \
			--data-dir {params.data_dir} 2>&1 | tee {log};

		tar czf {output} {wildcards.sample}.virsorter;
		rm -rf tmp.{wildcards.sample}.virsorter

		"""

