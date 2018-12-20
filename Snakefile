configfile: "config.yaml"
#This is simply a list of the outputs that are going to be needed for each element in the config file
#It is NOT the inputs to each job...
#For this example, there is no input

#we can use this to run the following for clustering:
#snakemake --jobs 1 --cluster-config cluster.json --cluster "qsub -l walltime={cluster.walltime},nodes=1:ppn={cluster.threads} -q {cluster.queue} -A {cluster.project} -j oe"
rule all:
    input:
        expand("{sample}.new.txt", sample=config["samples"])

rule setup:
    output:
        "{sample}.new.txt"
    conda:
        "envs/test.yaml"
    shell:
        "echo hello > {output}"
