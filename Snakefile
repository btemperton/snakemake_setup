configfile: "config.yaml"
#This is simply a list of the outputs that are going to be needed for each element in the config file
#It is NOT the inputs to each job...
#For this example, there is no input

#We need a cluster.json script to define the parameters to use for our cluster
#This defines the walltime, queue and research project for each rule

#Note, we also have to use a custom jobscript, called custom_jobscript.sh
#You are supposed to be able to use --use-conda to initiate a conda environment
#from the env/env.yaml file
#However, on our servers, it can't find conda if you do this. The only way
#I've found to get round it is to create a custom job script that as its first lines
#initiates conda and then activates the conda environment. This is not ideal.


#we can use this to run the following for clustering:
#snakemake --jobs 1 --jobscript custom_jobscript.sh --cluster-config cluster.json --cluster "msub -V -l walltime={cluster.walltime},nodes=1:ppn={cluster.threads} -q {cluster.queue} -A {cluster.project} -j oe"




rule all:
    input:
        expand("{sample}.minimap.txt", sample=config["samples"])

rule setup:
    output:
        "{sample}.new.txt"
    conda:
        "envs/test.yaml"
    shell:
        "echo hello > {output}"

rule test_conda:
    input:
        rules.setup.output
    output:
        "{sample}.minimap.txt",
    shell:
        "minimap2 --version > {output}"
