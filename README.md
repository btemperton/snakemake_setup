# A Test Case for Running Snakemake on ISCA

This git represents a simple test case to show how to get snakemake to play nice with ISCA.
It consists of 4 important files:

1. `Snakefile` - This is the file that contains the rules that are to be run
2. `cluster.json` - This describes the parameters for each rule and how they are to be run on the servers
3. `config.yaml` - This contains the "samples". The snakefile will run each rule for each sample.
4. `custom_jobscript.sh` - This is a hack file to make it play nice with ISCA.

The reason you need `custom_jobscript.sh` is that otherwise, the created jobs don't know where `conda` is, so they can't run the programs you want. You are supposed to be able to use the flag `--use-conda` to get snakemake to create an appropriate environment from the `envs/test.yaml` file assigned to a rule. But it doesn't work.

The downside is that you will need to change the `custom_jobscript.sh` to point at your `conda` install, and on the next line call the `conda activate` on the appropriate environment. Once that's done, it works.

You can then get it to submit jobs to the server by running:

```
snakemake --jobs 2 \
  --jobscript custom_jobscript.sh \
  --cluster-config cluster.json \
  --cluster "msub -V -l walltime={cluster.walltime},nodes=1:ppn={cluster.threads} -q {cluster.queue} -A {cluster.project} -j oe"
```

Snakemake will then figure out that for each of the 5 samples specified in the `config.yaml`, the `test_conda` rule is needed to generate the outputs specified in the `all` rule. Because there is a dependency between the `test_conda` rule and the `setup` rule (because the input for the former is the output for the latter), it also figures out it needs to run `setup` for each sample.

The `jobs` parameter in the above command determines how many active jobs it can have. As there are 5 samples and 2 jobs, it will submit the task defined in `setup` for the first two jobs to the single queue (as specfied in `cluster.json`). Once these jobs complete, it will then submit the `test_conda` job to the parallel queue for the first two samples.

Once the first two samples are complete, it will move on to the next two and repeat. If you changed the `--jobs` flag to `5`, it would run `setup` for all jobs at once.

You'll notice on line 22 of the `Snakefile`, you can specify which rules you want to run locally, and which ones you want to turn into cluster jobs. If you uncomment this line, the `setup` rule will be run locally and only the `test_conda` rule will create a cluster job.

It's a pretty cool tool for creating reproducible computational pipelines.
