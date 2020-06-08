# snakemaker
Snakemake pipeline to run genome annotation with MAKER

Clone the repository with:
```bash
git clone https://github.com/chrishah/snakemaker.git
```

## **Prerequisites**

- A Linux cluster
- globally installed SLURM 18.08.7.1
- globally installed singularity 3.4.1+ 
- installed snakemake 5.10.0 (eg. in an anaconda environment - see `data/HOWTO_setup_conda_environment`)


## **Run the pipeline**


First you need to set up maker. This assumes you have a tarball of a current maker version (which you've downloaded yourself) deposited somewhere and specified the path to this tarball in the `data/config.yaml` file. If you have Repeatlibraries from Repbase these can also be setup in this step. You'll need to specify the path to the respective tarball also in the config file.
```bash
snakemake --use-singularity --singularity-args "-B $(pwd)" -p --until setup_maker
```

Once this is done you can just start the pipeline.

External evidence (e.g., proteins, transcriptomes) can be deposited in the respective directories in `data/`)

Regular run on single node:
```bash
snakemake --use-singularity -s Snakefile
```

Distribute on cluster, e.g. on VSC4 (which uses SLURM):
```bash
snakemake -s Snakefile \
        --jobs 2000 --latency-wait 120 \
        --use-singularity --singularity-args "-B $(pwd) -B $(pwd)/bin/RepeatMasker:/usr/local/RepeatMasker" \
        --cluster-config data/cluster_config_vsc4.yaml \
        --cluster '$(pwd)/bin/immediate_submit.py {dependencies}' \
        --immediate-submit --notemp -r -p
```


If you want to include Genemark in your analyses:
 - deposit the Genemark software directory somewhere
 - the directory is expected to contain your licence key in a file called `gm_key` in this directory
 - add the full path to the Genemark directory to your Snakemake call (see below)

Call the pipeline, including Genemark:
```bash
#this is where Genemark directory is
genemark="$(pwd)/data/external/gmes_linux_64/"
#call pipeline
snakemake -s Snakefile \
        --jobs 2000 --latency-wait 120 \
        --use-singularity --singularity-args "-B $(pwd) -B $genemark:/usr/local/Genemark -B $(pwd)/bin/RepeatMasker:/usr/local/RepeatMasker" \
        --cluster-config data/cluster_config_vsc4.yaml \
        --cluster '$(pwd)/bin/immediate_submit.py {dependencies}' \
        --immediate-submit --notemp -r -p
```

Make the rulegraph:
```bash
snakemake --rulegraph > temp
cat temp | grep "digraph" -A 1000 | dot -Tpng > rulegraph.png 
rm temp
```

## **Rulegraph**

## Rulegraph

<img src="https://github.com/chrishah/snakemaker/blob/master/Snakefile/rulegraph.png" eight="500">
