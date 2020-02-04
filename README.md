# snakemaker
Snakemake pipeline to run genome annotation with MAKER

Clone the repository with:
```bash
git clone https://github.com/chrishah/snakemaker.git
```

Add evidence into the respective directories in `./data/`.

Regular run on single node:
```bash
snakemake --use-singularity -s Snakefile
```
There are submission scripts for the different clusters:
```bash
data/snakemake_submit_VSC3.slurm.sh
data/snakemake_submit_VSC4.slurm.sh
```

Distribute on cluster, e.g. on VSC4:
```bash
snakemake -s Snakefile_immediate_submit \
        --jobs 100 --latency-wait 30 \
        --use-singularity \
        --cluster-config data/cluster_config_vsc4.yaml \
        --cluster '$(pwd)/bin/immediate_submit.py {dependencies}' \
        --immediate-submit --notemp -pr
```
