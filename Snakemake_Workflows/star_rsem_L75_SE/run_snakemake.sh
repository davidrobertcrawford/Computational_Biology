#!/bin/bash
module load snakemake

snakemake --use-conda --cluster-config ./config/cluster.json --cluster "sbatch --partition=norm,ccr --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.nthreads} --gres=lscratch:{cluster.lscratch}" --jobs 1000 --latency-wait 60 --keep-going --local-cores 4 all

