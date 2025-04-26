#!/bin/bash
module load snakemake/7.7.0

snakemake --cluster-config ./config/cluster.json --cluster "sbatch --partition=norm,quick --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.nthreads} --gres=lscratch:{cluster.lscratch}" --jobs 1000 --latency-wait 60 --keep-going --rerun-incomplete --local-cores 4 all

