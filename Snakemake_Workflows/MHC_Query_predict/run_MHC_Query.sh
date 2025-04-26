#!/bin/bash

module load snakemake/7.7.0

snakemake \
--cluster-config ../config/cluster.json \
--cluster "sbatch --partition=norm,quick --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.nthreads} --gres=lscratch:100 " \
--jobs 1000 \
--latency-wait 240 \
--keep-going \
--max-jobs-per-second 1 \
--max-status-checks-per-second 0.1 \
--rerun-incomplete
