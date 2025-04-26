#!/bin/bash

module load snakemake

snakemake \
--use-conda \
--cluster-config ./config/cluster.json \
--cluster "sbatch --partition=norm,ccr --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.nthreads} --gres={cluster.gres}" \
--jobs 2500 \
--latency-wait 120 \
--keep-going \
--local-cores 4 all
