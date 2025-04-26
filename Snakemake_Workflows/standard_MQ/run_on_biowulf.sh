#!/bin/bash

sbatch \
--time=2:00:00 \
--mem=4g \
--partition=norm,ccr \
run-snakemake.sh
