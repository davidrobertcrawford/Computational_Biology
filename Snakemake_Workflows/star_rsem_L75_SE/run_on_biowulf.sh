#!/bin/bash

sbatch --time=24:00:00 --mem=4g --partition=norm,ccr run_snakemake.sh
