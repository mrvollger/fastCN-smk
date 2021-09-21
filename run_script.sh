#!/usr/bin/env bash
set -euo pipefail

mkdir dir -p logs/drmaa

snakemake \
    --drmaa " -l centos=7 -l h_rt=48:00:00 -l mfree={resources.mem}G -pe serial {threads} -V -cwd -S /bin/bash -w n" \
    --drmaa-log-dir logs/drmaa \
    --use-conda \
    --configfile config/config.yaml \
    --local-cores 20 \
    --cores 20 \
    --max-inventory-time 10000 \
    --resources load=1000 \
    --scheduler greedy \
    --latency-wait 60 \
    --restart-times 3 \
    "$@"
