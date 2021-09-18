#!/usr/bin/env bash
set -euo pipefail

mkdir dir -p logs/drmaa

snakemake \
    --drmaa " -l centos=7 -l h_rt=48:00:00 -l mfree={resources.mem}G -pe serial {threads} -V -cwd -S /bin/bash -w n" \
    --drmaa-log-dir logs/drmaa \
    --use-conda \
    "$@"
