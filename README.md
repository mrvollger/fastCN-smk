# A snakemake for fastCN reference setup

[![Actions Status](https://github.com/mrvollger/fastCN-smk/workflows/CI/badge.svg)](https://github.com/mrvollger/fastCN-smk/actions)
[![Actions Status](https://github.com/mrvollger/fastCN-smk/workflows/Linting/badge.svg)](https://github.com/mrvollger/fastCN-smk/actions)
[![Actions Status](https://github.com/mrvollger/fastCN-smk/workflows/black/badge.svg)](https://github.com/mrvollger/fastCN-smk/actions)
[![DOI](https://zenodo.org/badge/405398596.svg)](https://zenodo.org/badge/latestdoi/405398596)

The `Snakefile` is under `workflow`.

## Installing

Most installing is done by `Snakemake` but there is one small Makefile. Please refer to it so you can properly set your LD path before using.

## Configuration

See `config/config.yaml` for an example layout, and `.test/config.yaml` for a minimal test case.

## Workflow layout

![Workflow](./docs/dag.png)
