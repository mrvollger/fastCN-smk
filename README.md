# A snakemake for fastCN reference setup

[![Actions Status](https://github.com/mrvollger/fastCN-smk/workflows/CI/badge.svg)](https://github.com/mrvollger/fastCN-smk/actions)
[![Actions Status](https://github.com/mrvollger/fastCN-smk/workflows/Linting/badge.svg)](https://github.com/mrvollger/fastCN-smk/actions)
[![Actions Status](https://github.com/mrvollger/fastCN-smk/workflows/black/badge.svg)](https://github.com/mrvollger/fastCN-smk/actions)
[![DOI](https://zenodo.org/badge/405398596.svg)](https://zenodo.org/badge/latestdoi/405398596)

The `Snakefile` is under `workflow`.


# Install

Please start by installing [pixi](https://pixi.sh/latest/) which handles the environment of this Snakemake workflow.


You can then install the pixi environment by cloning this repository and running:

```
make && pixi install
```

# Usage
pixi handles the execution of the Snakemake workflows:
```
pixi run snakemake ...
```
And if you want to run this Snakemake from another directory you can do so with:
```
pixi run --manifest-path /path/to/snakemake/pixi.toml snakemake ...
```
where you update `/path/to/snakemake/pixi.toml` to the path of the `pixi.toml` you cloned.

And in place of ... use all the normal Snakemake arguments for your workflow.

## Configuration

See `config/config.yaml` for an example layout, and `.test/config.yaml` for a minimal test case.

## Workflow layout

![Workflow](./docs/dag.png)

## TODO

- add raw read depth counts
