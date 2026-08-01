# metagenomics-diversity-toolkit

A toolkit for computing alpha/beta diversity metrics and statistical analysis of metagenomic community composition data.

## Overview

This repository provides a command-line tool for computing standard ecological diversity metrics from an OTU or ASV abundance table: Shannon and Simpson alpha diversity per sample, and pairwise Bray-Curtis beta diversity dissimilarity between samples.

## Installation

```
pip install -r requirements.txt
```

## Usage

```
python scripts/diversity_metrics.py --table examples/example_otu_table.tsv --alpha-output alpha_diversity.tsv --beta-output beta_diversity.tsv
```

## Input format

The abundance table must be tab-separated with an OTU_ID column followed by one column per sample containing read counts or relative abundances.

## License

This project is licensed under the MIT License. See the LICENSE file for details.
