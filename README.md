# CNV analysis with Plotly

This repository contains scripts for analyzing and visualizing CNV (Copy Number Variation) data using Plotly.
The scripts are designed to work with CNVkit ratio files, containing the following columns:
- chromosome
- start
- end
- gene
- depth
- log2
- weight

## Installation
To set up and quickly start working with this repository is via Conda, which ensures installation of all dependencies. First clone the repository and then create conda environment with the following command:
```bash
git clone https://github.com/altuda/CNVplotter.git
cd CNVplotter
conda env create -f environment.yml
```


## Scripts
### `run_plotly_cnv.py`

This script is designed to analyze and visualize CNV data using Plotly. It reads CNVkit ratio files, calculates rolling means, and generates interactive plots with customizable features.

#### Usage

```bash
python run_plotly_cnv.py -i <input_directory> -c <chromosome_sizes_file> -o <output_directory> -p <prefix> -g <additional_genes>
```

## Argument Descriptions

- `-i` or `--input_dir`: **Required**. Specify the directory containing the CNV files.
- `-c` or `--chrom_sizes`: **Required**. Provide the path to the chromosome sizes file.
- `-o` or `--output_dir`: **Optional**. Set the directory for output files (default is `output`).
- `-p` or `--prefix`: **Optional**. Define a prefix for output file names (default is `CNV_plot`).
- `-g` or `--genes`: **Optional**. Provide a comma-separated list of additional genes to highlight in the plots.



