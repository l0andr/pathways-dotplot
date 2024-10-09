# Pathway-Dotplot

Pathway-Dotplots is a Python-based tool for generating dot plots from pathway analysis results, with hierarchical clustering support and several plot customization options. The tool is primarily intended for visualizing the enrichment of biological pathways across different samples using data such as p-values, normalized enrichment scores (NES), and more.

## Features

- **Hierarchical clustering**: Cluster pathways based on similarity in enrichment scores.
- **Customizable Dotplots**: Multiple plot styles including size-based, color-based, gray-scale, and more.
- **Input flexibility**: Accepts multiple pathway analysis result files in CSV format.
- **P-value thresholds**: Highlight pathways that pass a significance threshold.
- **High-quality output**: Generates publication-ready PDF plots.

## Requirements

- Python 3.7+
- Required packages: `matplotlib`, `numpy`, `pandas`, `scipy`

You can install the required packages by running:
```bash
pip install -r requirements.txt
