# Pathway-Dotplot

Pathway-Dotplots is a Python-based tool for generating dot plots from pathway analysis results, with hierarchical clustering support and several plot customization options. The tool is primarily intended for visualizing the enrichment of biological pathways across different samples using data such as p-values, normalized enrichment scores (NES), and more.

![image](https://github.com/user-attachments/assets/bcfa4daa-32d0-476d-a32f-14b5a435465c)



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
```

## Usage
  The tool processes input CSV files containing pathway enrichment results and generates dot plots based on the selected options.

## Example command:
```bash
python pathway_dotplots.py -indir /path/to/input/dir -outdir /path/to/output/dir --plot_type size --pvalue_threshold 0.05 --show  
```
## Command-line options:
-indir: Path to the directory containing input files (in CSV format). <br>
--input_file_mask: Mask for the input files (default: *.csv). <br>
-outdir: Path to the directory where output plots will be saved. <br>
--plot_type: Type of plot to generate. <br>
Options: <br> 
  * size: Dot size will be proportional to -log10(pvalue).  
  * invsize: Dot size will be proportional to the enrichment score.  
  * gray: Dots with p-value greater than the threshold will be plotted in gray.  
  * white: Dots with p-values below the threshold will be plotted in white.

--pathway_sort: Pathway sorting strategy. <br>
Options: <br>
* cluster: Sort pathways using hierarchical clustering.  
* first_enrich: Sort pathways based on the first sample's enrichment score.  
* enrich_threshold: Sort pathways using enrichment score and p-value threshold.  
--pvalue_threshold: Threshold for p-value (default: 0.1). <br>
--show: If specified, the plot will be displayed in the GUI. <br>
--imgname: Filename for the output plot (default: dot_plot). <br>
Input Format
Each input file should be a CSV containing pathway enrichment results for a specific sample. The following columns are required: <br>

pathway: The name of the pathway. <br>
NES: Normalized Enrichment Score (or any other enrichment score).<br>
padj: Adjusted p-value.<br>

Output
The tool generates a dot plot in PDF format for each sample, displaying pathways across samples, where the size or color of the dots reflects pathway enrichment or p-value.<br>

