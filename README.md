# Detect horizontal gene transfer events from metagenome data

Zeqian Li

Last updated: May 21st, 2023

Snakemake implementation of [Anvi'o's gene clustering](https://merenlab.org/2016/11/08/pangenomics-v2/) for pan-genome analysis.

1. From a fasta file, create a all-to-all blasta/MMSeq2. 
2. Format the file as an mcl input. 
3. Cluste by mcl.
4. Format the output.

## Usage 

Install the following conda environments
- TODO

Run:
- Change configurations in `config/config.yml`. See commments in the file. 
- Dry run for testing: 
    ``` snakemake --cores all -n ```
- Actual run
``` snakemake --cores all --conda-frontend conda --use-conda -k```
- Use a custom config file: 
``` snakemake --cores all --configfile config/config_sag.yml --conda-frontend conda --use-conda -k```
