# twist_myeloid_lymphoid
Repository for the analysis pipeline on the Twist Myeloid-lymphoid panel

### What is in there:
  - shell scripts to download raw fastq and bcl files
  - shell script to convert bcl files to fastq files if needed
  - individual shell and R scripts for running the analysis
  - A complete Snakemake workflow

## Preparing files
To download fastq files, a shell script is already prepared to make it easier. 


## Snakemake workflow
![description](code/snakefile/rulegraph.png)
## Running snakemake
Depending on what you want to run snakemake on, there are different options built in.
The workflow ignores undetermined and control fastq files.
### Running all rules for all files in a run
When running this workflow, certain things need to be specified in the command line:
```
snakemake --cores 4 --directory /base/working/directory --config run=run_number --configfile config.yaml
```

- cores: the number of cores specified for the run. Since it runs some jobs parallel, more cores results in less runtime, but takes up more resources. Recommended: 4
- directory: working directory for analysis.
- config: specifying analysis run number (e.g. run=run_30)
- configfile: .yaml file containing necessary files for analysis. In case a file needs to be replaced (for example a bed file), it should be changed in this file.

