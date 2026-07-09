# bbi-scirna-demux

## Intro

This *bbi-scirna-demux* pipeline runs the Illumina *bcl-convert* program to make *.bam* files that are demultiplexed by PCR primer pairs, followed by the *rna_rtlig_demux* program to make sample-specific pcr primer pair *.bam* files. Use the *main* branch.

## Installation

Install the following software

- Nextflow: this pipeline uses Nextflow DSL2 so you must install a recent version of Nextflow. I use version 24.10.2 successfully. If you need to run the *bbi-dmux* and *bbi-sci* pipelines too, you will need two different Nextflow version so install the new Nextflow in its own location because the recent versions no longer support DSL1.
- Rust: the *rna_rtlig_demux* program is written in Rust so you must install the Rust compiler.
- rna_rtlig_demux: this program is a compiled program written in Rust so you must compile it and copy the executable to the *bbi-scirna-demux/bin* directory.
- python3 interpreter: I use version 3.12.1 successfully.

### Install Nextflow

See the Nextflow installation instructions at

https://www.nextflow.io/docs/latest/install.html

### Install Rust

See the Rust installation instructions at

https://www.rust-lang.org/tools/install

### Build and install *rna_rtlig_demux*

Run the following commands

```
cd bbi-scirna-demux/src/rna_rtlig_demux
cargo build --release
cp target/release/rna_rtlig_demux ../../bin
```

I recommend that you build *rna_rtlig_demux* on a newer cluster node, for example, s020 on the Shendure cluster.

### Load Python3 on the Genome Sciences cluster

The Python3 interpreter is loaded by Nextflow. The module load command is in the file *bbi-scirna-demux/nextflow.config*.

## Run *bbi-scirna-demux*

I recommend that you use the *run.demux.sh* script in this repository. You must edit the script to use the correct Nextflow program and Nextflow main.nf script.

### Edit the *experiment.config* file.

Edit *experiment.config* to set the following parameters for your run:

```
params.samplesheet_json
params.illumina_run_dir
params.output_dir
```

*params.samplesheet_json* is the path to your samplesheet JSON file for the run. *params.illumina_run_dir* is the path to the Illumina flowcell directory, and *params.output_dir* is the path to the directory where the processing output is written to.

### Make the samplesheet JSON file.

The file bbi-scirna-demux/samplesheet/scirna_samplesheet.py has detailed information for making a samplesheet JSON file. The steps are

- run *bbi-scirna-demux/samplesheet/lims2scrunch.py* on a LIMS CSV manifest file to make a CSV samplesheet file where each row describes a sample. Run *lims2scrunch.py --help* command for more information.
- alternatively, run *bbi-scirna-demux/samplesheet/samplesheet_scrunch.py* on a samplesheet CSV file that is suitable for the *bbi-dmux* pipeline. In the simplest case, the input file has one row per RT well and the output CSV file has one row per sample. *samplesheet_scrunch.py* also adds columns that give the PCR primer wells or columns and rows. Run *samplesheet_scrunch.py --help* command for more information.
- edit the scrunched *.csv* file as needed using a spreadsheet program. This is needed to add the hash barcode file paths, for example. See the scirna_samplesheet.py* program for information about the input *.csv* file format.
- run *bbi-scirna-demux/samplesheet/scirna_samplesheet.py* to convert the scrunched, edited  CSV file to a JSON file. *scirna_samplesheet.py* requires a command line parameter that gives the number of lanes used in the sequencing run. Run *scirna_samplesheet.py -d* for detailed documentation. (At this early stage of the program's life, there may be omissions and errors in the documentation.)
- the ligation barcode file for the standard and jumbo/mega sci experiments differ. Use the *bbi-scirna-demux/data/ligation.txt* for standard sci and *bbi-scirna-demux/data/ligation_megasci.row_sorted.tsv* for jumbo sci.

### Run *bbi-scirna-demux*

Use the *run.demux.sh* bash script to start the pipeline run.

The output *.bam* files are in the directory demux_out. The *.bam* file contains unaligned read sequences and quality values as well as barcode data in the SAMtags.

## Overview of the demux and analyze pipeline workflows.

- bcl-convert makes fastq files demultiplexed by PCR barcode pairs. There is one fastq file for each acceptable PCR barcode pair (and the fastq filenames have the p7 and p5 barcode sequence indices, as well as the lane number, in them). These fastq files are used internally only, and are not returned to the user, i.e., 'published' by the Nextflow pipeline.
- the reads in these fastq files are demultiplexed by RT and ligation barcodes. The RT barcodes identify reads by sample. The resulting reads are written to unaligned BAM files where all reads in a BAM file have the same lane and rt, p7, and p5 barcodes. The barcode sequence indices and lane number are part of the BAM file name. These BAM files are 'published' to the 'demux_out' directory.
- these (sample, pcr pair, lane) BAM files are merged such that reads in the output BAM file belong to the same sample, process_group, and have the same p7 and p5 barcodes. This means that BAM files whose reads have the same sample_name, process_group, p7 barcode, and p5 barcode, but different lanes, are merged. If you have the same sample, with the same sample_name, in multiple libraries, each library run in a different lane, and you need the pipeline to process separately the reads from each library, assign a distinct process_group value to the sample entries for each library. These merged BAM files are not 'published'.
- the unaligned BAM files are processed by trimgalore to trim off adapter sequence. These trimmed BAM files are not published.
- the trimmed read BAM files are aligned to the reference genome by STARsolo and the aligned reads are assigned to cells. STARsolo output files are for reads that have the same sample_name, process_group, and p7 and p5 barcodes. These aligned read BAM files are not 'published'.
- these (sample, process_group, pcr pair) aligned output files are merged such that all reads with the same sample_name and process_group, are written to the same output BAM file. This means that BAM files with the same sample_name and process_group, but different p7 and/or p5 barcodes, are merged. If you have samples in a lane that have the same sample_name but are from different libraries, where these samples are distinguished using different p7 and p5 barcode combinations, you must assign the samples distinct process_group values to avoid merging their reads by this stage. The merged BAM, count matrix, and statistics files are 'published' to the 'analyze_out' directory.
- a cds file and umap.png file are made for each sample, process_group count matrix.
- when the hash_file value is set for a sample, the untrimmed BAM files are processed to find candidate hash reads and a cds is made with the hash read information.

## Notes

- the samplesheet.json file includes a list of lanes to which each sample is applied. If you make a samplesheet.json file for a machine with one lane and later run the library on a machine with a different number of lanes, you must regenerate the samplesheet.json file.


