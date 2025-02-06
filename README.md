# bbi-scirna-demux

## Intro

This *bbi-scirna-demux* pipeline runs the Illumina *bcl-convert* program to make *.bam* files that are demultiplexed by PCR primer pairs, followed by the *rna_rtlig_demux* program to make sample-specific pcr primer pair *.bam* files.

## Installation

Install the following software

- Nextflow: this pipeline uses Nextflow DSL2 so you must install a recent version of Nextflow. I use version 24.10.2 successfully. If you need to run the *bbi-dmux* and *bbi-sci* pipelines too, install the new Nextflow in a different location because the recent versions no longer support DSL1.
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

This is a bit lengthy at this time. The steps are

- run *bbi-scirna-demux/samplesheet/samplesheet_scrunch.py* on a samplesheet CSV file that is suitable for the *bbi-dmux* pipeline. The input file has one row per RT well and the output CSV file has one row per sample. *samplesheet_scrunch.py* also adds columns that give the PCR primer wells or columns and rows. Run *samplesheet_scrunch.py --help* command for more information.
- run *bbi-scirna-demux/samplesheet/scirna_samplesheet.py* to convert the scrunched CSV file to a JSON file. *scirna_samplesheet.py* requires a command line parameter that gives the number of lanes used in the sequencing run. Run *scirna_samplesheet.py -d* for detailed documentation. (At this early stage of the program's life, there may be errors and omissions in the documentation.) You may need to edit the scrunched CSV file in a spreadsheet program in order to add columns described in the *scirna_samplesheet.py* documentation.

### Run *bbi-scirna-demux*

Use the *run.demux.sh* bash script to start the pipeline run.

The output *.bam* files are in the directory demux_out. The *.bam* file contains unaligned read sequences and quality values as well as barcode data in the SAMtags.
