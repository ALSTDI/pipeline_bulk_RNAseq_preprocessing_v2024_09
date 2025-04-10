# pipeline_bulk_RNAseq_preprocessing_v2024_09

Pipeline to preprocess RNAseq fastqs to reads by genes or isoforms for secondary analyses. Pipeline automatically detects whether the reads are single-end or paired-end, and the strandedness of the reads.

# Quick start: 
To run the pipeline:  <br>
- Modify `sampleSheet.csv` <br>
- Create `YYMM_PROJECT.yaml` based on `example.yaml` and modify accordingly. We will keep a log of `sampleSheet.csv` and this `MMYY_PROJECT.yaml` <br>
- Initiate the pipeline on GCloud. You will not need the `local` profile unless you're debugging.  <br>
```
nextflow run -profile cloud main.nf -params-file MMYY_PROJECT.yaml -with-report
```
- After the pipeline is finished, edit the project location from `archiving_logs.sh` and push related files to PROJECT/_project_logs for book-keeping. 

# Pipeline file structure:
. <br>
├── README.md <br>
├── example.yaml <br>
├── flowchart.html <br>
├── main.nf <br>
├── modules <br>
│   ├── diag <br>
│   │   └── main.nf <br>
│   ├── paired_processes <br>
│   │   └── main.nf <br>
│   └── single_processes <br>
│       └── main.nf <br>
├── nextflow.config <br>
└── sampleSheet.csv <br>

## Static files:
This repository contains the following static components that you will not have to modify:  <br>
1. `README.md`: this file <br>
2. `flowchart.html`: workflow diagram for steps in this pipeline <br>
3. `main.nf`: the main executable that initiates the bulk RNAseq nextflow preprocessing workflow  <br>
4. `modules`: modules store processes that are used in main.nf workflow  <br>
5. `nextflow.config`: you usually will not have to modify this file unless you are an advanced user <br>

## Files need modifying:
1. `example.yaml`: Environmental variables are supplied in this example yaml file to run the pipeline. The example here points to the human genome for STAR, rsem index in `refDir`, and the salmon index in `salmonDir`. The example also sets the project directory `projectDir` to the input bucket - where all input fastqs live, and the output location `out_bucket` - where all pipeline results will be stored.  <br>
2. `sampleSheet.csv`: Supply the sample information for this project. You can reference this example `sampleSheet.csv` here to modify this file accordingly. 

# CHANGELOG: 
#### February 2025: 
Add publishDir for STAR aligned bam outputs 

# TODO:
Hard-coded variables that may need to be dynamic in the future: <br>
diag process is grabbing ${SAMPLE}_2.fq.gz - this is the deliv format from Novogene. (?) <br>
Currently, the pipeline is static in terms of CPUs and memory and bootdisk size needed for each process. We need to try to optimize this based on task.attempt at a later time point. Especially when bootdisk size is set to 500Gb now...  <br>

Add a clean up step at the end. <br>



