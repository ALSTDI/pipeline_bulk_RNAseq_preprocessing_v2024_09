#!/bin/bash
# Use this script to archive the logs from the RNAseq bulk nextflow pipeline to the appropriate 
# project folder on the nextflow-batch-output bucket 
gsutil cp 2501_ARC.yaml gs://nextflow-batch-output/bulkRNA/2501_ARC/_project_logs/
gsutil cp report-20250121-75320685.html gs://nextflow-batch-output/bulkRNA/2501_ARC/_project_logs/
gsutil cp sampleSheet.csv gs://nextflow-batch-output/bulkRNA/2501_ARC/_project_logs/