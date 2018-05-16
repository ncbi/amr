#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: DockerRequirement
    dockerPull: ncbi/amr:18.05

  
inputs:
  query: File
  fasta: File
  hmmdb: File
  fam: File
  gff: File?
  parse_deflines: boolean
  num_threads: int?
  cpu: int?
  
outputs:
  result:
    type: File
    outputSource: amr_report/output

steps:
  makeblastdb:
    run: wf_makeblastdb.cwl
    in:
      fasta: fasta
    out:
      [blastdb]
  
  blastp:
    run: blastp.cwl
    in:
      query: query
      db: makeblastdb/blastdb
      parse_deflines: parse_deflines
      num_threads: num_threads
    out:
      [output]

  hmmsearch:
    run: hmmsearch.cwl
    in:
      query: query
      db: hmmdb
      cpu: cpu
    out:
      [hmmsearch_out,hmmdom_out]

  amr_report:
    run: amr_report.cwl
    in:
      fam: fam
      blastp: blastp/output
      hmmdom: hmmsearch/hmmdom_out
      hmmsearch: hmmsearch/hmmsearch_out
      gff: gff
    out:
      [output]
      
