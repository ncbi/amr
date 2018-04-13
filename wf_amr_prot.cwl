#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  
inputs:
  query: File
  fasta: File
  gff: File?
  parse_deflines: boolean
  num_threads: int?
  
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
    out:
      [hmmsearch_out,hmmdom_out]

  amr_report:
    run: amr_report.cwl
    in:
      blastp: blastp/output
      hmmdom: hmmsearch/hmmdom_out
      hmmsearch: hmmsearch/hmmsearch_out
      gff: gff
    out:
      [output]
      
