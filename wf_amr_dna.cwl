#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
    
inputs:
  query: File
  fasta: File
  fam: File
  parse_deflines: boolean
  ident_min: float
  complete_cover_min: float
  query_gencode: int
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

  blastx:
    run: blastx.cwl
    in:
      query: query
      db: makeblastdb/blastdb
      parse_deflines: parse_deflines
      query_gencode: query_gencode
      num_threads: num_threads
    out:
      [output]

  amr_report:
    run: amr_report.cwl
    in:
      fam: fam
      blastx: blastx/output
      ident_min: ident_min
      complete_cover_min: complete_cover_min
    out:
      [output]
      
