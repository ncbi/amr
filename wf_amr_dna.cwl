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
  fasta_check_out:
    type: File
    outputSource: fasta_check/output
    
steps:
  fasta_check:
    run: fasta_check.cwl
    in:
      fasta: query
    out:
      [output]
  
  makeblastdb:
    run: wf_makeblastdb.cwl
    in:
      fasta_check_dummy: fasta_check/output
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
      
