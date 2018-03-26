#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

#requirements:
    
inputs:
  query: File
  parse_deflines: boolean
  ident_min: float
  cover_min: float
  query_gencode: int
  num_threads: int?
  
outputs:
  result:
    type: File
    outputSource: amr_report/output

steps:
  blastx:
    run: blastx.cwl
    in:
      query: query
      parse_deflines: parse_deflines
      query_gencode: query_gencode
      num_threads: num_threads
    out:
      [output]

  amr_report:
    run: amr_report.cwl
    in:
      blastx: blastx/output
      ident_min: ident_min
      cover_min: cover_min
    out:
      [output]
      
