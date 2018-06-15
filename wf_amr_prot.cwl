#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: DockerRequirement
    dockerPull: ncbi/amr:18.06

  
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
  fasta_check_out:
    type: File
    outputSource: fasta_check/output
  gff_check_out:
    type: File
    outputSource: gff_check/output

steps:
  fasta_check:
    run: fasta_check.cwl
    in:
      fasta: query
      aa:
        default: True
    out:
      [output]

  gff_check:
    run: gff_check.cwl
    in:
      gff: gff
      fasta: query
      #locus_tag:
      #  default: locus.tags
    out:
      [output]
      
  makeblastdb:
    run: wf_makeblastdb.cwl
    in:
      fasta_check_dummy: fasta_check/output
      gff_check_dummy: gff_check/output
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
      fasta_check_dummy: fasta_check/output
      gff_check_dummy: gff_check/output
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
      
