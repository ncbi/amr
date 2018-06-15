#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: ncbi/amr:18.06

baseCommand: blastp
stdout: blastp.out
inputs:
  query:
    type: File
    inputBinding:
      prefix: -query
  db:
    type: Directory
    inputBinding:
      prefix: -db
      valueFrom: $(self.path)/$(self.basename)
  task:
    type: string?
    default: blastp-fast
    inputBinding:
      prefix: -task
  outfmt:
    type: string?
    default: "6 qseqid sseqid length nident qstart qend qlen sstart send slen qseq"
    inputBinding:
      prefix: -outfmt
  show_gis:
    type: boolean
    default: true
    inputBinding:
      prefix: -show_gis
  word_size:
    type: int?
    default: 6
    inputBinding:
      prefix: -word_size
  threshold:
    type: int?
    default: 21
    inputBinding:
      prefix: -threshold
  evalue:
    type: double?
    default: 1e-20
    inputBinding:
      prefix: -evalue
  comp_based_stats:
    type: int?
    default: 0
    inputBinding:
      prefix: -comp_based_stats
  num_threads:
    type: int?
    inputBinding:
      prefix: -num_threads
  parse_deflines:
    type: boolean?
    default: true
    inputBinding:
      prefix: -parse_deflines

outputs:
  - id: output
    type: File
    outputBinding:
      glob: "blastp.out"
