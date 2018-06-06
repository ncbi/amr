#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: ncbi/amr:18.05

baseCommand: amr_report
stdout: output.txt
inputs:
  fam:
    type: File
    inputBinding:
      prefix: -fam
  blastp:
    type: File?
    inputBinding:
      prefix: -blastp
  blastx:
    type: File?
    inputBinding:
      prefix: -blastx
  hmmdom:
    type: File?
    inputBinding:
      prefix: -hmmdom
  hmmsearch:
    type: File?
    inputBinding:
      prefix: -hmmsearch
  gff:
    type: File?
    inputBinding:
      prefix: -gff
  outfile:
    type: string?
    #default: "results.sseqid"
    inputBinding:
      prefix: -out
  ident_min:
    type: float?
    inputBinding:
      prefix: -ident_min
  complete_cover_min:
    type: float?
    inputBinding:
      prefix: -complete_cover_min
  partial_cover_min:
    type: float?
    inputBinding:
      prefix: -partial_cover_min
  pseudo:
    type: boolean?
    default: true
    inputBinding:
      prefix: -pseudo
  qc:
    type: boolean?
    default: false
    inputBinding:
      prefix: -qc
  verbose:
    type: int?
    default: 0
    inputBinding:
      prefix: -verbose

outputs:
  output:
    type: stdout
