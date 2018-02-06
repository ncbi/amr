cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: ncbi/amr:18.02

baseCommand: amr_report

inputs:
  fam:
    type: string?
    default: "/fam.tab"
    inputBinding:
      prefix: -fam
  blastp:
    type: File
    inputBinding:
      prefix: -blastp
  hmmdom:
    type: File
    inputBinding:
      prefix: -hmmdom
  hmmsearch:
    type: File
    inputBinding:
      prefix: -hmmsearch
  outfile:
    type: string
    default: "results.sseqid"
    inputBinding:
      prefix: -out
  verbose:
    type: int?
    default: 0
    inputBinding:
      prefix: -verbose

outputs:
  - id: output
    type: File
    outputBinding:
      glob: "*.sseqid"

