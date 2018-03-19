cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: ncbi/amr_utils:18.02

baseCommand: amr_report

stdout: output.txt

inputs:
  fam:
    type: string?
    default: "/fam.tab"
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
  cover_min:
    type: float?
    inputBinding:
      prefix: -cover_min
  pseudo:
    type: boolean?
    default: false
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

