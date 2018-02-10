cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: ncbi/hmmer_amr:18.02

baseCommand: hmmsearch
#stdout: $(inputs.db.basename).output
#stdout: $(inputs.db).output
inputs:
  tblout:
    type: string?
    default: hmmsearch.out
    inputBinding:
      position: 1
      prefix: --tblout
  noali:
    type: boolean?
    default: true
    inputBinding:
      position: 2
      prefix: --noali
  domtblout:
    type: string?
    default: domtbl.out
    inputBinding:
      position: 3
      prefix: --domtblout
  cut_tc:
    type: boolean?
    default: true
    inputBinding:
      position: 4
      prefix: --cut_tc
  Z:
    type: int?
    default: 10000
    inputBinding:
      position: 5
      prefix: -Z
  db:
    type: string?
    default: /hmmer/hmmerdb/AMR.LIB
    inputBinding:
      position: 6
#      prefix: 
#      separate: false
  query:
    type: File
    inputBinding:
      position: 7

outputs:
  - id: hmmsearch_out
    type: File
    outputBinding:
      glob: "hmmsearch.out"
  - id: hmmdom_out
    type: File
    outputBinding:
      glob: "domtbl.out"
