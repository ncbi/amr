cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: ncbi/blast_amr:18.02

baseCommand: blastx
#stdout: $(inputs.db).out
stdout: blastx.out
inputs:
  query:
    type: File
    inputBinding:
      prefix: -query
  db:
    type: string
    default: AMRProt
    inputBinding:
      prefix: -db
  query_gencode:
    type: int?
    default: 11
    inputBinding:
      prefix: -query_gencode
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
    default: 3
    inputBinding:
      prefix: -word_size
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
  max_target_seqs:
    type: int?
    default: 10000
    inputBinding:
      prefix: -max_target_seqs
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
      glob: "blastx.out"
