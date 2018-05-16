#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: ncbi/amr:18.05

# Check the correctness of a FASTA file. Exit with an error if it is incorrect.
# Usage: fasta_check <in> [-qc] [-verbose 0] [-noprogress] [-profile] [-json ""] [-log ""] [-aa] [-hyphen]
# Help:  fasta_check -help|-h
# Parameters:
#   <in>: FASTA file
#   [-qc]: Integrity checks (quality control)
#   [-verbose 0]: Level of verbosity
#   [-noprogress]: Turn off progress printout
#   [-profile]: Use chronometers to profile
#   [-json ""]: Output file in Json format
#   [-log ""]: Error log file, appended
#   [-aa]: Amino acid sequenes, otherwise nucleotide
#   [-hyphen]: Hyphens are allowed
    
baseCommand: fasta_check
stdout: fasta_check.out
inputs:
  fasta:
    type: File
    inputBinding:
      position: 1
  qc:
    type: string?
    inputBinding:
      prefix: -qc
  verbose:
    type: int?
    inputBinding:
      prefix: -verbose
  noprogress:
    type: boolean?
    inputBinding:
      prefix: -noprogress
  profile:
    type: boolean?
    inputBinding:
      prefix: -profile
  json:
    type: string?
    inputBinding:
      prefix: -json
  log:
    type: string?
    inputBinding:
      prefix: -log
  aa:
    type: boolean?
    inputBinding:
      prefix: -aa
  hyphen:
    type: boolean?
    inputBinding:
      prefix: -hyphen
      
outputs:
  - id: output
    type: File
    outputBinding:
      glob: "fasta_check.out"
