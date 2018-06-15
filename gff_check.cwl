#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: ncbi/amr:18.06

# Check the correctness of a .gff-file. Exit with an error if it is incorrect.
# Usage: gff_check <gff> [-qc] [-verbose 0] [-noprogress] [-profile] [-threads 1] [-json ""] [-log ""] [-fasta ""] [-locus_tag ""]
# Help:  gff_check -help|-h
# Parameters:
#   <gff>: .gff-file, if an empty string then exit 0
#   [-qc]: Integrity checks (quality control)
#   [-verbose 0]: Level of verbosity
#   [-noprogress]: Turn off progress printout
#   [-profile]: Use chronometers to profile
#   [-threads 1]: Max. number of threads
#   [-json ""]: Output file in Json format
#   [-log ""]: Error log file, appended
#   [-fasta ""]: Protein FASTA file
#   [-locus_tag ""]: File with matches: "<FASTA id> <GFF id>", where <id> is from "[locus_tag=<id>]" in the FASTA comment and from the .gff-file

baseCommand: gff_check
stdout: gff_check.out
inputs:
  gff:
    type: File?
    default:
      class: File
      basename: "emptystring"
      contents: ""
    inputBinding:
      position: 1
  gff_file:
    type: File?
  qc:
    type: string?
    inputBinding:
      prefix: -qc
      position: 2
  verbose:
    type: int?
    inputBinding:
      prefix: -verbose
      position: 3
  noprogress:
    type: boolean?
    inputBinding:
      prefix: -noprogress
      position: 4
  profile:
    type: boolean?
    inputBinding:
      prefix: -profile
      position: 5
  threads:
    type: int?
    inputBinding:
      prefix: -threads
      position: 6
  json:
    type: string?
    inputBinding:
      prefix: -json
      position: 7
  log:
    type: string?
    inputBinding:
      prefix: -log
      position: 8
  fasta:
    type: File?
    default: ""
    inputBinding:
      prefix: -fasta
      position: 9
  locus_tag:
    type: File?
    inputBinding:
      prefix: -locus_tag
      position: 10
      
outputs:
  - id: output
    type: File
    outputBinding:
      glob: "gff_check.out"
