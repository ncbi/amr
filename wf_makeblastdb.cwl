#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: StepInputExpressionRequirement
    
inputs:
  fasta: File
  
outputs:
  blastdb:
    type: Directory
    outputSource: mkdir/blastdb

steps:
  makeblastdb:
    run:
      class: CommandLineTool
      hints:
        DockerRequirement:
          dockerPull: ncbi/amr:latest
      requirements:
        - class: InitialWorkDirRequirement
          listing:
            - entry: $(inputs.fasta)
              writable: False
      #makeblastdb -in AMRProt -dbtype prot
      baseCommand: makeblastdb
      inputs:
        fasta:
          type: File
          inputBinding:
            prefix: -in
        dbtype:
          type: string?
          default: prot
          inputBinding:
            prefix: -dbtype

      outputs:
        blastfiles:
          type: File[]
          outputBinding:
            glob: "*"
    in:
      fasta: fasta
    out:
      [blastfiles]

  mkdir:
    run:
      class: CommandLineTool
      requirements:
        - class: ShellCommandRequirement
      arguments:
        - shellQuote: false
          valueFrom: >-
            mkdir $(inputs.blastdir) && cp
      inputs:
        blastfiles:
          type: File[]
          inputBinding:
            position: 1
        blastdir:
          type: string
          inputBinding:
            position: 2

      outputs:
        blastdb:
          type: Directory
          outputBinding:
            glob: $(inputs.blastdir)
    in:
      blastfiles: makeblastdb/blastfiles
      blastdir:
        source: fasta
        valueFrom: $(self.basename)
    out:
      [blastdb]


