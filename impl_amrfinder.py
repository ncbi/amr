import argparse
import os
import subprocess
import yaml

def run(updater_parser):
    parser = argparse.ArgumentParser(
        parents=[updater_parser],
        description='Run (and update) the amr_finder pipeline.')
    parser.add_argument('-v', '--version', action='store_true', help='print version information')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-p', '--protein',    action='store_true',   help='Amino-acid sequences to search using BLASTP and HMMER')
    group.add_argument('-n', '--nucleotide', action='store_false',  help='genomic sequence to search using BLASTX')

    parser.add_argument('-npd', '--no_parse_deflines', action='store_true',
                        help='Do not use -parse_deflines option for blast (sometimes fixes issues with format of the input FASTA file being automatically parsed by BLAST)')
    parser.add_argument('-o',   '--output',            dest='outfile',    help='tabfile output to this file instead of STDOUT')

    # Options relating to protein input (-p):
    #parser.add_argument('-f <out.fa> FASTA file containing proteins identified as candidate AMR genes
    #parser.add_argument('-g <gff> GFF file indicating genomic location for proteins in -p <protein>
    # Options relating to nucleotide sequence input (-n)
    #parser.add_argument('-i <0.9> Minimum proportion identical translated AA residues
    #parser.add_argument('-c <0.5> Minimum coverage of reference protein sequence
    #parser.add_argument('-t <11> Translation table for blastx

    parser.add_argument('fasta', help='FASTA file containing the query sequence(s).')
    args = parser.parse_args()

    parse_deflines = True
    if args.no_parse_deflines:
        parse_deflines = False
    
    params =  {
        'query': {                                                      
            'class': 'File',                                               
            'location': args.fasta
        },
        'parse_deflines': parse_deflines
    }
    param_file = 'armfinder_params.yaml'
    stream = open(param_file, 'w')
    yaml.dump(params, stream)
    #print(yaml.dump(params))

    script_path = os.path.dirname(os.path.realpath(__file__))
    cwlscript = script_path + "/wf_amr.cwl"
    cwl = subprocess.run(['cwltool', cwlscript, param_file],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)

    for line in open('results.sseqid','r'):
        print(line, end='')
    
