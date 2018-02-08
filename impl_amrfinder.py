import argparse

def run(updater_parser):
    parser = argparse.ArgumentParser(
        parents=[updater_parser],
        description='Run (and update) the amr_finder pipeline.')
    parser.add_argument('-v',   '--version',           dest='version',    help='print version information', action='store_true')
    parser.add_argument('-p',   '--protein',           dest='protein',    help='Amino-acid sequences to search using BLASTP and HMMER')
    parser.add_argument('-n',   '--nucleotide',        dest='nucleotide', help='genomic sequence to search using BLASTX')
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

    args = parser.parse_args()


