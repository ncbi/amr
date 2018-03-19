from __future__ import print_function
import argparse
import os
import subprocess
import sys
import tempfile
import yaml
from distutils import spawn

def print_versions(spath):
    os.chdir(os.path.dirname(os.path.realpath(__file__)))
    revision = "Unknown"
    latest = "Unknown"
    try:
        err = open(os.devnull, "wb")
        r = subprocess.check_output("set -o pipefail; svn info | grep ^Revision | cut -d' ' -f2", shell=True, stderr=err)
        revision = r.strip()
        r = subprocess.check_output("set -o pipefail; svn info | grep ^URL | cut -d' ' -f2", shell=True, stderr=err)
        url = r.strip()
        r = subprocess.check_output("set -o pipefail; svn info {} | grep ^Revision | cut -d' ' -f2".format(url), shell=True, stderr=err)
        latest = r.strip()
    except subprocess.CalledProcessError:
        revision = "Unknown"
        latest = "Unknown"
    
    print("  Current CWL Subversion revision:", revision)
    print("   Latest CWL Subversion revision:", latest)
    print("Current Docker container versions:")
    subprocess.check_call("grep -hPo '(?<=dockerPull: )(.*)(?=$)' *.cwl | sort -u | awk '{printf(\"    %s\\n\", $1)}'", shell=True)

class cwlgen:
    def __init__(self, args):
        self.args = args
        self.parse_deflines = True if self.args.parse_deflines else False
        self.do_protein = True if self.args.protein else False
            
    def prot_params(self):
        p = {
            'query': {                                                      
                'class': 'File',
                'location': os.path.realpath(self.args.fasta)
            },
            'parse_deflines': self.parse_deflines
        }
        if self.args.gff:
            p['gff'] = {
                'class': 'File',
                'location': os.path.realpath(self.args.gff)
                }
        return p
    
    def dna_params(self):
        return {
            'query': {                                                      
                'class': 'File',
                'location': os.path.realpath(self.args.fasta)
            },
            'parse_deflines': self.parse_deflines,
            'ident_min': self.args.ident_min,
            'cover_min': self.args.coverage_min,
            'query_gencode': self.args.translation_table
        }
    
    def params(self):
        params = self.prot_params() if self.do_protein else self.dna_params()
        
        (fdstream, self.param_file) = tempfile.mkstemp(suffix=".cwl", prefix="amr_params_")
        stream = os.fdopen(fdstream, 'w')
        yaml.dump(params, stream)
        #print(self.param_file)
        #print(yaml.dump(params))

    def run(self):
        cwlcmd = []
        if spawn.find_executable("cwl-runner") != None:
            cwlcmd = ['cwl-runner']
        elif spawn.find_executable("cwltool") != None:
            cwlcmd = ['cwltool']
        else:
            print("No CWL platform found.", file=sys.stderr)
            sys.exit(1)
        docker_avail = spawn.find_executable("docker")
        if docker_avail == None:
            cwlcmd.extend(['--user-space-docker-cmd', 'udocker'])
        script_path = os.path.dirname(os.path.realpath(__file__))
        script_name = "/wf_amr_prot.cwl" if self.do_protein else "/wf_amr_dna.cwl"
        cwlscript = script_path + script_name
        cwlcmd.extend([cwlscript, self.param_file])
        
        try:
            out = None
            if not self.args.show_output:
                out = open(os.devnull, "wb")
            subprocess.check_call(cwlcmd, stdout=out, stderr=out)

            for line in open('output.txt','r'):
                print(line, end='')
        except subprocess.CalledProcessError:
            print(cwl.stdout)
        except OSError:
            print(cwl.stdout)
        finally:
            self.cleanup()    

    def cleanup(self):
        def safe_remove(f):
            if os.path.exists(f):
                os.remove(f)
        #safe_remove(self.param_file)
        safe_remove("output.txt")
        
        # Cleanup after cwltool's use of py2py3
        safe_remove('/tmp/futurized_code.py')
        safe_remove('/tmp/original_code.py')
        safe_remove('/tmp/py2_detection_code.py')
            
def run(updater_parser):
    parser = argparse.ArgumentParser(
        parents=[updater_parser],
        description='Run (and optionally update) the amr_finder pipeline.')

    parser.add_argument('-d', '--parse_deflines', action='store_true',
                        help='Use -parse_deflines option for blast. This sometimes fixes issues with format of the input FASTA file being automatically parsed by BLAST. (default: %(default)s)')
    parser.add_argument('-o',   '--output', dest='outfile',
                        help='tabfile output to this file instead of STDOUT')

    # Options relating to protein input (-p):
    #parser.add_argument('-f <out.fa> FASTA file containing proteins identified as candidate AMR genes
    parser.add_argument('-g', '--gff', help='GFF file indicating genomic location for proteins')
    # Options relating to nucleotide sequence input (-n)
    parser.add_argument('-i', '--ident_min', type=float,
                        help='Minimum proportion identical translated AA residues (default: %(default)s).')
    parser.add_argument('-c', '--coverage_min', type=float,
                        help='Minimum coverage of reference protein sequence (default: %(default)s).')
    parser.add_argument('-t', '--translation_table', type=int,
                        help='Translation table for blastx (default: %(default)s). More info may be found at https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c')
    parser.set_defaults(ident_min=0.9,
                        coverage_min=0.5,
                        translation_table=11)
    
    parser.add_argument('-s', '--show_output', action='store_true',
                        help='Show the stdout and stderr output from the pipeline execution (verbose mode, useful for debugging).')
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-p', '--protein',    action='store_true',
                       help='Amino-acid sequences to search using BLASTP and HMMER')
    group.add_argument('-n', '--nucleotide', action='store_false',
                       help='genomic sequence to search using BLASTX')
    parser.add_argument('fasta', help='FASTA file containing the query sequence(s).')
    args = parser.parse_args()

    if args.version:
        print_versions()
        sys.exit()
    
    g = cwlgen(args)
    g.params()
    g.run()
    
        
    
