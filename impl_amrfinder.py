from __future__ import print_function
import argparse
import errno
import os
#import re
import subprocess
import sys
import tempfile
import yaml
from distutils import spawn
from ftplib import FTP

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

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def available_cpu_count():
    try:
        import multiprocessing
        return multiprocessing.cpu_count()
    except (ImportError, NotImplementedError):
        pass

    return 1
        
def check_data(files = [ 'AMR.LIB', 'AMRProt', 'fam.tab' ]):
    pre = os.path.dirname(os.path.realpath(__file__)) + "/data/latest/"
    if not os.path.isdir(pre):
        return False
    for base in files:
        f = pre + base
        if not os.path.isfile(f):
            return False
    return True

def get_latest(ls):
    # Example results to be parsed
    #'modify=20180330183106;perm=fle;size=4096;type=dir;unique=28UFC8F9;UNIX.group=562;UNIX.mode=0444;UNIX.owner=14; 2018-03-30.1'
    #'modify=20180410155844;perm=adfr;size=12;type=OS.unix=symlink;unique=28UFC8FE;UNIX.group=562;UNIX.mode=0444;UNIX.owner=14; latest'
    #'modify=20180409170726;perm=fle;size=4096;type=dir;unique=28UFC8FE;UNIX.group=562;UNIX.mode=0444;UNIX.owner=14; 2018-04-09.1'
    dirs = {}
    latest = ''
    for f in ls:
        fields = f.split(';')
        typ = fields[3].split('=', 1)[1]
        if typ == 'dir':
            u = fields[4].split('=', 1)[1]
            n = fields[8].strip()
            dirs[u] = n
        elif typ == 'OS.unix=symlink':
            n = fields[8].strip()
            if n == 'latest':
                latest_uniq = fields[4].split('=', 1)[1]

    return dirs[latest_uniq]
    
def update_data():
    prevdir = os.getcwd()
    pre = os.path.dirname(os.path.realpath(__file__)) + "/data/"
    mkdir_p(pre)
    os.chdir(pre)
    ftp = FTP('ftp.ncbi.nlm.nih.gov')     # connect to host, default port
    ftp.login()                     # user anonymous, passwd anonymous@
    ftp.cwd('/pathogen/Antimicrobial_resistance/AMRFinder/data')               # change into "debian" directory
    ls = []
    ftp.retrlines('MLSD', ls.append)
    latest = get_latest(ls)
    #print("Latest = {}".format(latest))
    mkdir_p(latest)
    os.chdir(latest)
    ftp.cwd(latest)
    files = []
    ftp.retrlines('NLST', files.append)
    for f in files:
        print("  Fetching {}...".format(f), end='') 
        ftp.retrbinary('RETR {}'.format(f), open(f, 'wb').write)
        print("success!")
        
    os.chdir("..")
    if os.path.islink("latest"):
        os.unlink("latest")
    os.symlink(latest, "latest")
        
    os.chdir(prevdir)
    
class cwlgen:
    def __init__(self, args):
        self.args = args
        self.parse_deflines = True
        self.do_protein = True if self.args.protein else False
            
    def prot_params(self):
        p = {
        'query': {                                                      
            'class': 'File',
            'location': os.path.realpath(self.args.fasta)
            },
        'fasta': {
            'class': 'File',
            'location': os.path.realpath(self.args.fastadb)
            },
        'hmmdb': {
            'class': 'File',
            'location': os.path.realpath(self.args.hmmdb)
            },
        'fam': {
            'class': 'File',
            'location': os.path.realpath(self.args.fam)
            },
        'parse_deflines': self.parse_deflines
        }
        if self.args.gff:
            p['gff'] = {
                'class': 'File',
                'location': os.path.realpath(self.args.gff)
            }
        if self.args.num_threads:
            p['num_threads'] = self.args.num_threads
            p['cpu'] = self.args.num_threads
        return p
    
    def dna_params(self):
        p = {
            'query': {                                                      
                'class': 'File',
                'location': os.path.realpath(self.args.fasta)
                },
            'fasta': {
                'class': 'File',
                'location': os.path.realpath(self.args.fastadb)
                },
            'fam': {
                'class': 'File',
                'location': os.path.realpath(self.args.fam)
                },
            'parse_deflines': self.parse_deflines,
            'ident_min': self.args.ident_min,
            'complete_cover_min': self.args.coverage_min,
            'query_gencode': self.args.translation_table
            }
        if self.args.num_threads:
            p['num_threads'] = self.args.num_threads
        return p
    
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
        if self.args.parallel:
            cwlcmd.extend(['--parallel'])
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
        except subprocess.CalledProcessError as eCPE:
            print(eCPE.cmd)
            print("Return code:", eCPE.returncode)
            print(eCPE.output)
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

class FastaAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if option_string == '-p' or option_string == '--protein':
            setattr(namespace, 'protein', True)
        if option_string == '-n' or option_string == '--nucleotide':            
            setattr(namespace, 'protein', False)
        #print('%r %r %r' % (namespace, values, option_string))
        setattr(namespace, self.dest, values)
        
def run(updater_parser):
    parser = argparse.ArgumentParser(
        parents=[updater_parser],
        description='Run (and optionally update) the amr_finder pipeline.')
    parser.add_argument('-U', '--update-data', action='store_true', help='Update auxillary data from the ftp site (default: %(default)s)')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-p', '--protein', dest='fasta', action=FastaAction,
                       help='Amino-acid sequences to search using BLASTP and HMMER')
    group.add_argument('-n', '--nucleotide', dest='fasta', action=FastaAction,
                       help='Genomic sequence to search using BLASTX')
    
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

    parser.add_argument('--fastadb', type=str,
                        help='FASTA database to be searched (default: %(default)s).')
    parser.add_argument('--hmmdb', type=str,
                        help='HMMR database to be searched (default: %(default)s).')
    parser.add_argument('--fam', type=str,
                        help='FAM file for HMMR (default: %(default)s).')

    pre = os.path.dirname(os.path.realpath(__file__)) + "/data/latest/"
    
    parser.add_argument('-s', '--show_output', action='store_true',
                        help='Show the stdout and stderr output from the pipeline execution (verbose mode, useful for debugging).')

    parser.add_argument('-P', '--parallel', action='store_true',
                        help='[experimental] Run jobs in parallel. Does not currently keep track of ResourceRequirements like the number of cores or memory and can overload this system.')
    parser.add_argument('-N', '--num_threads', type=int,
                        help='Number of threads to use for blast/hmmr (default: %(default)s).')
    max_cpus = min(8, available_cpu_count())
    parser.set_defaults(ident_min=0.9,
                        coverage_min=0.9,
                        translation_table=11,
                        fastadb = pre + 'AMRProt',
                        hmmdb = pre + 'AMR.LIB',
                        fam = pre + 'fam.tab',
                        num_threads=max_cpus)
    
    args = parser.parse_args()

    if args.version:
        print_versions()
        sys.exit()

    has_data = check_data()
    if not has_data:
        print("Required supplementary data not present, downloading via ftp.")
        update_data()
        
    g = cwlgen(args)
    g.params()
    g.run()
    
        
    
