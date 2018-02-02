#!/bin/bash
#DEBUG=1
DEBUG=''

# paths
BLASTP=/usr/bin/blastp
HMMER=/opt/hmmer/3.1b2/bin
AMRDir=/panfs/pan1.be-md.ncbi.nlm.nih.gov/bacterial_pathogens/backup/data/amrtests/locked_db/170130
HMMLIB=$AMRDir/AMR.LIB
#FAMREPORT=~brovervv/code/prod/famReport
FAMREPORT=/panfs/pan1.be-md.ncbi.nlm.nih.gov/bacterial_pathogens/backup/packages/amrfinder/famReport
#EXTRACT_FASTA=~brovervv/code/genetics/extractFastaProt
EXTRACT_FASTA=/panfs/pan1.be-md.ncbi.nlm.nih.gov/bacterial_pathogens/backup/packages/amrfinder/extractFastaProt

function quit {
    if [ -e "$TmpFNam" ]
    then
        [ -n "$DEBUG" ] || rm -f $TmpFNam.*
    fi
    [ -n "$DEBUG" ] && echo "TmpFNam = $TmpFNam"
    exit $exitcode
}
function usage {
      echo "
Print an AMR report from protein sequence
Usage: amrfinder-prot.sh [-options] <protein.fa>
Options:
    -g <annotation.gff> - for contig/position
    -B  - no blasts in report (only used for testing)
    -E  - no exact/allele matches (only used for testing)
    -f <out.fa> - FASTA output file of reported AMR proteins
    -p  - use -parse_deflines for blast (sometimes fixes problems caused
        - by defline formats in the input file)
    -d <amr library directory> - Use alternate AMR database directory
        - Default is $AMRDir
    "
    exitcode=1
    quit
}

BLASTS=''
EXACT=''
PARSE_DEFLINES=''
GFF=''

if [ $# -eq 0 ]
then
    usage
fi

command_line="$0 $@"

while getopts ":g:BSf:phd:" opt; do
  case $opt in
    g)
      echo "option -g $OPTARG" >&2
      GFF="$OPTARG"
      ;;
    B)
      echo "option -B" >&2
      BLASTS='-noblast'
      ;;
    E)
      echo "option -E" >&2
      EXACT='-nosame'
      ;;
    f)
      echo "option -f $OPTARG" >&2
      output_file="$OPTARG"
      ;;
    p)
      echo "option -p" >&2
      PARSE_DEFLINES='-parse_deflines'
      ;;
    d)
      echo "option -d $OPTARG" >&2
      AMRDir="$OPTARG"
      ;;
    h)
      usage
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      ;;
    :)
      echo "option -$OPTARG requires a filename." >&2
      usage
      ;;
  esac
done

shift $((OPTIND-1))

PROTEIN_FILE=$1
shift 1

if [ "X$2" != "X" ]
then
    echo "Unknown options: $@" >&2
    usage
fi


# now run amrfinder
TmpFNam=`mktemp`
exitcode=1


### BLAST
$BLASTP  -task blastp-fast  -db $AMRDir/AMRProt  -query $PROTEIN_FILE  -show_gis  -word_size 6  -threshold 21  -evalue 1e-20  -comp_based_stats 0  $PARSE_DEFLINES  -outfmt '6 qseqid sseqid length nident qstart qend qlen sstart send slen qseq' > $TmpFNam.blastp 2> $TmpFNam.err \
    || {
    echo "Error running BLAST: " >&2
    cat $TmpFNam.err >&2
    quit
} &


### HMMER
$HMMER/hmmsearch --tblout $TmpFNam.hmmsearch  --noali  --domtblout $TmpFNam.dom  --cut_tc  -Z 10000  $HMMLIB $PROTEIN_FILE 2>&1 > $TmpFNam.out \
    || {
    echo ""
    cat $TmpFNam.hmmsearch
    echo ""
    cat $TmpFNam.dom
    quit
}

wait

### famReport

$FAMREPORT -fam $AMRDir/fam.tab  -aa  -in $TmpFNam.blastp  -gff "$GFF"  -hmmsearch $TmpFNam.hmmsearch  -hmmdom $TmpFNam.dom  -out $TmpFNam.sseqid  -verbose 0 $BLASTS $EXACT  

if [ ! "$?" -eq 0 ]
then
    quit
fi

### extractFastaProt
if [ "X$output_file" != "X" ]
then
    $EXTRACT_FASTA  -in $1  -target $TmpFNam.sseqid > $output_file
    if ($?)
    then
        quit
    fi
fi

exitcode=0

quit
