#!/bin/csh -f

set AMRDir_def = ~brovervv/work/AMR/AMRFinder/all


if ($# != 5) then
  echo "Print an AMR report"
  echo "#1: Protein multi-FASTA file"
  echo "#2: .gff file or ''. Protein id should be in the attribute 'Name=<id>' (9th field) of the rows with type 'CDS' or 'gene' (3rd field)."
  echo "#3: amr_report parameters for testing: '-nosame -noblast -skip_hmm_check' or ''"
  echo "#4: Output file with the reported AMR proteins or ''"
  echo "#5: Directory with AMR data files; '' - $AMRDir_def"
  exit 1
endif


set DEBUG = 0


set famReport_EXEC = ~brovervv/code/prod
set QC = ""
if ($DEBUG) then
  set famReport_EXEC = ~brovervv/code/genetics
  set QC = -qc
endif


set HMMER=/opt/hmmer/3.1b2/bin
if (! -d $HMMER/)  set HMMER=/usr/local/hmmer/3.1b2/bin

set AMRDir = $AMRDir_def
if ("$5" != "")  set AMRDir = "$5"

set HMMLIB=$AMRDir/AMR.LIB
if (! -e $HMMLIB) then
  echo "$HMMLIB is not found"
  exit 1
endif


set tmp = `mktemp`  
set exitcode = 1


$famReport_EXEC/fasta_check $1 -aa -hyphen
if ($?) goto quit

set gff_match = ""
if ("$2" != "") then
  set locus_tag = ""
  grep '^>' $1 | grep '\[locus_tag=' > /dev/null
  if ($? == 0) then
    set locus_tag = "-locus_tag $tmp.match"
    set gff_match = "-gff_match $tmp.match"
  endif
  $famReport_EXEC/gff_check $2  -fasta $1  $locus_tag
	if ($?) goto quit
endif

set parse_deflines = ""
#if (??)  set parse_deflines = -parse_deflines
(/usr/bin/blastp  -task blastp-fast  -db $AMRDir/AMRProt  -query $1  -show_gis  -word_size 6  -threshold 21  -evalue 1e-20  -comp_based_stats 0  \
  -num_threads 8  $parse_deflines  \
  -outfmt '6 qseqid sseqid length nident qstart qend qlen sstart send slen qseq' \
  > $tmp.blastp) \
  >& $tmp.err
if ($?) goto quit
# PD-726, SB-1643
if ($DEBUG)  cp $tmp.blastp blastp.out   
if (0) then
	if (! -z $tmp.err) then
	  cat $tmp.err
	  goto quit
	endif
endif

$HMMER/hmmsearch  --tblout $tmp.hmmsearch  --noali  --domtblout $tmp.dom  --cut_tc  -Z 10000  --cpu 8  $HMMLIB $1 >& $tmp.out
if ($?) then
#if (1) then   
 #cat $tmp.out
  echo ""
  cat $tmp.hmmsearch
  echo ""
  cat $tmp.dom
  goto quit
endif
if ($DEBUG) then
  cp $tmp.hmmsearch hmmsearch.out   
  cp $tmp.dom       dom.out   
endif

$famReport_EXEC/amr_report  -fam $AMRDir/fam.tab  -blastp $tmp.blastp  -gff "$2"  $gff_match  -hmmsearch $tmp.hmmsearch  -hmmdom $tmp.dom  -out $tmp.sseqid  $QC  -verbose $DEBUG  $3 
if ($?) goto quit
#cp $tmp.sseqid found.sseqid 

if ("$4" != "") then
  ~brovervv/code/genetics/extractFastaProt  -in $1  -target $tmp.sseqid > $4
  if ($?) goto quit
endif


set exitcode = 0
quit:
if (! $DEBUG)  rm -f $tmp.*  
exit $exitcode
