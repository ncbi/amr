#!/bin/csh -f

set AMRDir_def = ~brovervv/work/AMR/AMRFinder/publish

if ($# != 5) then
  echo "Find AMR genes in a DNA sequence"
  echo "#1: DNA sequence file"
  echo "#2: Genetic code of DNA"
  echo "#3: Min. identity to the reference protein (0..1)"
  echo "#4: Min. coverage of the reference protein (0.5..1) to report a match as complete"
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
set AMRDir = $AMRDir_def
if ("$5" != "")  set AMRDir = "$5"


set TmpFNam = `mktemp`
set exitcode = 1


$famReport_EXEC/fasta_check $1 -hyphen
if ($?) goto quit

if (0) then
  cp $1 $TmpFNam.fa
  
  #formatdb  -i $TmpFNam.fa  -p F  -l /dev/null
  /usr/bin/makeblastdb  -in $TmpFNam.fa  -dbtype prot  -logfile /dev/null
  if ($?)  exit 1
endif

/usr/bin/blastx  -db $AMRDir/AMRProt  -query $1  \
  -show_gis  -word_size 3  -evalue 1e-20  -query_gencode $2  \
  -seg no  -comp_based_stats 0  -max_target_seqs 10000  \
  -outfmt '6 qseqid sseqid length nident qstart qend qlen sstart send slen qseq' \
  > $TmpFNam.out
if ($?) exit 1
# -parse_deflines  
# -query $TmpFNam.fa
if ($DEBUG)  cp $TmpFNam.out blastx.out  
#wc -l $TmpFNam.out

$famReport_EXEC/amr_report  -fam $AMRDir/fam.tab  -blastx $TmpFNam.out  -ident_min $3  -complete_cover_min $4  -pseudo  $QC  -verbose $DEBUG
if ($?) exit 1


set exitcode = 0
quit:
rm -f $TmpFNam.*
exit $exitcode



