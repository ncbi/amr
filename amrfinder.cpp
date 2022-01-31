// amrfinder.cpp

/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE                          
*               National Center for Biotechnology Information
*                                                                          
*  This software/database is a "United States Government Work" under the   
*  terms of the United States Copyright Act.  It was written as part of    
*  the author's official duties as a United States Government employee and 
*  thus cannot be copyrighted.  This software/database is freely available 
*  to the public for use. The National Library of Medicine and the U.S.    
*  Government have not placed any restriction on its use or reproduction.  
*                                                                          
*  Although all reasonable efforts have been taken to ensure the accuracy  
*  and reliability of the software and data, the NLM and the U.S.          
*  Government do not and cannot warrant the performance or results that    
*  may be obtained by using this software or data. The NLM and the U.S.    
*  Government disclaim all warranties, express or implied, including       
*  warranties of performance, merchantability or fitness for any particular
*  purpose.                                                                
*                                                                          
*  Please cite NCBI in any work or product based on this material.   
*
* ===========================================================================
*
* Author: Vyacheslav Brover
*
* File Description:
*   AMRFinder
*
* Dependencies: NCBI BLAST, HMMer
*               cat, cp, cut, head, ln, mv, sort, tail
*
* Release changes:
*   3.10.21 01/31/2022 PD-4071  https --> ftp
*   3.10.20 01/13/2022 PD-4069  "Vulnerability Disclosure" at the AMRFinder download web page
*   3.10.19 12/08/2021          $PATH is printed is case of error
*   3.10.18 11/01/2021 PD-4013  multi-domain protein hits
*   3.10.17 09/27/2021          program "dna_mutation" is found if the DNA file is empty
*   3.10.16 09/22/2021 PD-3958  point mutations oiverride HMM hits
*   3.10.15 08/23/2021          -mt_mode 1 is used only for multi-FASTA files: SB-3162
*   3.10.14 08/19/2021 PD-3918  BLAST output to stderr is not checked
*   3.10.13 08/19/2021 PD-3918  BLAST -mt_mode 1 is turned off for Mac, see SB-3163
*   3.10.12 08/19/2021 PD-3918  BLAST output to stderr is reported as error, except for BLASTN due to SB-3162
*   3.10.11 08/18/2021 PD-3826  dashes in a protein FASTA file are removed with a warning (crashes with some versions of HMMer)
*   3.10.10 08/16/2021 PD-3910  alien organism's proteins are removed from processing in amr_report.cpp (point mutations, susceptible)
*   3.10.9  08/13/2021 PD-3888  temporary files are named "amrfinder.XXXXXX"
*   3.10.8  07/06/2021 PD-3865  creating a directory does not break if upper-level directories are not readable
*                      PD-3867  "BLAST -mt_mode 1" is used
*                               makeblastdb and tblastn are used for large sequences 
*                               hmmsearch and tbasltn: query file is split into threads_max files
*                              "blastp -seg no" is used
*   3.10.7  05/18/2021 PD-3820  message for missing AMRProt blast database; https://github.com/ncbi/amr/issues/53
*   3.10.6  05/07/2021 PD-3796  for POINTN reported "target length" column = targetEnd - targetStart 
*   3.10.5  04/12/2021 PD-3772  --report_equidistant
*   3.10.4  03/24/2021 PD-3761  amrfinder --help will not break if /tmp is full
*   3.10.3  03/15/2021 PD-3749  --nucleotide_flank5_output, --nucleotide_flank5_size
*   3.10.2  03/03/2021 PD-3729  Neighboring point mutations are reported
*   3.10.1  02/17/2021 PD-3679  AMRProt-susceptible.tab
*   3.9.10  02/16/2021 PD-3694  message about missing "latest/" symbolic link; amrfinder_update.cpp: createLatestLink()
*   3.9.9   01/27/2021 PD-3674  crash for a custom database
*   3.9.8   01/04/2021 PD-3613  --dir is removed
*   3.9.7   12/03/2020 PD-3292  dependence on "mkdir" is removed
*   3.9.6   11/20/2020 PD-3613  --dir
*                               prepare_fasta_extract()
*   3.9.5   11/18/2020 PD-3292  dependence on awk is removed
*                               --help prints instruction on $TMPDIR
*   3.9.4   11/16/2020 PD-3609  ($TMPDIR or "/tmp") + "/XXXXXX"
*   3.9.3   11/05/2020 PD-3577  Merge lines for bifunctional proteins
*   3.9.2   11/04/2020 PD-3590  AMRProt has new fields #9 and #10: "subclass" and "class"
*   3.9.1   10/27/2020 PD-3583  AMRProt has a new field #8 "reportable"
*           09/30/2020 PD-2407  option --type is removed
*   3.8.28  09/29/2020 PD-3292  dependence on "uniq" is removed
*   3.8.27  09/28/2020 PD-2381  non-standard start codons are not changed in fusion proteins
*   3.8.26  09/25/2020 PD-2381  proteins with non-standard start codons that are extended in the N-terminal direction are EXACTP
*   3.8.25  09/25/2020 PD-3547  identification of frameshifts is disabled
*                               POINTX method with more SNPs is preferred over POINTP method
*   3.8.24  09/21/2020 PD-3536  --pointmut_all reports all SNPs in a reference gene repetition
*   3.8.23  09/16/2020 PD-3536  simplifying point mutations preference
*   3.8.22  09/15/2020 PD-3470  frameshift detection bug; preference of point mutation reference proteins 
*   3.8.21  09/14/2020 PD-3536  point mutations merging bug
*                      PD-3469  --force_update implies --update; -U
*   3.8.20  09/14/2020 PD-3531  "--parm -print_fam" bug
*   3.8.19  09/04/2020 PD-3292  removed the dependence on "grep"
*   3.8.18  09/03/2020 PD-3292  removed the dependence on "which"
*   3.8.17  09/02/2020 PD-3528  ordering of rows in the report is broken with parameter --name
*   3.8.16  09/01/2020 PD-2322  a complete nucleotide hit is not preferred to a partial protein hit; stopCodon field is borrowed from BLASTX to BPASTP
*   3.8.15  08/28/2020 PD-3475  Return BLAST alignment parameters for HMM-only hits where available
*   3.8.14  08/27/2020 PD-3470  method FRAME_SHIFT, amr_report is faster
*   3.8.13  08/25/2020 PD-2322  a complete nucleotide hit is preferred to a partial protein hit
*   3.8.12  08/24/2020 PD-2394  fusion genes are reported to include both gene symbols on each line
*   3.8.11  08/21/2020 PD-2407  --type
*   3.8.10  08/20/2020 PD-3469  --force_update
*   3.8.9   08/13/2020          BLAST -show_gis parameter is removed, more mutations are reported for --mutation_all
*   3.8.8   08/04/2020          bug in fasta_extract.cpp, more output in --nucleotide_output
*   3.8.7   08/03/2020 PD-3504  --protein_output, --nucleotide_output options by fasta_extract.cpp
*   3.8.6   07/29/2020 PD-3468  --name option
*           07/13/2020 PD-3484  -l for old database versions
*   3.8.5   07/10/2020 PD-3482  --ident_min instruction
*   3.8.4   05/13/2020 PD-3447  Custom point mutation does not match the reference sequence
*                               Text "*** ERROR ***" is not repeated twice
*   3.8.3   05/01/2020          WILDTYPE mutations were reported as 0-based
*   3.8.2   05/01/2020 PD-3419  taxgroup is removed from the DNA files, dna_mutation parameter: organism
*                      PD-3437  --mutation_all requires --organism
*                               all warnings are printed to stderr
*                               warnings are printed in bright yellow color; ERROR is printed in bright red color
*                      PD-3363  WILDTYPE mutations map on the reference gene with offset
*                               NOVEL is changed to UNKNOWN
*   3.8.1   04/30/2020 PD-3419  dna_mutation: reporting gene symbol for novel mutations; taxgroup and genesymbol are added to the DNA files
*   3.7.6   04/29/2020 PD-3419  dna_mutation: reporting gene symbol for novel mutations
*   3.7.5   04/22/2020 PD-3427  -h prints the help message
*   3.7.4   04/14/2020 PD-3391  Mac Conda installation
*   3.7.3   04/09/2020 PD-3416  Redundant QC check in alignment.cpp
*   3.7.2   04/08/2020 PD-3363  "WILDTYPE" was not reported
*   3.7.1   04/02/2020 PD-3154  GIs may be 0, accessions are main identifiers; file "AMRProt-suppress" is added accessions; DATA_VER_MIN is "2020-04-02.1"
*   3.6.19  03/24/2020          Check of ">lcl|" is done only for the first sequence in FASTA
*           03/24/2020 PD-3347  -lcl parameter in gff_check and amr_report
*   3.6.18  03/17/2020 PD-3396  amr_report.cpp prints a better error message on missing sublcass in data
*   3.6.17  03/12/2020          Software version is printed after software directory
*   3.6.16  03/06/2020 PD-3363  --mutation_all: UNKNOWN are not reported
*                      PD-2328  Last 2 columns of report are real HMM hits
*   3.6.15  02/24/2020          "database" is printed to stderr in one line in a canonical form (without links)
*   3.6.14  02/19/2020 PD-3363   --mutation_all: 4 types of mutations, adding DNA mutations
*   3.6.13  02/13/2020 PD-3359,issue#23   ln -s <db>: uses path2canonical()
*   3.6.12  02/13/2020 PD-3359,issue#23   AMRFinder database directory may contain spaces
*   3.6.11  02/13/2020 PD-3359,issue#23   AMRFinder code directory may contain spaces
*   3.6.10  02/06/2020 PD-3357,issue#21  --mutation_all bug
*           01/24/2020 PD-3345   Improved error message for "GFF file mismatch"
*   3.6.9   01/13/2020           "Database directory" is printed to stederr
*           01/10/2020 PD-3329   ln -s .../amrfinder abc: abc calls the right executables
*           01/20/2020           'rm" dependence is removed
*   3.6.8   01/10/2020           'gnl|' processing is simplified
*           01/09/2020 PD-3327   allow empty input files
*   3.6.7   01/09/2020           do not remove 'lcl|' from DNA FASTA
*   3.6.6   01/09/2020 PD-3326   'gnl|' is added only for gnl|PROJECT|ACC accessions if --pgap
*           01/09/2020 PD-3324   pipefail requires bash
*           01/08/2020 GP-28123  'gnl|' is added to report if --pgap
*   3.6.5                        sorting of reported rows: gene symbol is used as the last sorting column if contig is available
*   3.6.4   01/03/2020 PD-3230   sorting of reported rows: protein accession is ignored if contig is available
*   3.6.3   01/03/2020 PD-3230   sorting of reported rows
*           12/28/2019           QC in dna_mutation
*   3.6.2   12/27/2019 PD-3230   Redundant reported lines are removed for mutated reference proteins
*                                Reports are sorted by sort
*   3.6.1   12/27/2019 PD-3230   Mutated proteins are added to AMRProt
*   3.5.10  12/20/2019           --log
*   3.5.9   12/19/2019 PD-3294   blastx parameters: space added
*   3.5.8   12/18/2019 issues/19 changed message if db path is bad
*   3.5.7   12/18/2019 PD-3289   improved message for gff_check failure
*   3.5.6   12/18/2019 PD-3269   --gpipe is removed, --pgapx is replaced by --pgap
*   3.5.5   12/17/2019 PD-3287   short proteins at an end of a contig are reported
*   3.5.4   12/17/2019 PD-3287   truncated short proteins are not reported
*   3.5.3   12/16/2019 PD-3279   GPipe-GenColl assemblies, --gpipe_org
*                      GP-28025
*   3.5.2   12/13/2019 PD-3269   New flag --pgapx
*   3.5.1   12/12/2019 PD-3277   Files AMRProt-mutation.tab, AMRProt-suppress, AMR_DNA-<TAXGROUP>.tab and taxgroup.tab have headers
*   3.4.3   12/11/2019 PD-2171   --mutation_all bug
*                                --debug does not imply "-verbose 1"
*   3.4.2   12/10/2019 PD-3209   alignment correction for mutations
*                                point_mut.{hpp,cpp} -> alignment.{hpp,cpp}
*                                dna_point_mut.cpp -> dna_mutation.cpp
*                                AMRProt-point_mut.tab -> AMRProt-mutation.tab
*                                protein resistance: "point_mutation" -> "mutation"
*                                amrfinder: --point_mut_all -> --mutation_all
*                      PD-3232   mutation detection redesign
*                      PD-3267   mutation in a mutated context
*   3.4.1   12/03/2019 PD-3193   AMR_DNA-*.tab: column "genesymbol" is removed
*                                product name is fixed for point mutations
*                                point_mut.cpp -> dna_point_mut.cpp
*   3.3.2   11/26/2019 PD-3193   Indel mutations: partially implemented
*                                Bug fixed: protein point mutations were reported incorrectly if there was an offset w.r.t. the reference sequence
*                                Files AMRProt-point_mut.tab and AMR_DNA-<taxgroup>.tab: columns allele, symbol are removed
*                                Files taxgroup.list and gpipe.tab are replaced by taxgroup.tab
*   3.3.1   11/22/2019 PD-3206   New files: taxgroup.list, gpipe.tab; new option --list_organisms
*   3.2.4   11/15/2019 PD-3191   dna_mutation.cpp: neighborhoodMismatch <= 0.04; good(): length >= min (refLen, 2 * flankingLen + 1)
*   3.2.3   11/14/2019 PD-3192   Fixed error made by PD-3190
*   3.2.3   11/13/2019 PD-3190   organisms for --gpipe
*   3.2.3   11/12/2019 PD-3187   Sequence name is always from AMRProt, not from fam.tab
*   3.2.2   11/06/2019 PD-2244   Added "LANG=C" before "sort"
*
*/


#ifdef _MSC_VER
  #error "UNIX is required"
#endif
   
#undef NDEBUG 
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;



#undef DIR  // PD-3613


// PAR!
// PD-3051
#define DATA_VER_MIN "2021-02-18.1"  



namespace 
{


// PAR
constexpr size_t threads_max_min = 1;  
constexpr size_t threads_def = 4;
// Cf. amr_report.cpp
constexpr double ident_min_def = 0.9;
constexpr double partial_coverage_min_def = 0.5;

    
#define HELP  \
"Identify AMR and virulence genes in proteins and/or contigs and print a report\n" \
"\n" \
"DOCUMENTATION\n" \
"    See https://github.com/ncbi/amr/wiki for full documentation\n" \
"\n" \
"UPDATES\n" \
"    Subscribe to the amrfinder-announce mailing list for database and software update notifications:\n" \
"    https://www.ncbi.nlm.nih.gov/mailman/listinfo/amrfinder-announce"




struct ThisApplication : ShellApplication
{
  ThisApplication ()
    : ShellApplication (HELP, true, true, true)
    {
    	addFlag ("update", "Update the AMRFinder database", 'u');  // PD-2379
    	addFlag ("force_update", "Force updating the AMRFinder database", 'U');  // PD-3469
    #ifdef DIR
      addKey ("dir", "Common directory of the --protein, --nucleotide and --gff files", "", '\0', "DIRECTORY");
    #endif
    	addKey ("protein", "Input protein FASTA file", "", 'p', "PROT_FASTA");
    	addKey ("nucleotide", "Input nucleotide FASTA file", "", 'n', "NUC_FASTA");
    	addKey ("gff", "GFF file for protein locations. Protein id should be in the attribute 'Name=<id>' (9th field) of the rows with type 'CDS' or 'gene' (3rd field).", "", 'g', "GFF_FILE");
      addFlag ("pgap", "Input files PROT_FASTA, NUC_FASTA and GFF_FILE are created by the NCBI PGAP");
    	addKey ("database", "Alternative directory with AMRFinder database. Default: $AMRFINDER_DB", "", 'd', "DATABASE_DIR");
    	addKey ("ident_min", "Minimum proportion of identical amino acids in alignment for hit (0..1). -1 means use a curated threshold if it exists and " + toString (ident_min_def) + " otherwise", "-1", 'i', "MIN_IDENT");
    	  // "PD-3482
    	addKey ("coverage_min", "Minimum coverage of the reference protein (0..1)", toString (partial_coverage_min_def), 'c', "MIN_COV");
      addKey ("organism", "Taxonomy group. To see all possible taxonomy groups use the --list_organisms flag", "", 'O', "ORGANISM");
      addFlag ("list_organisms", "Print the list of all possible taxonomy groups for mutations identification and exit", 'l');
    	addKey ("translation_table", "NCBI genetic code for translated BLAST", "11", 't', "TRANSLATION_TABLE");
    	addFlag ("plus", "Add the plus genes to the report");  // PD-2789
      addFlag ("report_common", "Report proteins common to a taxonomy group");  // PD-2756
    	addKey ("mutation_all", "File to report all mutations", "", '\0', "MUT_ALL_FILE");
    //addKey ("type", "Limit search to specific element types: " + all_types. toString (",") + ". A comma delimited list, case-insensitive", "", '\0', "TYPE");
    	  // "Element type" is a column name in the report
    	addKey ("blast_bin", "Directory for BLAST. Deafult: $BLAST_BIN", "", '\0', "BLAST_DIR");
    //addKey ("hmmer_bin" ??
      addFlag ("report_all_equal", "Report all equally-scoring BLAST and HMM matches");  // PD-3772
      addKey ("name", "Text to be added as the first column \"name\" to all rows of the report, for example it can be an assembly name", "", '\0', "NAME");
      addKey ("output", "Write output to OUTPUT_FILE instead of STDOUT", "", 'o', "OUTPUT_FILE");
      addKey ("protein_output", "Output protein FASTA file of reported proteins", "", '\0', "PROT_FASTA_OUT");
      addKey ("nucleotide_output", "Output nucleotide FASTA file of reported nucleotide sequences", "", '\0', "NUC_FASTA_OUT");
      addKey ("nucleotide_flank5_output", "Output nucleotide FASTA file of reported nucleotide sequences with 5' flanking sequences", "", '\0', "NUC_FLANK5_FASTA_OUT");
      addKey ("nucleotide_flank5_size", "5' flanking sequence size for NUC_FLANK5_FASTA_OUT", "0", '\0', "NUC_FLANK5_SIZE");
      addFlag ("quiet", "Suppress messages to STDERR", 'q');
      addFlag ("gpipe_org", "NCBI internal GPipe organism names");
    	addKey ("parm", "amr_report parameters for testing: -nosame -noblast -skip_hmm_check -bed", "", '\0', "PARM");
	    version = SVN_REV;  
    }



  void initEnvironment () final
  {
    ShellApplication::initEnvironment ();
    var_cast (name2arg ["threads"] -> asKey ()) -> defaultValue = to_string (threads_def);  
  }
  
  
  
  bool fastaCheck (const string &fName, 
                   bool prot, 
                   const string &qcS, 
                   const string &logFName, 
                   size_t &nSeq, 
                   size_t &len_max,
                   size_t &len_total) const
  // Input: fName: quoted
  // Return: false <=> hyphen in protein FASTA
  {
    try
    {
	    exec (fullProg ("fasta_check") + fName + "  " + (prot ? "-aa" : "-len "+ tmp + ".len") + (prot ? "" : "  -hyphen") + qcS + "  -log " + logFName + " > " + tmp + ".nseq", logFName); 
	  }
	  catch (...)
	  {
  	  if (prot)
  	  {
  	    LineInput f (logFName);
  	    while (f. nextLine ())
    	    if (contains (f. line, "hyphen in the sequence"))  // Cf. fasta_check.cpp
    	      return false;
    	}
    	throw;
	  }
  	const StringVector vec (tmp + ".nseq", (size_t) 10, true); 
  	if (vec. size () != 3)
      throw runtime_error (string (prot ? "Protein" : "DNA") + " fasta_check failed: " + vec. toString ("\n"));
    nSeq      = str2<size_t> (vec [0]);
    len_max   = str2<size_t> (vec [1]);
    len_total = str2<size_t> (vec [2]);
    QC_ASSERT (nSeq);
    QC_ASSERT (len_max);
    QC_ASSERT (len_total);
    return true;
  }

  
  
  string get_num_threads_param (const string &blast,
                                size_t threads_max_max) const
  {
    const size_t t = min (threads_max, threads_max_max);
    if (t <= 1)  // One thread is main
      return string ();
    
		bool num_threadsP = false;
		bool mt_modeP = false;
		{
      exec (fullProg (blast) + " -help > " + tmp + ".blast_help");
      LineInput f (tmp + ".blast_help");
      while (f. nextLine ())
      {
        trim (f. line);
        if (contains (f. line, "-num_threads "))
          num_threadsP = true;
        if (contains (f. line, "-mt_mode "))
          mt_modeP = true;
      }
    }
    
    if (! num_threadsP)
      return string ();
    
	  string s ("  -num_threads " + to_string (t));
	#ifndef __APPLE__
	  if (mt_modeP)  
	    s += "  -mt_mode 1";
	#endif
	    
	  return s;
  }



#if 0
  void checkBlastErr (const string &errFName) const
  {
  	const StringVector blastErr (errFName, (size_t) 10, true);  // PAR
  	if (! blastErr. empty ())
		  throw runtime_error (blastErr. toString ("\n"));
  }
#endif
  


  StringVector db2organisms () const
  {
		checkFile (tmp + ".db/taxgroup.tab");
		checkFile (tmp + ".db/AMRProt-mutation.tab");
		checkFile (tmp + ".db/AMRProt-susceptible.tab");
    exec ("tail -n +2 " + tmp + ".db/taxgroup.tab" + "            | cut -f 1 > " + tmp + ".tax_org");
    exec ("tail -n +2 " + tmp + ".db/AMRProt-mutation.tab" + "    | cut -f 1 > " + tmp + ".prot_org");
    exec ("tail -n +2 " + tmp + ".db/AMRProt-susceptible.tab" + " | cut -f 1 > " + tmp + ".susc_org");
    exec ("cat " + tmp + ".tax_org " + tmp + ".prot_org " + tmp + ".susc_org | sort -u > " + tmp + ".org");
    return StringVector (tmp + ".org", (size_t) 100, true);  // PAR
  }
  
  
  
  string col2num (const string &colName) const
  // Return: number
  // Input: tmp + ".amr": must have the header line
  {
    LineInput f (tmp + ".amr");
    EXEC_ASSERT (f. nextLine ());
    const List<string> columns (str2list (f. line, '\t'));
    size_t n = 1;
    for (const string& column : columns)
      if (column == colName)
        return to_string (n);
      else
        n++;
    throw runtime_error ("Column " + strQuote (colName) + " not found in " + tmp + ".amr");    
  }
  
  
  
  struct SortField : Named
  {
    bool numeric {false};
    
    explicit SortField (const string &name_arg,
                        bool numeric_arg = false)
      : Named (name_arg)
      , numeric (numeric_arg)
      {
        ASSERT (str2<int> (name) > 0);
      }
    void saveText (ostream &os) const
      { os << "-k" << name << ',' << name;
        if (numeric)
          os << 'n';
      }
  };



  void prepare_fasta_extract (StringVector &&columns,
                              const string &tmpSuf,
                              bool saveHeader) const
  // Input: tmp + ".amr"
  {
    TextTable t (tmp + ".amr");
    t. qc ();
    t. filterColumns (move (columns));
    t. rows. filterValue ([] (const StringVector& row) { return row [0] == "NA"; });
    t. rows. sort ();
    t. rows. uniq ();
    t. saveHeader = saveHeader;
    t. qc ();
    t. saveFile (tmp + "." + tmpSuf);
  }



  void shellBody () const final
  {
    const bool   force_update    =             getFlag ("force_update");
    const bool   update          =             getFlag ("update") || force_update;
    const string dir             =    
  #ifdef DIR
                                   appendS (getArg ("dir"), "/"); 
  #else
                                   "";
  #endif
    const string prot            = shellQuote (prependS (getArg ("protein"),    dir));
    const string dna             = shellQuote (prependS (getArg ("nucleotide"), dir));
    const string gff             = shellQuote (prependS (getArg ("gff"),        dir));
          string db              =             getArg ("database");
    const bool   pgap            =             getFlag ("pgap");
    const double ident           =             arg2double ("ident_min");
    const double cov             =             arg2double ("coverage_min");
    const string organism        = shellQuote (getArg ("organism"));   
    const bool   list_organisms  =             getFlag ("list_organisms");
    const uint   gencode         =             arg2uint ("translation_table"); 
    const bool   add_plus        =             getFlag ("plus");
    const bool   report_common   =             getFlag ("report_common");
    const string mutation_all    = shellQuote (getArg ("mutation_all"));  
  //const string type            =             getArg ("type");
          string blast_bin       =             getArg ("blast_bin");
    const bool   equidistant     =             getFlag ("report_all_equal");
    const string input_name      = shellQuote (getArg ("name"));
    const string parm            =             getArg ("parm");  
    const string output          = shellQuote (getArg ("output"));
    const string prot_out        = shellQuote (getArg ("protein_output"));
    const string dna_out         = shellQuote (getArg ("nucleotide_output"));
    const string dnaFlank5_out   = shellQuote (getArg ("nucleotide_flank5_output"));
    const uint   dnaFlank5_size  =             arg2uint ("nucleotide_flank5_size");
    const bool   quiet           =             getFlag ("quiet");
    const bool   gpipe_org       =             getFlag ("gpipe_org");
    
    
		const string logFName (tmp + ".log");  // Command-local log file


    Stderr stderr (quiet);
    stderr << "Running: "<< getCommandLine () << '\n';
    stderr << "Software directory: " << shellQuote (execDir) << "\n";
	  stderr << "Software version: " << version << '\n'; 
    
    if (contains (input_name, '\t'))
      throw runtime_error ("NAME cannot contain a tab character");

    if (threads_max < threads_max_min)
      throw runtime_error ("Number of threads cannot be less than " + to_string (threads_max_min));
    
		if (ident != -1.0 && (ident < 0.0 || ident > 1.0))
		  throw runtime_error ("ident_min must be between 0 and 1");
		
		if (cov < 0.0 || cov > 1.0)
		  throw runtime_error ("coverage_min must be between 0 and 1");
		  
	  if (report_common && emptyArg (organism))
		  throw runtime_error ("--report_common requires --organism");
	  if (report_common && ! add_plus)
		  throw runtime_error ("--report_common requires --plus");
		  
		  		  
		// PD-3437
	  if (! emptyArg (mutation_all) && emptyArg (organism))
	  {
	    const Warning warning (stderr);
		  stderr << "--mutation_all option used without -O/--organism option. No point mutations will be screened";
		}

  #if 0		
    StringVector typeVec;
		if (! type. empty ())
		{
		  const List<string> typeList (str2list (type, ','));
		  for (string s: typeList)
		  {
		    trim (s);
		    if (s. empty ())
		      continue;
		    strUpper (s);
		    if (! all_types. contains (s))
		      throw runtime_error ("Unknown element type " + strQuote (s));
		    typeVec << s;
		  }
		}
  #endif

		if (! emptyArg (output))
		  try { OFStream f (unQuote (output)); }
		    catch (...) { throw runtime_error ("Cannot create output file " + output); }

    
    // For timing... 
    const time_t start = time (NULL);
    
    
    const size_t threads_max_max = get_threads_max_max ();
    if (threads_max > threads_max_max)
    {
      stderr << "The number of threads cannot be greater than " << threads_max_max << " on this computer" << '\n'
             << "The current number of threads is " << threads_max << ", reducing to " << threads_max_max << '\n';
      threads_max = threads_max_max;
    }


    string defaultDb;
    #ifdef CONDA_DB_DIR
    // we're in condaland
      if (const char* s = getenv("CONDA_PREFIX")) {
        defaultDb = string (s) + "/share/amrfinderplus/data/latest";
      } else if (const char* s = getenv("PREFIX")) {
        const Warning warning (stderr);
        stderr << "This was compiled for running under bioconda, but $CONDA_PREFIX was not found" << "\n";
        defaultDb = string (s) + "/share/amrfinderplus/data/latest";
        stderr << "Reverting to $PREFIX: " << defaultDb;
      } else {
        const Warning warning (stderr);
        stderr << "This was compiled for running under bioconda, but $CONDA_PREFIX was not found" << "\n";
        stderr << "Reverting to hard coded directory: " << CONDA_DB_DIR "/latest";
        defaultDb = CONDA_DB_DIR "/latest";
      }
    #else
    // not in condaland
      defaultDb = execDir + "data/latest";
    #endif
    ASSERT (isRight (defaultDb, "/latest"));
        
		// db
		if (db. empty ())
		{
    	if (const char* s = getenv ("AMRFINDER_DB"))
    		db = s;
    	else
			  db = defaultDb;
		}
		ASSERT (! db. empty ());		  


		if (update)
    {
      // PD-2447
      if (! emptyArg (prot) || ! emptyArg (dna))
        throw runtime_error ("AMRFinder -u/--update option cannot be run with -n/--nucleotide or -p/--protein options");
      if (! getArg ("database"). empty ())
        throw runtime_error ("AMRFinder update option (-u/--update) only operates on the default database directory. The -d/--database option is not permitted");
      if (getenv ("AMRFINDER_DB"))
      {
        const Warning warning (stderr);
        stderr << "AMRFINDER_DB is set, but AMRFinder auto-update only downloads to the default database directory";
        db = defaultDb;
      }
  		const Dir dbDir (db);
      if (! dbDir. items. empty () && dbDir. items. back () == "latest")
      {
        prog2dir ["amrfinder_update"] = execDir;
  		  exec (fullProg ("amrfinder_update") + " -d " + shellQuote (dbDir. getParent ()) + ifS (force_update, " --force_update") 
  		          + ifS (quiet, " -q") + ifS (qc_on, " --debug") + " > " + logFName, logFName);
      }
      else
      {
        const Warning warning (stderr);
        stderr << "Updating database directory works only for databases with the default data directory format." << "\n"
               << "         Please see https://github.com/ncbi/amr/wiki for details." << "\n"
               << "         Current database directory is: " << dbDir. get () << "\n"
               << "         New database directories will be created as subdirectories of " << dbDir. getParent ();
      }
		}


    const string downloadLatestInstr ("\nTo download the latest version to the default directory run: amrfinder -u");
    
		if (! directoryExists (db))  // PD-2447
		  throw runtime_error ("No valid AMRFinder database found.\nSymbolic link is not found: " + db + ifS (! update, downloadLatestInstr));
		stderr << "Database directory: " << shellQuote (path2canonical (db)) << "\n";		
    exec ("ln -s " + shellQuote (path2canonical (db)) + " " + tmp + ".db");


		// PD-3051
		{
  	  istringstream versionIss (version);
  		const SoftwareVersion softwareVersion (versionIss);
  		const SoftwareVersion softwareVersion_min (db + "/database_format_version.txt");
  	//stderr << "Software version: " << softwareVersion. str () << '\n'; 
  		const DataVersion dataVersion (db + "/version.txt");
  		istringstream dataVersionIss (DATA_VER_MIN); 
  		const DataVersion dataVersion_min (dataVersionIss);  
      stderr << "Database version: " << dataVersion. str () << '\n';
      if (softwareVersion < softwareVersion_min)
        throw runtime_error ("Database requires sofware version at least " + softwareVersion_min. str ());
      if (dataVersion < dataVersion_min)
        throw runtime_error ("Software requires database version at least " + dataVersion_min. str () + downloadLatestInstr);
    }


    if (list_organisms)
    {
      const StringVector organisms (db2organisms ());
      cout << endl << "Available --organism options: " + organisms. toString (", ") << endl;
      return;
    }    		  

		  
    {
      string searchMode;
      StringVector includes;
      if (emptyArg (prot))
      {
        if (emptyArg (dna))
        {
          if (update)
            return;
  	  	  throw runtime_error ("Parameter --protein or --nucleotide must be present");
    		}
        else
        {
      		if (! emptyArg (gff))
            throw runtime_error ("Parameter --gff is redundant");
          searchMode = "translated nucleotide";
        }
      }
      else
      {
        searchMode = "protein";
        if (emptyArg (dna))
        {
          searchMode += "-only";
          includes << key2shortHelp ("nucleotide") + " and " + key2shortHelp ("gff") + " options to add translated searches";
        }
        else
        {
      		if (emptyArg (gff))
            throw runtime_error ("If parameters --protein and --nucleotide are present then parameter --gff must be present");
          searchMode = "combined translated and protein";
        }
      }
      if (emptyArg (prot) && ! emptyArg (prot_out))
        throw runtime_error ("Parameter --protein must be present for --protein_output");
      if (emptyArg (dna) && ! emptyArg (dna_out))
        throw runtime_error ("Parameter --nucleotide must be present for --nucleotide_output");
      if (emptyArg (dna) && ! emptyArg (dnaFlank5_out))
        throw runtime_error ("Parameter --nucleotide must be present for --nucleotide_flank5_output");
      if (emptyArg (dnaFlank5_out) && dnaFlank5_size > 0)
        throw runtime_error ("Parameter --nucleotide_flank5_output must be present for --nucleotide_flank5_size");
      if (! emptyArg (dnaFlank5_out) && dnaFlank5_size == 0)
        throw runtime_error ("Parameter --nucleotide_flank5_size must be present with a positive value for --nucleotide_flank5_output");
      ASSERT (! searchMode. empty ());
      if (emptyArg (organism))
        includes << key2shortHelp ("organism") + " option to add mutation searches and suppress common proteins";
      else
        searchMode += " and mutation";
      stderr << "AMRFinder " << searchMode << " search\n";

      for (const string& include : includes)
        stderr << "  - include " << include << '\n';

      StringVector emptyFiles;
      if (! emptyArg (prot) && ! getFileSize (unQuote (prot)))  emptyFiles << prot;
      if (! emptyArg (dna)  && ! getFileSize (unQuote (dna)))   emptyFiles << dna;
      if (! emptyArg (gff)  && ! getFileSize (unQuote (gff)))   emptyFiles << gff;      
      for (const string& emptyFile : emptyFiles)
      {
        const Warning warning (stderr);
        stderr << "Empty file: " << emptyFile;
      }
    }
      

    // blast_bin
    if (blast_bin. empty ())
    	if (const char* s = getenv ("BLAST_BIN"))
    		blast_bin = string (s);
    if (! blast_bin. empty ())
    {
	    if (! isRight (blast_bin, "/"))
	    	blast_bin += "/";
	    prog2dir ["blastp"]      = blast_bin;
	    prog2dir ["blastx"]      = blast_bin;
	    prog2dir ["tblastn"]     = blast_bin;
	    prog2dir ["blastn"]      = blast_bin;
      prog2dir ["makeblastdb"] = blast_bin;
	  }
	  
	  // organism --> organism1
	  string organism1;
	  bool suppress_common = false;	  
	  if (! emptyArg (organism))
	  {
	  	organism1 = unQuote (organism);
 	  	replace (organism1, ' ', '_');
 	  	ASSERT (! organism1. empty ());
      if (gpipe_org)
      {
        LineInput f (db + "/taxgroup.tab");
        Istringstream iss;
        bool found = false;
        while (f. nextLine ())
        {
	  	    if (isLeft (f. line, "#"))
	  	      continue;
          iss. reset (f. line);
          string org, gpipeOrg;
          int num = -1;
          iss >> org >> gpipeOrg >> num;
          QC_ASSERT (! org. empty ());
          QC_ASSERT (num >= 0);
          QC_ASSERT (iss. eof ());
          if (organism1 == gpipeOrg)
          {
            organism1 = org;
            found = true;
            break;
          }
        }
        if (! found)
          organism1. clear ();
      }
      if (! organism1. empty ())
      {
        const StringVector organisms (db2organisms ());
        if (! organisms. contains (organism1))
          throw runtime_error ("Possible organisms: " + organisms. toString (", "));  
      }
 	  }
	  if (! organism1. empty ())
 	  	if (! report_common)
 	  	  suppress_common = true;
 	  ASSERT (! contains (organism1, ' '));


    const string qcS (qc_on ? " -qc" : "");
		const string force_cds_report (! emptyArg (dna) && ! organism1. empty () ? "-force_cds_report" : "");  // Needed for dna_mutation
		
								  
    prog2dir ["fasta_check"]   = execDir;
    prog2dir ["fasta2parts"]   = execDir;
    prog2dir ["amr_report"]    = execDir;	
		prog2dir ["dna_mutation"]  = execDir;
    prog2dir ["fasta_extract"] = execDir;
    
    
    bool lcl = false;
    if (pgap && ! emptyArg (dna))  // PD-3347
    {
      LineInput f (unQuote (dna));
      while (f. nextLine ())
        if (isLeft (f. line, ">"))
        {
          lcl = isLeft (f. line, ">lcl|");
          break;
        }
    }
    

    // Create files for amr_report    
    string amr_report_blastp;	
 		string amr_report_blastx;
 		bool hmmChunks = false;
 		bool tblastnChunks = false;
	  const string pgapS (ifS (pgap, " -pgap" + ifS (lcl, " -lcl")));
    {
  		#define BLAST_FMT    "-outfmt '6 qseqid sseqid qstart qend qlen sstart send slen qseq sseq'"
  		#define TBLASTN_FMT  "-outfmt '6 sseqid qseqid sstart send slen qstart qend qlen sseq qseq'"

  		
  		// PD-2967
  		const string blastp_par ("-comp_based_stats 0  -evalue 1e-10  -seg no");  
  		  // was: -culling_limit 20  // PD-2967
  		if (! emptyArg (prot))
  		{
  			string gff_match;
  			if (getFileSize (unQuote (prot)))
  			{
    			findProg ("blastp");  			
    			findProg ("hmmsearch");
    			
    			string prot1 (prot);  // Protein FASTA with no dashes in the sequences
          size_t nProt = 0;
          size_t protLen_max = 0;
          size_t protLen_total = 0;
          if (! fastaCheck (prot, true, qcS, logFName, nProt, protLen_max, protLen_total))
          {
            prot1 = shellQuote (tmp + ".prot");
            OFStream outF (unQuote (prot1));
            LineInput f (unQuote (prot)); 
            while (f. nextLine ())
            {
              trimTrailing (f. line);
              if (f. line. empty ())
              	continue;
            	if (f. line [0] != '>')
            	  replaceStr (f. line, "-", "");
            	outF << f. line << endl;
            }
            {
        	    const Warning warning (stderr);
        		  stderr << "Ignoring dash '-' characters in the sequences of the protein file " << prot;
        		}
        		EXEC_ASSERT (fastaCheck (prot1, true, qcS, logFName, nProt, protLen_max, protLen_total));
          }
    			
   			  // gff_check
    			if (! emptyArg (gff) && ! contains (parm, "-bed"))
    			{
    			  string locus_tag;
    			  {
      			  bool locus_tagP = false;
      			  {
        			  LineInput f (unQuote (prot1));
        			  while (f. nextLine ())
        			    if (   ! f. line. empty () 
        			        && f. line [0] == '>'
        			        && contains (f. line, "[locus_tag=")
        			       )
        			    {
        			      locus_tagP = true;
        			      break;
        			    }
        			}
      			  if (locus_tagP /*|| gpipe*/)
      			  {
      			    locus_tag = " -locus_tag " + tmp + ".match";
      			    gff_match = " -gff_match " + tmp + ".match";
      			  }
      			}
    			  prog2dir ["gff_check"] = execDir;		
    			  string dnaPar;
    			  if (! emptyArg (dna))
    			    dnaPar = " -dna " + dna;
    			  try 
    			  {
    			    exec (fullProg ("gff_check") + gff + " -prot " + prot1 + dnaPar + pgapS + locus_tag + qcS + " -log " + logFName, logFName);
    			  }
    			  catch (...)
    			  {
    			    StringVector vec (logFName, (size_t) 10, false);  // PAR
    			    if (! vec. empty ())
    			      if (vec [0]. empty ())
    			        vec. eraseAt (0);
    			    if (! vec. empty ())
    			      if (vec [0] == "*** ERROR ***")
    			        vec. eraseAt (0);
    			    throw runtime_error ("GFF file mismatch.\n" + vec. toString ("\n"));  // PD-3289, PD-3345
    			  } 
    			}
    			
    			if (! fileExists (db + "/AMRProt.phr"))
    				throw runtime_error ("The BLAST database for AMRProt was not found. Use amrfinder -u to download and prepare database for AMRFinderPlus");
    				 // "BLAST database " + shellQuote (db + "/AMRProt.phr") + " does not exist");
    			
    			stderr << "Running blastp...\n";
    			{
      			const Chronometer_OnePass cop ("blastp", cerr, false, qc_on && ! quiet);
      			// " -task blastp-fast -word_size 6  -threshold 21 "  // PD-2303
      			exec (fullProg ("blastp") + " -query " + prot1 + " -db " + tmp + ".db/AMRProt" + "  " 
      			      + blastp_par + get_num_threads_param ("blastp", min (nProt, protLen_total / 10000)) + " " BLAST_FMT " -out " + tmp + ".blastp > /dev/null 2> " + tmp + ".blastp-err", logFName);
      		//checkBlastErr (tmp + ".blastp-err");
      		}
    			  
    			stderr << "Running hmmsearch...\n";
    			{
       			const Chronometer_OnePass cop ("hmmsearch", cerr, false, qc_on && ! quiet);
      			ASSERT (threads_max >= 1);
      			if (threads_max > 1 && nProt > threads_max / 2)  // PAR
      			{
        		  createDirectory (tmp + ".hmm_chunk");
        		  exec (fullProg ("fasta2parts") + prot1 + " " + to_string (threads_max) + " " + tmp + ".hmm_chunk" + qcS + " -log " + logFName, logFName);
        		  createDirectory (tmp + ".hmmsearch_dir");
        		  createDirectory (tmp + ".dom_dir");
              Threads th (threads_max - 1, true);  
        		  FileItemGenerator fig (false, true, false, tmp + ".hmm_chunk", false);
        		  string item;
        		  while (fig. next (item))
          			th. exec (fullProg ("hmmsearch") 
      			              + "  --tblout "    + tmp + ".hmmsearch_dir/" + item + "  --noali"
      			              + "  --domtblout " + tmp + ".dom_dir/"       + item + "  --cut_tc  -Z 10000  --cpu 0  " + tmp + ".db/AMR.LIB" + " " + tmp + ".hmm_chunk/" + item + " > /dev/null 2> /dev/null"
      			             );
        		  hmmChunks = true;
        	  }
        	  else
        			exec (fullProg ("hmmsearch") + " --tblout " + tmp + ".hmmsearch  --noali  --domtblout " + tmp + ".dom  --cut_tc  -Z 10000  --cpu " + to_string (threads_max - 1) + "  " + tmp + ".db/AMR.LIB" + " " + prot1 + " > /dev/null 2> /dev/null", logFName);
        	}
  		  }
  		  else
  		  {
  		    OFStream::create (tmp + ".blastp");
  		    OFStream::create (tmp + ".hmmsearch");
  		    OFStream::create (tmp + ".dom");
  		  }  

  		  amr_report_blastp = "-blastp " + tmp + ".blastp  -hmmsearch " + tmp + ".hmmsearch  -hmmdom " + tmp + ".dom";
  			if (! emptyArg (gff))
  			  amr_report_blastp += "  -gff " + gff + gff_match;
  		}  		

  		
  		if (! emptyArg (dna))
  		{
  		  const bool blastn = ! organism1. empty () && fileExists (db + "/AMR_DNA-" + organism1);
  		  if (getFileSize (unQuote (dna)))
    		{
          size_t nDna = 0;
          size_t dnaLen_max = 0;
          size_t dnaLen_total = 0;
          EXEC_ASSERT (fastaCheck (dna, false, qcS, logFName, nDna, dnaLen_max, dnaLen_total));
          const string blastx (dnaLen_max > 100000 ? "tblastn" : "blastx");  // PAR

    			stderr << "Running " << blastx << "...\n";
    			findProg (blastx);
          {
       			const Chronometer_OnePass cop (blastx, cerr, false, qc_on && ! quiet);
        		const string tblastn_par (blastp_par + "  -word_size 3  -max_target_seqs 10000");
        		const string blastx_par (tblastn_par + "  -query_gencode " + to_string (gencode));
      			ASSERT (threads_max >= 1);
      			if (blastx == "blastx")
      			{
        			exec (fullProg ("blastx") + "  -query " + dna + " -db " + tmp + ".db/AMRProt" + "  "
            			  + blastx_par + " " BLAST_FMT " " + get_num_threads_param ("blastx", min (nDna, dnaLen_total / 10002))
            			  + " -out " + tmp + ".blastx > /dev/null 2> " + tmp + ".blastx-err", logFName);
        		//checkBlastErr (tmp + ".blastx-err");
            }
            else
            {
              ASSERT (blastx == "tblastn");
        			findProg ("makeblastdb");
        			exec ("cp " + dna + " " + tmp + ".nucl");
           	  exec (fullProg ("makeblastdb") + " -in " + tmp + ".nucl" + "  -dbtype nucl  -logfile /dev/null");  
        			if (threads_max > 1)
        			{
          		  createDirectory (tmp + ".AMRProt_chunk");
          		  exec (fullProg ("fasta2parts") + tmp + ".db/AMRProt " + to_string (threads_max) + " " + tmp + ".AMRProt_chunk" + qcS + " -log " + logFName, logFName);
          		  createDirectory (tmp + ".tblastn_dir");
          		  createDirectory (tmp + ".tblastn_dir.err");
                Threads th (threads_max - 1, true);  
          		  FileItemGenerator fig (false, true, false, tmp + ".AMRProt_chunk", false);
          		  string item;
          		  while (fig. next (item))
            			th. exec (fullProg ("tblastn") + "  -db " + tmp + ".nucl  -query " + tmp + ".AMRProt_chunk/" + item + "  "
              			        + tblastn_par + " " TBLASTN_FMT "  -out " + tmp + ".tblastn_dir/" + item + " > /dev/null 2> " + tmp + ".tblastn_dir.err/" + item);
          		  tblastnChunks = true;
          	  }
          	  else
          	  {
          			exec (fullProg ("tblastn") + "  -db " + tmp + ".nucl  -query " + tmp + ".db/AMRProt  "
              			  + tblastn_par + " " TBLASTN_FMT "  -out " + tmp + ".blastx > /dev/null 2> " + tmp + ".tblastn-err", logFName);
          		//checkBlastErr (tmp + ".tblastn-err");
              }
            }
          }

          if (blastn)
      		{
      			findProg ("blastn");
      			stderr << "Running blastn...\n";
       			const Chronometer_OnePass cop ("blastn", cerr, false, qc_on && ! quiet);
      			exec (fullProg ("blastn") + " -query " + dna + " -db " + tmp + ".db/AMR_DNA-" + organism1 + " -evalue 1e-20  -dust no  " 
      			      + get_num_threads_param ("blastn", min (nDna, dnaLen_total / 2500000)) + " " BLAST_FMT " -out " + tmp + ".blastn > " + logFName + " 2> " + tmp + ".blastn-err", logFName);
      		//checkBlastErr (tmp + ".blastn-err");  // SB-3162 ??
      		}
    		}
    		else
    		{
  		    OFStream::create (tmp + ".blastx");
  		    OFStream::create (tmp + ".len");
      		if (blastn)
    		    OFStream::create (tmp + ".blastn");
  		  }
   		  amr_report_blastx = "-blastx " + tmp + ".blastx  -dna_len " + tmp + ".len";
  		}
  	}
  	
  	
  	if (hmmChunks)
  	{
  	  exec ("cat " + tmp + ".hmmsearch_dir/* > " + tmp + ".hmmsearch");
  	  exec ("cat " + tmp + ".dom_dir/*       > " + tmp + ".dom");
  	}

  	if (tblastnChunks)
  	{
  	  exec ("cat " + tmp + ".tblastn_dir/* > " + tmp + ".blastx");
  	  exec ("cat " + tmp + ".tblastn_dir.err/* > " + tmp + ".tblastn-err");
 		//checkBlastErr (tmp + ".tblastn-err");
  	}


  	if (suppress_common)
  	{
			OFStream outF (tmp + ".suppress_prot");
			LineInput f (db + "/AMRProt-suppress");
			while (f. nextLine ())
			  if (! isLeft (f. line, "#"))
  			{
  			  string org, accver;
  			  istringstream iss (f. line);
  			  iss >> org >> accver;
  			  QC_ASSERT (! accver. empty ());
  			  if (org == organism1)
  			    outF << accver << endl;
  			}
	  }
		

    // tmp + ".amr", tmp + ".mutation_all"
		stderr << "Making report...\n";
    const string nameS (emptyArg (input_name) ? "" : " -name " + input_name);
    {
 			const Chronometer_OnePass cop ("amr_report", cerr, false, qc_on && ! quiet);
      const string mutation_allS (emptyArg (mutation_all) ? "" : ("-mutation_all " + tmp + ".mutation_all"));      
      const string coreS (add_plus ? "" : " -core");
      const string equidistantS (equidistant ? " -report_equidistant" : "");
  		exec (fullProg ("amr_report") + " -fam " + shellQuote (db + "/fam.tab") + "  " + amr_report_blastp + "  " + amr_report_blastx
      		  + "  -organism " + strQuote (organism1) 
      		  + "  -mutation "    + shellQuote (db + "/AMRProt-mutation.tab") 
      		  + "  -susceptible " + shellQuote (db + "/AMRProt-susceptible.tab") 
      		  + " " + mutation_allS + " "
      		  + force_cds_report + " -pseudo" + coreS + equidistantS
      		  + (ident == -1 ? string () : "  -ident_min "    + toString (ident)) 
      		  + "  -coverage_min " + toString (cov)
      		  + ifS (suppress_common, " -suppress_prot " + tmp + ".suppress_prot") + pgapS
      		  + nameS + qcS + " " + parm + " -log " + logFName + " > " + tmp + ".amr", logFName);
  	}
		if (   ! emptyArg (dna) 
		    && ! organism1. empty ()
		    && fileExists (db + "/AMR_DNA-" + organism1)
		   )
		{
 			const Chronometer_OnePass cop ("dna_mutation", cerr, false, qc_on && ! quiet);
      const string mutation_allS (emptyArg (mutation_all) ? "" : ("-mutation_all " + tmp + ".mutation_all.dna")); 
			exec (fullProg ("dna_mutation") + tmp + ".blastn " + shellQuote (db + "/AMR_DNA-" + organism1 + ".tab") + " " + strQuote (organism1) + " " + mutation_allS 
			      + nameS + qcS + " -log " + logFName + " > " + tmp + ".amr-snp", logFName);
			exec ("tail -n +2 " + tmp + ".amr-snp >> " + tmp + ".amr");
      if (! emptyArg (mutation_all))
  			exec ("tail -n +2 " + tmp + ".mutation_all.dna >> " + tmp + ".mutation_all");
	  }

    // Column names are from amr_report.cpp

    // Sorting AMR report
    // PD-2244, PD-3230
    string sortS;
    {
      Vector<SortField> sortFields;
      if (! (emptyArg (dna) && emptyArg (gff)))
        sortFields << SortField (col2num ("Contig id"))
                   << SortField (col2num ("Start"), true)
                   << SortField (col2num ("Stop"), true)
                   << SortField (col2num ("Strand"));
      sortFields << SortField (col2num ("Protein identifier"))
                 << SortField (col2num ("Gene symbol"));
      for (const SortField& sf : sortFields)
        sortS += " " + sf. str ();
    }
		exec ("head -1 "              + tmp + ".amr                      >  " + tmp + ".amr-out");
		exec ("LANG=C && tail -n +2 " + tmp + ".amr | sort " + sortS + " >> " + tmp + ".amr-out");
 		exec ("mv " + tmp + ".amr-out " + tmp + ".amr");

    // Sorting mutation_all
    if (! emptyArg (mutation_all))
    {
  		exec ("head -1 "              + tmp + ".mutation_all                                >  " + tmp + ".mutation_all-out");
  		exec ("LANG=C && tail -n +2 " + tmp + ".mutation_all | sort -u | sort " + sortS + " >> " + tmp + ".mutation_all-out");  
  		  // "sort -u | sort <sortS>" replaces "sort <sortS> | uniq"
   		exec ("mv " + tmp + ".mutation_all-out " + mutation_all);
    }


		if (emptyArg (output))
		  exec ("cat " + tmp + ".amr");
		else
		  exec ("cp " + tmp + ".amr " + output);
		  		  

    if (! emptyArg (prot_out))
    {
      prepare_fasta_extract (StringVector {"Protein identifier", "Gene symbol", "Sequence name"}, "prot_out", false);
      exec (fullProg ("fasta_extract") + prot + " " + tmp + ".prot_out -aa" + qcS + " -log " + logFName + " > " + prot_out, logFName);  
    }
    if (! emptyArg (dna_out))
    {
      prepare_fasta_extract (StringVector {"Contig id", "Start", "Stop", "Strand", "Gene symbol", "Sequence name"}, "dna_out", false);
      exec (fullProg ("fasta_extract") + dna + " " + tmp + ".dna_out" + qcS + " -log " + logFName + " > " + dna_out, logFName);  
    }
    if (! emptyArg (dnaFlank5_out))
    {
      prepare_fasta_extract (StringVector {"Contig id", "Start", "Stop", "Strand", "Gene symbol", "Sequence name"}, "dna_out", true);
      //                                    0            1        2       3
      TextTable t (tmp + ".dna_out");
      t. qc ();
      for (StringVector& row : t. rows)
        if (row [3] == "+")
          row [1] = to_string (max (1, stoi (row [1]) - (int) dnaFlank5_size));
        else
          row [2] = to_string (stoi (row [2]) + (int) dnaFlank5_size);
      t. saveHeader = false;
      t. qc ();
      t. saveFile (tmp + ".dnaFlank5_out");
      exec (fullProg ("fasta_extract") + dna + " " + tmp + ".dnaFlank5_out" + qcS + " -log " + logFName + " > " + dnaFlank5_out, logFName);  
    }

		
    // timing the run
    const time_t end = time (NULL);
    stderr << "AMRFinder took " << end - start << " seconds to complete\n";
  }
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



