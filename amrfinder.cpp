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
*  Please cite the author in any work or product based on this material.   
*
* ===========================================================================
*
* Author: Vyacheslav Brover
*
* File Description:
*   AMRFinder
*
*/


#ifdef _MSC_VER
  #error "UNIX is required"
#endif
   
#undef NDEBUG 
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;



namespace 
{
  

// PAR
constexpr size_t threads_max_min = 4;
// Cf. amr_report.cpp
constexpr double ident_min_def = 0.9;
constexpr double partial_coverage_min_def = 0.5;
  

#define ORGANISMS "Campylobacter|Escherichia|Salmonella"
  
		

// ThisApplication

struct ThisApplication : ShellApplication
{
  ThisApplication ()
    : ShellApplication ("Identify AMR genes in proteins and/or contigs and print a report", true, true, true)
    {
    	addKey ("protein", "Protein FASTA file to search", "", 'p', "PROT_FASTA");
    	addKey ("nucleotide", "Nucleotide FASTA file to search", "", 'n', "NUC_FASTA");
    	addKey ("database", "Alternative directory with AMRFinder database. Default: $AMRFINDER_DB", "", 'd', "DATABASE_DIR");
    	addFlag ("update", "Update the AMRFinder database", 'u');
    	addKey ("gff", "GFF file for protein locations. Protein id should be in the attribute 'Name=<id>' (9th field) of the rows with type 'CDS' or 'gene' (3rd field).", "", 'g', "GFF_FILE");
    	addKey ("ident_min", "Minimum identity for nucleotide hit (0..1). -1 means use a curated threshold if it exists and " + toString (ident_min_def) + " otherwise", "-1", 'i', "MIN_IDENT");
    	addKey ("coverage_min", "Minimum coverage of the reference protein (0..1)", toString (partial_coverage_min_def), 'c', "MIN_COV");
      addKey ("organism", "Taxonomy group for point mutation assessment\n    " ORGANISMS, "", 'O', "ORGANISM");
    	addKey ("translation_table", "NCBI genetic code for translated blast", "11", 't', "TRANSLATION_TABLE");
    	addKey ("parm", "amr_report parameters for testing: -nosame -noblast -skip_hmm_check -bed", "", '\0', "PARM");
    	addKey ("point_mut_all", "File to report all target positions of reference point mutations", "", '\0', "POINT_MUT_ALL_FILE");
    	addKey ("blast_bin", "Directory for BLAST. Deafult: $BLAST_BIN", "", '\0', "BLAST_DIR");
      addKey ("output", "Write output to OUTPUT_FILE instead of STDOUT", "", 'o', "OUTPUT_FILE");
      addFlag ("quiet", "Suppress messages to STDERR", 'q');
	  #ifdef SVN_REV
	    version = SVN_REV;
	  #endif
	  #if 0
	    setRequiredGroup ("protein",    "Input");
	    setRequiredGroup ("nucleotide", "Input");
	  #endif
    }



  void initEnvironment () final
  {
    ShellApplication::initEnvironment ();
    var_cast (name2arg ["threads"] -> asKey ()) -> defaultValue = toString (threads_max_min);  
  }



  void shellBody () const final
  {
    const string prot          = shellQuote (getArg ("protein"));
    const string dna           = shellQuote (getArg ("nucleotide"));
          string db            =             getArg ("database");
    const bool   update        =             getFlag ("update");
    const string gff           = shellQuote (getArg ("gff"));
    const double ident         =             arg2double ("ident_min");
    const double cov           =             arg2double ("coverage_min");
    const string organism      = shellQuote (getArg ("organism"));   
    const uint   gencode       =             arg2uint ("translation_table"); 
    const string parm          =             getArg ("parm");  
    const string point_mut_all =             getArg ("point_mut_all");  
          string blast_bin     =             getArg ("blast_bin");
    const string output        = shellQuote (getArg ("output"));
    const bool   quiet         =             getFlag ("quiet");
    
    
		const string logFName (tmp + ".log");


    Stderr stderr (quiet);
    stderr << "Running "<< getCommandLine () << '\n';
    const Verbose vrb (qc_on);
    
    if (threads_max < threads_max_min)
      throw runtime_error ("Number of threads cannot be less than " + toString (threads_max_min));
    
    time_t start, end;  // For timing... 
    start = time (NULL);

    string mode;
    StringVector includes;
    if (emptyArg (prot))
      if (emptyArg (dna))
  		  throw runtime_error ("Parameter --prot or --nucleotide must be present");
      else
      {
    		if (! emptyArg (gff))
          throw runtime_error ("Parameter --gff is redundant");
        mode = "translated nucleotide";
      }
    else
    {
      mode = "protein";
      if (emptyArg (dna))
      {
        mode += "-only";
        includes << key2shortHelp ("nucleotide") + " and " + key2shortHelp ("gff") + " options to add translated searches";
      }
      else
      {
    		if (emptyArg (gff))
          throw runtime_error ("If parameters --prot and --nucleotide are present then parameter --gff must be present");
        mode = "combined translated plus protein";
      }
    }
    if (emptyArg (organism))
      includes << key2shortHelp ("organism") + " option to add point-mutation searches";
    else
      mode += " and point-mutation";
      
      
    const string dbSuff ("/data/latest");  // "/../data" 
    
		// db
		if (db. empty ())
		{
    	if (const char* s = getenv ("AMRFINDER_DB"))
    		db = string (s);
    	else
			  db = execDir + dbSuff;
		}
		ASSERT (! db. empty ());
		  
		if (update)
    {
      findProg ("amrfinder_update");	
		  exec (fullProg ("amrfinder_update") + " -d " + db + "/.." + ifS (quiet, " -q") + ifS (qc_on, " --debug") + " > " + logFName, logFName);
		}

		if (! directoryExists (db))
		  throw runtime_error ("Directory with data \"" + db + "\" does not exist");


    stderr << "AMRFinder " << mode << " search with database " << db << "\n";
    for (const string& include : includes)
      stderr << "  - include " << include << '\n';


		if (ident != -1 && (ident < 0 || ident > 1))
		  throw runtime_error ("ident_min must be between 0 and 1");
		
		if (cov < 0 || cov > 1)
		  throw runtime_error ("coverage_min must be between 0 and 1");
		  

		if (! emptyArg (output))
		  try { OFStream f (unQuote (output)); }
		    catch (...) { throw runtime_error ("Cannot open output file " + output); }

    
    // blast_bin
    if (blast_bin. empty ())
    	if (const char* s = getenv ("BLAST_BIN"))
    		blast_bin = string (s);
    if (! blast_bin. empty ())
    {
	    if (! isRight (blast_bin, "/"))
	    	blast_bin += "/";
	    prog2dir ["blastp"] = blast_bin;
	    prog2dir ["blastx"] = blast_bin;
	    prog2dir ["blastn"] = blast_bin;
	  }
	  
	  string organism1;
	  if (! emptyArg (organism))
	  {
	  	organism1 = unQuote (organism);
 	  	replace (organism1, ' ', '_');
 	  }


	  if (! emptyArg (organism))
	  {
 	  	string errMsg;
			try { exec ("grep -w ^" + organism1 + " " + db + "/AMRProt-point_mut.tab &> /dev/null"); }
			  catch (const runtime_error &)
			  { 
			  	errMsg = "No protein point mutations for organism " + organism;
			  }
  		if (! emptyArg (dna))
  			if (! fileExists (db + "/AMR_DNA-" + organism1))
  	      errMsg = "No DNA point mutations for organism " + organism;
  		if (! errMsg. empty ())
  		  throw runtime_error (errMsg + "\nPossible organisms: " ORGANISMS);
	  }
        

    const string qcS (qc_on ? "-qc  -verbose 1" : "");
    const string point_mut_allS (point_mut_all. empty () ? "" : ("-point_mut_all " + point_mut_all));
		const string force_cds_report (! emptyArg (dna) && ! emptyArg (organism) ? "-force_cds_report" : "");  // Needed for point_mut
		
								  
    findProg ("fasta_check");
    findProg ("fasta2parts");
    findProg ("amr_report");	
    
    
    string blastp_par;	
 		string blastx_par;
 		bool blastxChunks = false;
    {
      Threads th (threads_max - 1, true);  

  		if ( ! emptyArg (prot))
  		{
  			findProg ("blastp");
  			findProg ("hmmsearch");

  		  exec (fullProg ("fasta_check") + prot + " -aa -hyphen  -log " + logFName, logFName);  
  			
  			string gff_match;
  			if (! emptyArg (gff) && ! contains (parm, "-bed"))
  			{
  			  string locus_tag;
  			  const int status = system (("grep '^>.*\\[locus_tag=' " + prot + " > /dev/null"). c_str ());
  			  if (status == 0)
  			  {
  			    locus_tag = "-locus_tag " + tmp + ".match";
  			    gff_match = "-gff_match " + tmp + ".match";
  			  }
  			  findProg ("gff_check");		
  			  string dnaPar;
  			  if (! emptyArg (dna))
  			    dnaPar = " -dna " + dna;
  			  exec (fullProg ("gff_check") + gff + " -prot " + prot + dnaPar + " " + locus_tag + " -log " + logFName, logFName);
  			}
  			
  			if (! fileExists (db + "/AMRProt.phr"))
  				throw runtime_error ("BLAST database " + shellQuote (db + "/AMRProt") + " does not exist");
  			
  			stderr << "Running blastp...\n";
  			// " -task blastp-fast -word_size 6  -threshold 21 "  // PD-2303
  			th << thread (exec, fullProg ("blastp") + " -query " + prot + " -db " + db + "/AMRProt  -show_gis  -evalue 1e-20  -comp_based_stats 0  "
  			  "-num_threads 6  "
  			  "-outfmt '6 qseqid sseqid length nident qstart qend qlen sstart send slen qseq sseq' "
  			  "-out " + tmp + ".blastp &> /dev/null", string ());
  			stderr << "Running hmmsearch...\n";
  			th << thread (exec, fullProg ("hmmsearch") + " --tblout " + tmp + ".hmmsearch  --noali  --domtblout " + tmp + ".dom  --cut_tc  -Z 10000  --cpu 8  " + db + "/AMR.LIB " + prot + "&> " + tmp + ".out", string ());

  		  blastp_par = "-blastp " + tmp + ".blastp  -hmmsearch " + tmp + ".hmmsearch  -hmmdom " + tmp + ".dom";
  			if (! emptyArg (gff))
  			  blastp_par += "  -gff " + gff + " " + gff_match;
  		}  		
  		
  		if (! emptyArg (dna))
  		{
  			stderr << "Running blastx...\n";
  			findProg ("blastx");
  		  exec (fullProg ("fasta_check") + dna + " -hyphen  -len "+ tmp + ".len  -log " + logFName, logFName); 
  		  const size_t threadsAvailable = th. getAvailable ();
  		  ASSERT (threadsAvailable);
  		  if (threadsAvailable >= 2)
  		  {
    		  exec ("mkdir " + tmp + ".chunk");
    		  exec (fullProg ("fasta2parts") + dna + " " + toString (threadsAvailable) + " " + tmp + ".chunk  -log " + logFName, logFName);   // PAR
    		  exec ("mkdir " + tmp + ".blastx_dir");
    		  FileItemGenerator fig (false, true, tmp + ".chunk");
    		  string item;
    		  while (fig. next (item))
      			th << thread (exec, fullProg ("blastx") + "  -query " + tmp + ".chunk/" + item + " -db " + db + "/AMRProt  "
      			  "-show_gis  -word_size 3  -evalue 1e-20  -query_gencode " + toString (gencode) + "  "
      			  "-seg no  -comp_based_stats 0  -max_target_seqs 10000  "
      			  "-outfmt '6 qseqid sseqid length nident qstart qend qlen sstart send slen qseq sseq' "
      			  "-out " + tmp + ".blastx_dir/" + item + " &> /dev/null", string ());
    		  blastxChunks = true;
  		  }
  		  else
    			th << thread (exec, fullProg ("blastx") + "  -query " + dna + " -db " + db + "/AMRProt  "
    			  "-show_gis  -word_size 3  -evalue 1e-20  -query_gencode " + toString (gencode) + "  "
    			  "-seg no  -comp_based_stats 0  -max_target_seqs 10000  -num_threads 6 "
    			  "-outfmt '6 qseqid sseqid length nident qstart qend qlen sstart send slen qseq sseq' "
    			  "-out " + tmp + ".blastx &> /dev/null", string ());
  		  blastx_par = "-blastx " + tmp + ".blastx  -dna_len " + tmp + ".len";
  		}

  		if (! emptyArg (dna) && ! emptyArg (organism))
  		{
  			ASSERT (fileExists (db + "/AMR_DNA-" + organism1));
  			findProg ("blastn");
  			findProg ("point_mut");
  			stderr << "Running blastn...\n";
  			exec (fullProg ("blastn") + " -query " +dna + " -db " + db + "/AMR_DNA-" + organism1 + " -evalue 1e-20  -dust no  "
  			  "-outfmt '6 qseqid sseqid length nident qstart qend qlen sstart send slen qseq sseq' -out " + tmp + ".blastn &> /dev/null");
  		}
  	}
  	
  	
  	if (blastxChunks)
  	  exec ("cat " + tmp + ".blastx_dir/* > " + tmp + ".blastx");
		

		exec (fullProg ("amr_report") + " -fam " + db + "/fam.tab  " + blastp_par + "  " + blastx_par
		  + "  -organism " + organism + "  -point_mut " + db + "/AMRProt-point_mut.tab " + point_mut_allS + " "
		  + force_cds_report + " -pseudo"
		  + (ident == -1 ? string () : "  -ident_min "    + toString (ident)) 
		  + "  -coverage_min " + toString (cov)
		  + " " + qcS + " " + parm + " -log " + logFName + "> " + tmp + ".amr-raw", logFName);

		if (! emptyArg (dna) && ! emptyArg (organism))
		{
			ASSERT (fileExists (db + "/AMR_DNA-" + organism1));
			exec (fullProg ("point_mut") + tmp + ".blastn " + db + "/AMR_DNA-" + organism1 + ".tab " + qcS + " -log " + logFName + " > " + tmp + ".amr-snp", logFName);
			exec ("tail -n +2 " + tmp + ".amr-snp >> " + tmp + ".amr-raw");
		}
		
		// $tmp.amr-raw --> $tmp.amr
    
    // timing the run
    end = time(NULL);
    stderr << "AMRFinder took " << end - start << " seconds to complete\n";

    string sort_cols;
    if (   ! force_cds_report. empty ()
        || ! blastx_par. empty ()
        || ! emptyArg (gff)
       )
      sort_cols = " -k2 -k3n -k4n -k5";
		exec ("head -1 " + tmp + ".amr-raw > " + tmp + ".amr");
		exec ("tail -n +2 " + tmp + ".amr-raw | sort" + sort_cols + " -k1 >> " + tmp + ".amr");


		if (emptyArg (output))
		  exec ("cat " + tmp + ".amr");
		else
		  exec ("cp " + tmp + ".amr " + output);
  }
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



