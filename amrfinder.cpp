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
	
	
	
string shellQuote (string s)
{
	replaceStr (s, "\'", "\'\"\'\"\'");
	return "\'" + s + "\'";
}



bool emptyArg (const string &s)
{
	return s. empty () || s == "\'\'";
}



struct Stderr : Singleton<Stderr>
{
  bool quiet {false};
  
  explicit Stderr (bool quiet_arg)
    : quiet (quiet_arg)
    {}
    
  template <typename T>
    Stderr& operator<< (const T& t) 
      { if (quiet)
          return *this;
        cerr << t;
        return *this;
      }
};



string tmp;



string which (const string &progName)
// Return: isRight(,"/") or empty()
{
	ASSERT (! tmp. empty ());
	
	try { exec ("which " + progName + " 1> " + tmp + ".src 2> /dev/null"); }
	  catch (const runtime_error &)
	    { return string (); }
	LineInput li (tmp + ".src");
	const string s (li. getString ());
	return getDirName (s);
}

	

string execDir;
map<string,string> prog2dir;
	
	
	
void findProg (const string &progName)
// Output: prog2dir
{
	ASSERT (! progName. empty ());
	ASSERT (! contains (progName, '/'));
	ASSERT (isRight (execDir, "/"));
	
	string dir;
	if (! find (prog2dir, progName, dir))
	{
		dir = fileExists (execDir + progName)
		        ? execDir
		        : which (progName);
	  if (dir. empty ())
	  {
	  	cout << "AMRFinder binary " << shellQuote (progName) << " is not found." << endl;
		  cout << "Please make sure that " << shellQuote (progName) << " is in the same directory as " + shellQuote (Common_sp::programName) + " or is in your $PATH." << endl;
	  	exit (1);
	  }	
	  prog2dir [progName] = dir;
	}
	  
	ASSERT (isRight (dir, "/"));
}



string fullProg (const string &progName)
// Requires: After findProg(progName)
{
	string dir;
	EXEC_ASSERT (find (prog2dir, progName, dir));
	ASSERT (isRight (dir, "/"));
	return dir + progName + " ";
}




// ThisApplication

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Identify AMR genes in proteins and/or contigs and print a report", true, true)
    {
    	addKey ("protein", "Protein FASTA file to search", "", 'p', "PROT_FASTA");
    	addKey ("nucleotide", "Nucleotide FASTA file to search", "", 'n', "NUC_FASTA");
    	addKey ("database", "Alternative directory with AMRFinder database. Default: $AMRFINDER_DB ", "", 'd', "DATABASE_DIR");
    //addKey ("fasta_prot", "Create FASTA file containing protein sequence for identified proteins", "", 'f', "FASTA_PROT_OUT");
    	addKey ("gff", "GFF file for protein locations. Protein id should be in the attribute 'Name=<id>' (9th field) of the rows with type 'CDS' or 'gene' (3rd field).", "", 'g', "GFF_FILE");
    	addKey ("ident_min", "Minimum identity for nucleotide hit (0..1)", "0.9", 'i', "MIN_IDENT");
    	addKey ("coverage_min", "Minimum coverage of the reference protein to report a match as complete (0..1)", "0.9", 'c', "MIN_COV");
        addKey ("organism", "Taxonomy group for point mutation assessment\n    Campylobacter|Escherichia|Salmonella", "", 'O', "ORGANISM");
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



  void body () const final
  {
    const string prot          = shellQuote (getArg ("protein"));
    const string dna           = shellQuote (getArg ("nucleotide"));
          string db            =             getArg ("database");
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
    
    
    Stderr stderr (quiet);
    
    stderr << "Running "<< getCommandLine () << '\n';
    
    const Verbose vrb (qc_on);
    
    
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
    
    stderr << "AMRFinder " << mode << " search\n";
    for (const string& include : includes)
      stderr << "  - include " << include << '\n';


		if (ident <= 0 || ident > 1)
		  throw runtime_error ("ident_min must be between 0 and 1");
		
		if (cov <= 0 || cov > 1)
		  throw runtime_error ("coverage_min must be between 0 and 1");
		  
    
    // tmp
		tmp = tmpnam (NULL);
		if (tmp. empty ())
			throw runtime_error ("Cannot generate a temporary file");
		if (qc_on)
			cout << tmp << endl;  
		
    // execDir
		execDir = getProgramDirName ();
		if (execDir. empty ())
			execDir = which (programName);
		ASSERT (isRight (execDir, "/"));
				
		// db
		if (db. empty ())
		{
    	if (const char* s = getenv ("AMRFINDER_DB"))
    		db = string (s);
    	else
			  db = execDir + "/../data";
		}
		ASSERT (! db. empty ());
		if (! directoryExists (db))
		  throw runtime_error ("Directory with data \"" + db + "\" does not exist");

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
	  
	  if (! emptyArg (organism))
	  {
	  	string organism1 (getArg ("organism"));
 	  	replace (organism1, ' ', '_');
 	  	bool found = true;
			try { exec ("grep -w ^" + organism1 + " " + db + "/AMRProt-point_mut.tab &> /dev/null"); }
			  catch (const runtime_error &)
			  { 
			  	cout << "No protein point mutations for organism " + organism << endl;
			  	found = false;
			  }
  		if (! emptyArg (dna))
  			if (! fileExists (db + "/AMR_DNA-" + organism1))
  			{
  	      cout << "No DNA point mutations for organism " + organism << endl;
			  	found = false;
  			}
  		if (! found)
  		{
		  	cout << "Possible organisms:" << endl;
		  	exec ("cut -f 1 " + db + "/AMRProt-point_mut.tab | sort | uniq | sed 's/_/ /g'");
 				exit (1);
  		}
	  }
        

    const string qcS (qc_on ? "-qc  -verbose 1" : "");
    const string point_mut_allS (point_mut_all. empty () ? "" : ("-point_mut_all " + point_mut_all));
		const string force_cds_report (! emptyArg (dna) && ! emptyArg (organism) ? "-force_cds_report" : "");  // Needed for point_mut
		const string logFName (tmp + ".log");
			
						  
    findProg ("fasta_check");
    findProg ("amr_report");	
    
    
    string blastp_par;	
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
			{
				cout << "BLAST database " << shellQuote (db + "/AMRProt") << " does not exist" << endl;
				exit (1);
			}
			
			stderr << "Running blastp...\n";
			exec (fullProg ("blastp") + " -task blastp-fast  -query " + prot + " -db " + db + "/AMRProt  -show_gis  -word_size 6  -threshold 21  -evalue 1e-20  -comp_based_stats 0  "
			  "-num_threads 6  "
			  "-outfmt '6 qseqid sseqid length nident qstart qend qlen sstart send slen qseq sseq' "
			  "-out " + tmp + ".blastp &> /dev/null");
			stderr << "Running hmmsearch...\n";
			exec (fullProg ("hmmsearch") + " --tblout " + tmp + ".hmmsearch  --noali  --domtblout " + tmp + ".dom  --cut_tc  -Z 10000  --cpu 8  " + db + "/AMR.LIB " + prot + "&> " + tmp + ".out");

		  blastp_par = "-blastp " + tmp + ".blastp  -hmmsearch " + tmp + ".hmmsearch  -hmmdom " + tmp + ".dom";
			if (! emptyArg (gff))
			  blastp_par += "  -gff " + gff + " " + gff_match;
		}
		
		
		string blastx_par;
		if (! emptyArg (dna))
		{
			findProg ("blastx");
		  exec (fullProg ("fasta_check") + dna + " -hyphen  -len "+ tmp + ".len  -log " + logFName, logFName); 
			stderr << "Running blastx...\n";
			exec (fullProg ("blastx") + "  -query " + dna + " -db " + db + "/AMRProt  "
			  "-show_gis  -word_size 3  -evalue 1e-20  -query_gencode " + toString (gencode) + "  "
			  "-seg no  -comp_based_stats 0  -max_target_seqs 10000  "
			  "-outfmt '6 qseqid sseqid length nident qstart qend qlen sstart send slen qseq sseq' "
			  "-out " + tmp + ".blastx &> /dev/null");
		  blastx_par = "-blastx " + tmp + ".blastx  -ident_min " + getArg ("ident_min") 
		               + "  -complete_cover_min " + getArg ("coverage_min")
		               + "  -dna_len " + tmp + ".len";
		}
		

		exec (fullProg ("amr_report") + " -fam " + db + "/fam.tab  " + blastp_par + "  " + blastx_par
		  + "  -organism " + organism + "  -point_mut " + db + "/AMRProt-point_mut.tab " + point_mut_allS + " "
		  + force_cds_report + " -pseudo"
		  + " " + qcS + " " + parm + " -log " + logFName + "> " + tmp + ".amr-raw", logFName);

		
		if (! emptyArg (dna) && ! emptyArg (organism))
		{
			string organism1 = getArg ("organism");
			replace (organism1, ' ', '_');
			ASSERT (fileExists (db + "/AMR_DNA-" + organism1));
			findProg ("blastn");
			findProg ("point_mut");
			stderr << "Running blastn...\n";
			exec (fullProg ("blastn") + " -query " +dna + " -db " + db + "/AMR_DNA-" + organism1 + " -evalue 1e-20  -dust no  "
			  "-outfmt '6 qseqid sseqid length nident qstart qend qlen sstart send slen qseq sseq' -out " + tmp + ".blastn &> /dev/null");
			exec (fullProg ("point_mut") + tmp + ".blastn " + db + "/AMR_DNA-" + organism1 + ".tab " + qcS + " -log " + logFName + " > " + tmp + ".amr-snp", logFName);
			exec ("tail -n +2 " + tmp + ".amr-snp >> " + tmp + ".amr-raw");
		}
		
		
		// $tmp.amr-raw --> $tmp.amr
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
		
		
		if (! qc_on)
		  exec ("rm -f " + tmp + "*");  
  }
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



