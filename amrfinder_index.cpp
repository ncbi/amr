// amrfinder_index.cpp

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
*   Indexing of AMRFinder data
*
* Dependencies: NCBI BLAST, HMMer
*
* Release changes: see amrfinder.cpp
*
*/




#ifdef _MSC_VER
  #error "UNIX is required"
#endif
   
#undef NDEBUG 
#include "common.hpp"
using namespace Common_sp;

#include "common.inc"




namespace 
{


	
// ThisApplication

struct ThisApplication : ShellApplication
{
  ThisApplication ()
    : ShellApplication ("Index the database for AMRFinder", true, false, true, true)
    {
    	addPositional ("DATABASE", "Directory with AMRFinder database");
    	addKey ("blast_bin", "Directory for BLAST", "", '\0', "BLAST_DIR");
    	addKey ("hmmer_bin", "Directory for HMMer", "", '\0', "HMMER_DIR");
	    version = SVN_REV;
    }



  void shellBody () const final
  {
    string dbDir     = getArg ("DATABASE");
    string blast_bin = getArg ("blast_bin");
    string hmmer_bin = getArg ("hmmer_bin");

    addDirSlash (dbDir);
    addDirSlash (blast_bin);
    addDirSlash (hmmer_bin);
    
        
    const Verbose vrb (qc_on);
    
    
    if (! directoryExists (dbDir))
      throw runtime_error ("Database directory " + dbDir + " does not exist");
                            
    if (! blast_bin. empty ())
      prog2dir ["makeblastdb"] = blast_bin;
    findProg ("makeblastdb");    

    if (! hmmer_bin. empty ())
      prog2dir ["hmmpress"] = hmmer_bin;
    findProg ("hmmpress");
    

    // Cf. amrfinder_update.cpp
    StringVector dnaPointMuts;
    {
      LineInput f (dbDir + "taxgroup.tsv");
      while (f. nextLine ())
      {
   	    if (isLeft (f. line, "#"))
	 	      continue;
        string taxgroup, gpipe;
        int n = -1;
        istringstream iss (f. line);
        iss >> taxgroup >> gpipe >> n;
        QC_ASSERT (n >= 0);
        if (n)
          dnaPointMuts << taxgroup;
      }
    }
    
    
    stderr. section ("Indexing");
    exec (fullProg ("hmmpress") + " -f " + shellQuote (dbDir + "AMR.LIB") + " > /dev/null 2> " + tmp + "/hmmpress.err", tmp + "/hmmpress.err");
    setSymlink (dbDir, tmp + "/db", true);
	  exec (fullProg ("makeblastdb") + " -in " + tmp + "/db/AMRProt.fa" + "  -dbtype prot  -logfile " + tmp + "/makeblastdb.AMRProt", tmp + "/makeblastdb.AMRProt");  
	  exec (fullProg ("makeblastdb") + " -in " + tmp + "/db/AMR_CDS.fa" + "  -dbtype nucl  -logfile " + tmp + "/makeblastdb.AMR_CDS", tmp + "/makeblastdb.AMR_CDS");  
    for (const string& dnaPointMut : dnaPointMuts)
  	  exec (fullProg ("makeblastdb") + " -in " + tmp + "/db/AMR_DNA-" + dnaPointMut + ".fa  -dbtype nucl  -logfile " + tmp + "/makeblastdb.AMR_DNA-" + dnaPointMut, tmp + "/makeblastdb.AMR_DNA-" + dnaPointMut);
  }
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



