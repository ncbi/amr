// amrfinder_update.cpp

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
*   Updating of AMRFinder data
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
	
	
	
// ThisApplication

struct ThisApplication : ShellApplication
{
  ThisApplication ()
    : ShellApplication ("Identify AMR genes in proteins and/or contigs and print a report", true, true, true)
    {
    	addKey ("database", "Directory for AMRFinder database", "", 'd', "DATABASE_DIR");
      addFlag ("quiet", "Suppress messages to STDERR", 'q');
	  #ifdef SVN_REV
	    version = SVN_REV;
	  #endif
    }


  
  void download (const string &dbDir,
                 const string &fName) const
  {
    ASSERT (isDirName (dbDir));
    ASSERT (! fName. empty ());    
    exec ("wget 'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinder/data/latest/" + fName + "' -O " + dbDir + fName + " -o /dev/null");
  }



  void body () const final
  {
          string dbDir = getArg ("database");
    const bool   quiet = getFlag ("quiet");
    
    
    Stderr stderr (quiet);
    stderr << "Running "<< getCommandLine () << '\n';
    const Verbose vrb (qc_on);
    
    
    // dbDir
    if (! isRight (dbDir, "/"))
      dbDir += "/";    
    if (! directoryExists (dbDir))
      throw runtime_error ("Directory " + strQuote (dbDir) + " does not exist");
    
    findProg ("makeblastdb");
    
    
    string latest_version;
    exec ("wget 'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinder/data/' -O " + tmp + ".html  -o /dev/null");
    exec ("grep '<a href=\"latest/\">latest/</a>' -B 1 " + tmp + ".html | head -1 | sed 's/^[^>]*>//1' | sed 's|/<.*$||1' > " + tmp + ".txt");
    {
      LineInput li (tmp + ".txt");
      latest_version = li. getString ();
      if (latest_version. empty ())
        throw runtime_error ("Cannot get the latest version");
    }
    
    stderr << "Dowloading AMRFinder database version " << latest_version << "\n";
    download (dbDir, "AMR.LIB");
    download (dbDir, "AMRProt");
    download (dbDir, "AMR_CDS");
    download (dbDir, "fam.tab");
    download (dbDir, "changes.txt");

	  exec (fullProg ("makeblastdb") + " -in " + dbDir + "AMRProt  -dbtype prot  -logfile /dev/null");  
  }
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



