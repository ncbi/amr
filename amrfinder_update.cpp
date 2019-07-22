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

#include <curl/curl.h>

#include "common.hpp"
using namespace Common_sp;

#include "amrfinder.inc"



namespace 
{



struct Curl
{
  CURL* eh {nullptr};


  Curl ()
    { eh = curl_easy_init ();
      ASSERT (eh);
    }
 ~Curl ()
   { curl_easy_cleanup (eh); }


  void download (const string &url,
                 const string &fName);
  string read (const string &url);
};

	


size_t write_stream_cb (char* ptr,
                        size_t size, 
                        size_t nMemb, 
                        void* userData)
{
  ASSERT (ptr);
  ASSERT (size == 1);
  ASSERT (userData);
  
  OFStream& f = * static_cast <OFStream*> (userData);
  FOR (size_t, i, nMemb)
    f << ptr [i];;
  
  return nMemb;
}


 	
void Curl::download (const string &url,
                     const string &fName) 
{
  ASSERT (! url. empty ());  
  ASSERT (! fName. empty ());  
  
  OFStream f (fName);
  curl_easy_setopt (eh, CURLOPT_URL, url. c_str ());
  curl_easy_setopt (eh, CURLOPT_WRITEFUNCTION, write_stream_cb);
  curl_easy_setopt (eh, CURLOPT_WRITEDATA, & f);
  EXEC_ASSERT (curl_easy_perform (eh) == 0);
}



size_t write_string_cb (char* ptr,
                        size_t size, 
                        size_t nMemb, 
                        void* userData)
{
  ASSERT (ptr);
  ASSERT (size == 1);
  ASSERT (userData);
  
  string& s = * static_cast <string*> (userData);
  FOR (size_t, i, nMemb)
    s += ptr [i];;
  
  return nMemb;
}


 	
string Curl::read (const string &url)
{
  ASSERT (! url. empty ());  
  
  string s;  s. reserve (1024);  // PAR  
  curl_easy_setopt (eh, CURLOPT_URL, url. c_str ());
  curl_easy_setopt (eh, CURLOPT_WRITEFUNCTION, write_string_cb);
  curl_easy_setopt (eh, CURLOPT_WRITEDATA, & s);
  EXEC_ASSERT (curl_easy_perform (eh) == 0);
  
  return s;
}

//



// #define URL "https://ftp.ncbi.nlm.nih.gov/pathogen/Technical/AMRFinder_technical/v2.data/"
#define URL "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/data/"

string getLatestVersion (Curl &curl)
// Return: empty() <=> failure
{
  string prevLine;
  const StringVector dir (curl. read (URL), '\n');
  if (verbose ())
    cout << dir << endl;
  for (const string& line : dir)
  {
    if (contains (line, "<a href=\"latest/\">latest/</a>"))
      break;
    prevLine = line;      
  }
  if (prevLine. empty ())
    return string ();
    
  const size_t pos1 = prevLine. find ('>');
  if (pos1 == string::npos)
    return string ();
  prevLine. erase (0, pos1 + 1);

  const size_t pos2 = prevLine. find ("/<");
  if (pos2 == string::npos)
    return string ();
  prevLine. erase (pos2);
  
  return prevLine;
}



void fetchAMRFile (Curl &curl,
                   const string &dir,
                   const string &fName) 
{
  ASSERT (isDirName (dir));
  ASSERT (! fName. empty ());  
  curl. download (string (URL "latest/") + fName, dir + fName);
}



	
// ThisApplication

struct ThisApplication : ShellApplication
{
  ThisApplication ()
    : ShellApplication ("Identify AMR genes in proteins and/or contigs and print a report", false, true, true)
    {
    	addKey ("database", "Directory for all versions of AMRFinder databases", "$BASE/data", 'd', "DATABASE_DIR");
    	  // Symbolic link ??
      addFlag ("quiet", "Suppress messages to STDERR", 'q');
	  #ifdef SVN_REV
	    version = SVN_REV;
	  #endif
    }



  void shellBody () const final
  {
          string mainDirOrig = getArg ("database");
    const bool   quiet       = getFlag ("quiet");
    
        
    Stderr stderr (quiet);
    stderr << "Running "<< getCommandLine () << '\n';
    const Verbose vrb (qc_on);
    
    // mainDir
    const Dir mainDir (mainDirOrig);
    string mainDirS (mainDir. get ());
    if (! isRight (mainDirS, "/"))
      mainDirS += "/";    

    findProg ("makeblastdb");
    findProg ("hmmpress");
    
    Curl curl;    
    
    
    const string latest_version (getLatestVersion (curl));
    if (latest_version. empty ())
      throw runtime_error ("Cannot get the latest version");
      
    const string latestDir (mainDirS + latest_version + "/");
    const string latestLink (mainDirS + "latest");
      
    if (! directoryExists (mainDirS))
      exec ("mkdir " + mainDirS);
    
    if (directoryExists (latestDir))
      stderr << latestDir << " already exists, overwriting what was there\n";
    else
      exec ("mkdir " + latestDir);
    
    if (directoryExists (latestLink))
      exec ("rm " + latestLink);
    exec ("ln -s " + latest_version + " " + latestLink);
    
    StringVector dnaPointMuts (ORGANISMS, '|');
    
    stderr << "Dowloading AMRFinder database version " << latest_version << " into " << latestDir << "\n";
    fetchAMRFile (curl, latestDir, "AMR.LIB");
    fetchAMRFile (curl, latestDir, "AMRProt");
    fetchAMRFile (curl, latestDir, "AMRProt-point_mut.tab");
    fetchAMRFile (curl, latestDir, "AMR_CDS");
    for (const string& dnaPointMut : dnaPointMuts)
    {
      fetchAMRFile (curl, latestDir, "AMR_DNA-" + dnaPointMut);
      fetchAMRFile (curl, latestDir, "AMR_DNA-" + dnaPointMut + ".tab");
    }
    fetchAMRFile (curl, latestDir, "fam.tab");
    fetchAMRFile (curl, latestDir, "changes.txt");
    
    stderr << "Indexing" << "\n";
    exec (fullProg ("hmmpress") + " -f " + latestDir + "AMR.LIB > /dev/null 2> /dev/null");
	  exec (fullProg ("makeblastdb") + " -in " + latestDir + "AMRProt  -dbtype prot  -logfile /dev/null");  
	  exec (fullProg ("makeblastdb") + " -in " + latestDir + "AMR_CDS  -dbtype nucl  -logfile /dev/null");  
    for (const string& dnaPointMut : dnaPointMuts)
  	  exec (fullProg ("makeblastdb") + " -in " + latestDir + "AMR_DNA-" + dnaPointMut + "  -dbtype nucl  -logfile /dev/null");  
  }
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



