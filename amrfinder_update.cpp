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
* Dependencies: NCBI BLAST, HMMer
*               ln
*               curl.{h,c}
*
* Release changes: see amrfinder.cpp
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



string curMinor;




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
  
  {
    OFStream f (fName);
    curl_easy_setopt (eh, CURLOPT_URL, url. c_str ());
    curl_easy_setopt (eh, CURLOPT_WRITEFUNCTION, write_stream_cb);
    curl_easy_setopt (eh, CURLOPT_WRITEDATA, & f);
    EXEC_ASSERT (curl_easy_perform (eh) == 0);
  }
  
  ifstream f (fName);
  string s;
  f >> s;
  if (s == "<?xml")
    throw runtime_error ("Cannot download " + strQuote (fName));
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



#undef HTTPS  // Otherwise: ftp


#if 0
  #define URL "https://ftp.ncbi.nlm.nih.gov/pathogen/Technical/AMRFinder_technical/test_database/"
#else
  #ifdef HTTPS
    #define URL "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/"  
  #else
    #define URL "ftp://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/"  
  #endif
#endif



string getLatestMinor (Curl &curl)
// Return: empty() <=> failure
{
  StringVector dir (curl. read (URL), '\n', true);
  if (verbose ())
  {
    save (cout, dir, '\n'); 
    cout << endl;
  }
    
  Vector<SoftwareVersion> vers;  
  for (string& line : dir)
  #ifdef HTTPS
    if (isLeft (line, "<a href="))
  	  try 
  	  {
        const size_t pos1 = line. find ('>');
        QC_ASSERT (pos1 != string::npos);
        line. erase (0, pos1 + 1);

        const size_t pos2 = line. find ("/<");
        QC_ASSERT (pos2 != string::npos);
        line. erase (pos2);
        
    	  istringstream iss (line);
  		  SoftwareVersion ver (iss, true);
  		  vers << move (ver);
  		}
  		catch (...) {}
  #else
    if (! contains (line, " -> "))
    {
      trimTrailing (line);
      const size_t pos = line. rfind (' ');
      if (pos != string::npos)
      {
    	  istringstream iss (line. substr (pos + 1));
    	  try 
    	  {
    		  SoftwareVersion ver (iss, true);
    		  vers << move (ver);
    		}
  		  catch (...) {}    		
      }
    }
  #endif
  if (vers. empty ())
    return string ();
    
  vers. sort ();
  return vers. back (). getMinor ();
}



string getLatestDataVersion (Curl &curl,
                             const string &minor)
// Return: empty() <=> failure
{
  StringVector dir (curl. read (URL + minor + "/"), '\n', true);
  if (verbose ())
  {
    save (cout, dir, '\n');
    cout << endl;
  }
    
  Vector<DataVersion> dataVersions;  
  for (string& line : dir)
  #ifdef HTTPS
    if (isLeft (line, "<a href="))
      try
      {
        const size_t pos1 = line. find ('>');
        QC_ASSERT (pos1 != string::npos);
        line. erase (0, pos1 + 1);

        const size_t pos2 = line. find ("/<");
        QC_ASSERT (pos2 != string::npos);
        line. erase (pos2);
        
    	  istringstream iss (line);
  		  DataVersion dv (iss);
  		  dataVersions << move (dv);
    	}
  		catch (...) {}
  #else
    if (! contains (line, " -> "))
    {
      trimTrailing (line);
      const size_t pos = line. rfind (' ');
      if (pos != string::npos)
      {
    	  istringstream iss (line. substr (pos + 1));
    	  try 
    	  {
    		  DataVersion dv (iss);
    		  dataVersions << move (dv);
    		}
  		  catch (...) {}    		
      }
    }
  #endif
  if (dataVersions. empty ())
    return string ();
    
  dataVersions. sort ();
  return dataVersions. back (). str ();
}



void fetchAMRFile (Curl &curl,
                   const string &urlDir,
                   const string &localDir,
                   const string &fName) 
{
  ASSERT (isDirName (urlDir));
  ASSERT (isDirName (localDir));
  ASSERT (! fName. empty ());  
  curl. download (urlDir + fName, localDir + fName);
}



	
// ThisApplication

struct ThisApplication : ShellApplication
{
  ThisApplication ()
    : ShellApplication ("Update the database for AMRFinder from " URL "\n\
Requirement: the database directory contains subdirectories named by database versions.\
", false, true, true)
    {
    	addKey ("database", "Directory for all versions of AMRFinder databases", "$BASE/data", 'd', "DATABASE_DIR");
    	  // Symbolic link ??
    	addFlag ("force_update", "Force updating the AMRFinder database");  // PD-3469
      addFlag ("quiet", "Suppress messages to STDERR", 'q');
	    version = SVN_REV;

      // curMinor
      {
    	  istringstream versionIss (version);
    		const SoftwareVersion softwareVersion (versionIss);
        curMinor = softwareVersion. getMinor ();
      }
    }



  void createLatestLink (const string &mainDirS,
                         const string &latestDir) const
  {   
    ASSERT (! mainDirS. empty ()); 
    ASSERT (! latestDir. empty ()); 
    const string latestLink (mainDirS + "latest");
    if (directoryExists (latestLink))
      removeFile (latestLink);
    exec ("ln -s " + shellQuote (path2canonical (latestDir)) + " " + shellQuote (latestLink));
  }



  void shellBody () const final
  {
    const string mainDirOrig  = getArg ("database");
    const bool   force_update = getFlag ("force_update");
    const bool   quiet        = getFlag ("quiet");
    
        
    Stderr stderr (quiet);
    stderr << "Running: "<< getCommandLine () << '\n';
	//stderr << "Current software minor version: " << curMinor << '\n'; 
    const Verbose vrb (qc_on);
    

    Curl curl;    
        
    
    stderr << "Looking up databases at " << URL << '\n';

    // FTP site files
    const string latest_minor (getLatestMinor (curl));
    if (latest_minor. empty ())
      throw runtime_error ("Cannot get the latest software minor version");
  //stderr << "Latest software minor version: " << latest_minor << "\n";
    
    const string latest_data_version (getLatestDataVersion (curl, curMinor));
    if (latest_data_version. empty ())
      throw runtime_error ("Cannot get the latest database version for the current software");
  //stderr << "Latest database version: " << latest_data_version << "\n";
      
    const string cur_latest_data_version (getLatestDataVersion (curl, latest_minor));
    if (cur_latest_data_version. empty ())
      throw runtime_error ("Cannot get the latest database version for the latest software (" + latest_minor + ")");

    if (latest_data_version != cur_latest_data_version)  
    {   
      stderr << "\n";
      const Warning w (stderr);
      stderr << "A newer version of the database exists (" << cur_latest_data_version << "), but it requires "
                "a newer version of the software (" << latest_minor << ") to install.\n"
                "See https://github.com/ncbi/amr/wiki/Upgrading for more information.\n\n";
    }
                      
    
    findProg ("makeblastdb");
    findProg ("hmmpress");
    

    // Users's files  
    string mainDirS;
    {
      const Dir mainDir (mainDirOrig);
      mainDirS = mainDir. get ();
    }
    if (! isRight (mainDirS, "/"))
      mainDirS += "/";    

  //Dir (mainDirS). create ();
    
    const string versionFName ("version.txt");
    const string urlDir (URL + curMinor + "/" + latest_data_version + "/");
    
    const string latestDir (mainDirS + latest_data_version + "/");
    if (directoryExists (latestDir))
    {
      if (force_update)
        stderr << shellQuote (latestDir) << " already exists, overwriting what was there\n";
      else
      {
        curl. download (urlDir + versionFName, tmp);
        const StringVector version_old (latestDir + versionFName, (size_t) 100, true);
        const StringVector version_new (tmp, (size_t) 100, true);
        if (   ! version_old. empty () 
            && ! version_new. empty ()
            && version_old. front () == version_new. front ()
           )
        {
          const Warning w (stderr);
          stderr << shellQuote (latestDir) << " contains the latest version: " << version_old. front () << '\n';
          stderr << "Skipping update, use amrfinder --force_update to overwrite the existing database\n";
          createLatestLink (mainDirS, latestDir);
          return;
        }
      }
    }
    else
      Dir (latestDir). create ();
    
    stderr << "Downloading AMRFinder database version " << latest_data_version << " into " << shellQuote (latestDir) << "\n";
    fetchAMRFile (curl, urlDir, latestDir, "AMR.LIB");
    fetchAMRFile (curl, urlDir, latestDir, "AMRProt");
    fetchAMRFile (curl, urlDir, latestDir, "AMRProt-mutation.tab");
    fetchAMRFile (curl, urlDir, latestDir, "AMRProt-suppress");
    fetchAMRFile (curl, urlDir, latestDir, "AMRProt-susceptible.tab");
    fetchAMRFile (curl, urlDir, latestDir, "AMR_CDS");
    fetchAMRFile (curl, urlDir, latestDir, "database_format_version.txt");  // PD-3051 
    fetchAMRFile (curl, urlDir, latestDir, "fam.tab");
    fetchAMRFile (curl, urlDir, latestDir, "taxgroup.tab");
    fetchAMRFile (curl, urlDir, latestDir, versionFName);
    
    StringVector dnaPointMuts;
    {
      LineInput f (latestDir + "taxgroup.tab");
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
    for (const string& dnaPointMut : dnaPointMuts)
    {
      fetchAMRFile (curl, urlDir, latestDir, "AMR_DNA-" + dnaPointMut);
      fetchAMRFile (curl, urlDir, latestDir, "AMR_DNA-" + dnaPointMut + ".tab");
    }

    fetchAMRFile (curl, urlDir, latestDir, "changes.txt");
    
    stderr << "Indexing" << "\n";
    exec (fullProg ("hmmpress") + " -f " + shellQuote (latestDir + "AMR.LIB") + " > /dev/null 2> /dev/null");
    exec ("ln -s " + shellQuote (path2canonical (latestDir)) + " " + tmp + ".db");
	  exec (fullProg ("makeblastdb") + " -in " + tmp + ".db/AMRProt" + "  -dbtype prot  -logfile /dev/null");  
	  exec (fullProg ("makeblastdb") + " -in " + tmp + ".db/AMR_CDS" + "  -dbtype nucl  -logfile /dev/null");  
    for (const string& dnaPointMut : dnaPointMuts)
  	  exec (fullProg ("makeblastdb") + " -in " + tmp + ".db/AMR_DNA-" + dnaPointMut + "  -dbtype nucl  -logfile /dev/null");

    createLatestLink (mainDirS, latestDir);
  }
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



