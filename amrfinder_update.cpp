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
*               mkdir, ln, rm
*
* Release changes:
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



#if 0
  #define URL "https://ftp.ncbi.nlm.nih.gov/pathogen/Technical/AMRFinder_technical/test_database/"
#else
  #define URL "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/"
#endif



string getLatestMinor (Curl &curl)
// Return: empty() <=> failure
{
  StringVector dir (curl. read (URL), '\n');
  if (verbose ())
    cout << dir << endl;
    
  Vector<SoftwareVersion> vers;  
  for (string& line : dir)
    if (isLeft (line, "<a href="))
    {
      const size_t pos1 = line. find ('>');
      QC_ASSERT (pos1 != string::npos);
      line. erase (0, pos1 + 1);

      const size_t pos2 = line. find ("/<");
      QC_ASSERT (pos2 != string::npos);
      line. erase (pos2);
      
  	  istringstream iss (line);
  	  try 
  	  {
  		  SoftwareVersion ver (iss, true);
  		  vers << move (ver);
  		}
  		catch (...) {}
    }
  if (vers. empty ())
    return string ();
    
  vers. sort ();
  return vers. back (). getMinor ();
}



string getLatestDataVersion (Curl &curl,
                             const string &minor)
// Return: empty() <=> failure
{
  StringVector dir (curl. read (URL + minor + "/"), '\n');
  if (verbose ())
    cout << dir << endl;
    
  Vector<DataVersion> dataVersions;  
  for (string& line : dir)
    if (isLeft (line, "<a href="))
    {
      const size_t pos1 = line. find ('>');
      QC_ASSERT (pos1 != string::npos);
      line. erase (0, pos1 + 1);

      const size_t pos2 = line. find ("/<");
      QC_ASSERT (pos2 != string::npos);
      line. erase (pos2);
      
  	  istringstream iss (line);
  	  try 
  	  {
  		  DataVersion dv (iss);
  		  dataVersions << move (dv);
  		}
  		catch (...) {}
    }
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
    : ShellApplication ("Update the data for AMRFinder from " URL "\n\
Requirements:\n\
- the data/ directory contains subdirectories named by \"minor\" software versions (i.e., <major>.<minor>/);\n\
- the \"minor\" directories contain subdirectories named by data versions.\
", false, true, true)
    {
    	addKey ("database", "Directory for all versions of AMRFinder databases", "$BASE/data", 'd', "DATABASE_DIR");
    	  // Symbolic link ??
      addFlag ("quiet", "Suppress messages to STDERR", 'q');
	    version = SVN_REV;

      // curMinor
      {
    	  istringstream versionIss (version);
    		const SoftwareVersion softwareVersion (versionIss);
        curMinor = softwareVersion. getMinor ();
      }
    }



  void shellBody () const final
  {
    const string mainDirOrig = getArg ("database");
    const bool   quiet       = getFlag ("quiet");
    
        
    Stderr stderr (quiet);
    stderr << "Running "<< getCommandLine () << '\n';
    const Verbose vrb (qc_on);
    
    string mainDirS;
    {
      const Dir mainDir (mainDirOrig);
      mainDirS = mainDir. get ();
    }
    if (! isRight (mainDirS, "/"))
      mainDirS += "/";    

    findProg ("makeblastdb");
    findProg ("hmmpress");
    
    Curl curl;    
    
    
    // FTP site files
    const string latest_minor (getLatestMinor (curl));
    if (latest_minor. empty ())
      throw runtime_error ("Cannot get the latest software version");
    
    const string latest_version (getLatestDataVersion (curl, curMinor));
    if (latest_version. empty ())
      throw runtime_error ("Cannot get the latest data version for the current software");
      
    const string cur_latest_version (getLatestDataVersion (curl, latest_minor));
    if (cur_latest_version. empty ())
      throw runtime_error ("Cannot get the latest data version for the latest software (" + latest_minor + ")");

    if (latest_version != cur_latest_version)     
      stderr << "\nWARNING: A newer version of the database exists (" << cur_latest_version << "), but it requires "
                "a newer version of the software (" << latest_minor << ") to install.\n"
                "See https://github.com/ncbi/amr/wiki/Upgrading for more information.\n\n";
                      
    
    // Users's files  
    if (! directoryExists (mainDirS))
      exec ("mkdir -p " + mainDirS);
    
    const string latestDir (mainDirS + latest_version + "/");
    if (directoryExists (latestDir))
      stderr << latestDir << " already exists, overwriting what was there\n";
    else
      exec ("mkdir -p " + latestDir);

    {    
      const string latestLink (mainDirS + "latest");
      if (directoryExists (latestLink))
        exec ("rm " + latestLink);
      exec ("ln -s " + latest_version + " " + latestLink);
    }
    
    
    stderr << "Downloading AMRFinder database version " << latest_version << " into " << latestDir << "\n";
    const string urlDir (URL + curMinor + "/" + latest_version + "/");
    fetchAMRFile (curl, urlDir, latestDir, "AMR.LIB");
    fetchAMRFile (curl, urlDir, latestDir, "AMRProt");
    fetchAMRFile (curl, urlDir, latestDir, "AMRProt-mutation.tab");
    fetchAMRFile (curl, urlDir, latestDir, "AMRProt-suppress");
    fetchAMRFile (curl, urlDir, latestDir, "AMR_CDS");
    fetchAMRFile (curl, urlDir, latestDir, "database_format_version.txt");  // PD-3051 
    fetchAMRFile (curl, urlDir, latestDir, "fam.tab");
    fetchAMRFile (curl, urlDir, latestDir, "taxgroup.tab");
    fetchAMRFile (curl, urlDir, latestDir, "version.txt");
    
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



