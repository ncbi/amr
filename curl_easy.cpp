// curl_easy.cpp

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
*   curl_easy functions
*
*/


#undef NDEBUG

#include "curl_easy.hpp"
using namespace Common_sp;

#include "common.inc"



namespace CURL_sp
{
  
  
  
SoftwareVersion getLibVersion ()
{
  if (const curl_version_info_data* ver = curl_version_info (CURLVERSION_NOW))
  {
    const uint major = (ver->version_num >> 16) & 0xff;
    const uint minor = (ver->version_num >> 8)  & 0xff;
    const uint patch = ver->version_num         & 0xff;
    return SoftwareVersion (major, minor, patch);
  }
  return SoftwareVersion ();
}


  

// Curl

namespace 
{
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
      f << ptr [i];
    
    return nMemb;
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
      s += ptr [i];
    
    return nMemb;
  }
}


 	
void Curl::download (const string &url,
                     const string &fName) 
{
  ASSERT (! fName. empty ());  
  
  {
    OFStream f (fName);
    curl_easy_setopt (eh, CURLOPT_WRITEFUNCTION, write_stream_cb);
    curl_easy_setopt (eh, CURLOPT_WRITEDATA, & f);
    process (url, "download");
  }
  
  ifstream f (fName);
  string s;
  f >> s;
  if (s == "<?xml")
    throw runtime_error ("Cannot download " + strQuote (fName));
}



string Curl::read (const string &url)
{
  string s;  s. reserve (1024);  // PAR  
  curl_easy_setopt (eh, CURLOPT_WRITEFUNCTION, write_string_cb);
  curl_easy_setopt (eh, CURLOPT_WRITEDATA, & s);
  process (url, "read");
  
  return s;
}



void Curl::process (const string &url,
                    const string &error_msg_action)
{
  QC_ASSERT (! url. empty ());  

  char err [CURL_ERROR_SIZE + 1] = "";
  curl_easy_setopt (eh, CURLOPT_ERRORBUFFER, err);
  
  curl_easy_setopt (eh, CURLOPT_URL, url. c_str ());
  if (isLeft (url, "ftp://"))
    curl_easy_setopt (eh, CURLOPT_FTP_USE_EPSV, 0);

  const CURLcode cc = curl_easy_perform (eh);
  if (cc)
  {
    const SoftwareVersion ver (CURL_sp::getLibVersion ());
    throw runtime_error ("CURL: Cannot " + error_msg_action
                         + "\n  from " + url 
                         + "\n  code=" + to_string (cc) 
                         + "\n  error: " + err
                         + "\n  version: " + ver. str ()
                         );
  }
}




}  // namespace


