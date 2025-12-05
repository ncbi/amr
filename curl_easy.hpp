// curl_easy.hpp

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
* Dependencies: curl.{h,c}
*
* File Description:
*   curl_easy functions
*
*/


#ifdef _MSC_VER
  #error "UNIX is required"
#endif


#include <unistd.h>
extern "C" 
{
  #include <curl/curl.h>
    // Linking requires:  -lcurl
}

#include "common.hpp"
using namespace Common_sp;




namespace CURL_sp
{
  
  
SoftwareVersion getLibVersion ();

  

struct Curl
{
  CURL* eh {nullptr};


  Curl ()
    : eh (curl_easy_init ())
    { if (! eh)
        throw runtime_error ("Cannot initialize curl_easy");        
      // Override the libcurl system-wide default
      // PD-5495 / https://github.com/ncbi/amr/issues/170
      if (const char *env_ca_bundle = getenv ("CURL_CA_BUNDLE"))          
        curl_easy_setopt (eh, CURLOPT_CAINFO, env_ca_bundle);  
    }
 ~Curl ()
   { curl_easy_cleanup (eh); }


  void download (const string &url,
                 const string &fName);
  string read (const string &url);
private:
  void process (const string &url,
                const string &error_msg_action);
};

	


}  // namespace


