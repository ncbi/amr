// fasta2parts.cpp

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
*   Split the sequences a FASTA file into chunks without breaking sequences
*
*/
   
   
#undef NDEBUG 
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;



namespace 
{



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Split the sequences a FASTA file into parts without breaking sequences")
    {
      addPositional ("in", "FASTA file");
      addPositional ("parts_max", "Max. number of parts (>= 2)");
      addPositional ("dir", "Output directory where chunks are saved named by integers starting with 1");
    }



  void body () const final
  {
    const string fName     =               getArg ("in");
    const size_t parts_max = str2<size_t> (getArg ("parts_max"));
    const string dirName   =               getArg ("dir");
    
    if (parts_max <= 1)
      throw runtime_error ("Number of parts must be >= 2");
      

    const size_t chunk_min = (size_t) getFileSize (fName) / parts_max + 1;

    size_t part = 0;
    OFStream* out = nullptr;
    size_t seqSize = 0;
    LineInput f (fName); 
    while (f. nextLine ())
    {
      if (f. line. empty ())
      	continue;
    	if (   f. line [0] == '>'
    	    && seqSize >= chunk_min
    	   )
    	{
    	  delete out;
    	  out = nullptr;
    	  seqSize = 0;
    	}
    	if (! out)
    	{
    	  part++;
    	  ASSERT (part <= parts_max);
    	  out = new OFStream (dirName, toString (part), "");
    	}
   		*out << f. line << endl;
   		seqSize += f. line. size ();
	  }
	  delete out;
  }
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



