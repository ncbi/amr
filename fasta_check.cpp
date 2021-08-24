// fasta_check.cpp

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
*   Check the correctness of a FASTA file
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
    : Application ("Check the correctness of a FASTA file. Exit with an error if it is incorrect. Print the number of sequences, max. sequence length and total sequence length")
    {
      addPositional ("in", "FASTA file");
      addFlag ("aa", "Amino acid sequenes, otherwise nucleotide");
      addFlag ("hyphen", "Hyphens are allowed");
      addKey ("len", "Output file with lines: <sequence id> <length>");
	    version = SVN_REV;
    }



  void body () const final
  {
    const string fName    = getArg ("in");
    const bool aa         = getFlag ("aa");
    const bool hyphen     = getFlag ("hyphen");
    const string lenFName = getArg ("len");
    

    unique_ptr<OFStream> lenF;
    if (! lenFName. empty ())
      lenF. reset (new OFStream (lenFName));
    size_t lines = 0;
    bool first = true;
    StringVector ids;  ids. reserve (100000);  // PAR
    size_t seqSize_max = 0;
    size_t seqSize = 0;
    size_t allSize = 0;
    size_t nuc = 0;   
    {
      LineInput f (fName); 
      string id;
      string errorS;
      while (f. nextLine ())
      {
        trimTrailing (f. line);
        if (f. line. empty ())
        	continue;
      	errorS = "File " + fName + ", line " + toString (f. lineNum) + ": ";
      	lines++;
      	if (f. line [0] == '>')
      	{
      		size_t pos = 1;
      		while (pos < f. line. size () && ! isspace (f. line [pos]))
      		  pos++;
      		id = f. line. substr (1, pos - 1);
      		if (id. empty ())
      			throw runtime_error (errorS + "Empty sequence identifier");
        #if 0
      		if (id. size () > 1000)  // PAR
      			throw runtime_error (errorS + "Too long sequence identifier");
        #endif
      	  for (const char c : id)
      	  	if (! printable (c))
      	  		throw runtime_error (errorS + "Non-printable character in the sequence identifier: " + c);
      	  if (! first && seqSize == 0)
     	  		throw runtime_error (errorS + "Empty sequence");
     	  	if (lenF. get () && ! ids. empty ())
     	  	  *lenF << ids. back () << '\t' << seqSize << endl;
      	  ids << id;
   	  	  maximize (seqSize_max, seqSize);
      	  seqSize = 0;
      	}
      	else 
      	{
      		if (first)
      			throw runtime_error (errorS + "FASTA should start with '>'");
      		seqSize += f. line. size ();
      		allSize += f. line. size ();
      	  for (const char c : f. line)
      	  	if (c == '-')
      	  		if (hyphen)
      	  			;
      	  		else
  	    	  		throw runtime_error (errorS + "hyphen in the sequence");
  	    	  else
  	    	  {
  	    	  	const char c1 = toLower (c);
  	    	  	if (aa)
  	    	  	{
  		    	  	if (! charInSet (c1, "acdefghiklmnpqrstvwyxbzjuoacdefghiklmnpqrstvwyxbzjuo*"))
  		    	  		throw runtime_error (errorS + "Wrong amino acid character: '" + c + "'");
  		    	    if (charInSet (c1, "acgt"))
  		    	    	nuc++;
  		    	  }
  	    	  	else
  		    	  	if (! charInSet (c1, "acgtbdhkmnrsvwyacgtbdhkmnrsvwy"))
  		    	  		throw runtime_error (errorS + "Wrong nucleotide character: '" + c + "'");
  		    	}
      	}
      	first = false;
  	  }
  	}
  	if (lenF. get () && ! ids. empty ())
  	  *lenF << ids. back () << '\t' << seqSize << endl;
 	  maximize (seqSize_max, seqSize);

	  if (! lines)
	  	throw runtime_error ("Empty file"); 
	  if (! first && seqSize == 0)
  		throw runtime_error ("Empty sequence");
  	if (aa && (double) nuc / (double) allSize > 0.9)  // PAR
  		throw runtime_error ("Protein sequence looks like a nucleotide sequence");
  		
	  ids. sort ();
	  const size_t index = ids. findDuplicate ();
	  if (index != no_index)
	  	throw runtime_error ("Duplicate identifier: " + ids [index]);
	  	
	  cout << ids. size () << endl
	       << seqSize_max << endl
	       << allSize << endl;
  }
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



