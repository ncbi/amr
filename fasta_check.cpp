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

#include "common.hpp"
using namespace Common_sp;

#include "common.inc"



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
      addFlag ("ambig", "Ambiguous characters are allowed");
      addKey ("ambig_max", "Max. number of ambiguous characters in sequences", "0");
      addFlag ("stop_codon", "Stop codons ('*') in the protein sequence are allowed");
      addKey ("len", "Output file with lines: <sequence id> <length>");
      addKey ("out", "Output FASTA file with some of the issues fixed");
	    version = SVN_REV;
    }



  void body () const final
  {
    const string fName     = getArg ("in");
    const bool aa          = getFlag ("aa");
    const bool hyphen      = getFlag ("hyphen");
    const bool ambig       = getFlag ("ambig");
    const size_t ambig_max = str2<size_t> (getArg ("ambig_max"));
    const bool stop_codon  = getFlag ("stop_codon");
    const string lenFName  = getArg ("len");
    const string outFName  = getArg ("out");
    
    QC_IMPLY (stop_codon, aa);
    

    unique_ptr<OFStream> lenF;
    if (! lenFName. empty ())
      lenF. reset (new OFStream (lenFName));
    unique_ptr<OFStream> outF;
    if (! outFName. empty ())
      outF. reset (new OFStream (outFName));
    size_t lines = 0;
    StringVector ids;  ids. reserve (100000);  // PAR
    size_t seqSize_max = 0;
    size_t seqSize_sum = 0;
  //string errorS;
    // One sequence
    size_t xs = 0;
    string header;
    string seq;
    
    auto processSeq = [&] () 
  	  {
  	    if (! lines)
  	      return;
 	  		ASSERT (! header. empty ());
 	  		ASSERT (! ids. empty ());
 	  		const string id (ids. back ());
  	    if (aa && ! stop_codon)
  	    {  	    
    	    while (! seq. empty () && seq. back () == '*')
   	  		  if (outF)
   	  		    seq. erase (seq. size () - 1);
   	  		  else
    		      throw runtime_error (id + ": '*' at the sequence end");
    		}
  	    if (seq. empty ())
 	  		  throw runtime_error (id + ": Empty sequence");
 	  		bool skip = false;
  	    if (! ambig && xs > ambig_max)
  	    {
 	  		  if (outF)
 	  		    skip = true;
 	  		  else
  		      throw runtime_error (id + ": Too many ambiguities");
  		  }
  		  if (skip)
  		    { LOG ("Skipping " + id); }
  		  else
  		  {
     	  	if (lenF. get ())
     	  	  *lenF << id << '\t' << seq. size () << endl;
  	      if (outF)
  	        *outF << header << endl << seq << endl;
   	  	  maximize (seqSize_max, seq. size ());
   	  	  seqSize_sum += seq. size ();
   	  	}
    	  xs = 0;
    	  header. clear ();
    	  seq. clear ();
 	    };
    
    size_t nuc = 0;   
    {
      LineInput f (fName); 
      string id;
      while (f. nextLine ())
      {
        trimTrailing (f. line);
        if (f. line. empty ())
        	continue;
      	const string errorS ("File " + fName + ", " + f. lineStr (false) + ": ");
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
      	  		throw runtime_error (errorS + "Non-printable character in the sequence identifier: " + to_string ((int) c));
      	  // BLAST: PD-4548
      	  if (! aa)
      	  {
        	  if (id. front () == '?')
       	  		throw runtime_error (errorS + "Sequence identifier starts with '?'");
       	  	for (const char c : {',', ';', '.', '~'})
          	  if (id. back () == c)
         	  		throw runtime_error (errorS + "Sequence identifier ends with " + strQuote (string (1, c)));
        	  if (contains (id, "\\t"))
       	  		throw runtime_error (errorS + "Sequence identifier contains '\\t'");
        	  if (contains (id, ",,"))
       	  		throw runtime_error (errorS + "Sequence identifier contains ',,'");
       	  }
     	    processSeq ();
      	  header = f. line;
      	  ids << id;
      	}
      	else 
      	{
      		if (! lines)
      			throw runtime_error (errorS + "FASTA should start with '>'");
      	  for (const char c : f. line)
      	  {
      	    bool skip = false;
      	  	if (c == '-')
      	  		if (hyphen)
      	  			;
      	  		else
      	  		{
      	  		  if (outF)
      	  		    skip = true;
      	  		  else
  	    	  		  throw runtime_error (errorS + "Hyphen in the sequence");  	    	  		  
  	    	    }
  	    	  else
  	    	  {
  	    	  	const char c1 = toLower (c);
  	    	  	if (aa)
  	    	  	{
  		    	  	if (! charInSet (c1, "acdefghiklmnpqrstvwyxbzjuoacdefghiklmnpqrstvwyxbzjuo*"))
  		    	  		throw runtime_error (errorS + "Wrong amino acid character: (code = " + to_string ((int) c) + ") '" + c + "'");
  		    	    if (charInSet (c1, "acgt"))
  		    	    	nuc++;
  		    	    if (charInSet (c1, "xbzjuo"))
  		    	      xs++;
  		    	  }
  	    	  	else
  	    	  	{
  		    	  	if (! charInSet (c1, "acgtbdhkmnrsvwyacgtbdhkmnrsvwy"))
  		    	  		throw runtime_error (errorS + "Wrong nucleotide character: (code = " + to_string ((int) c) + ") '" + c + "'");
  		    	    if (charInSet (c1, "bdhkmnrsvwy"))
  		    	      xs++;
  		    	  }
  		    	}
  		    	if (! skip)
  		    	  seq += c;
  		    }
      	}
      	lines++;
  	  }
  	}
    processSeq ();	// Last sequence
	  if (! lines)
	  	throw runtime_error ("Empty file"); 
  	if (aa && (double) nuc / (double) seqSize_sum > 0.9)  // PAR
  		throw runtime_error ("Protein sequences looks like a nucleotide sequences");
  		
	  ids. sort ();
	  const size_t index = ids. findDuplicate ();
	  if (index != no_index)
	  	throw runtime_error ("Duplicate identifier: " + ids [index]);
	  	
	  cout << ids. size () << endl
	       << seqSize_max << endl
	       << seqSize_sum << endl;
  }
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



