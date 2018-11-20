// gff.cpp

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
*   .gff file reader
*
*/
   

#undef NDEBUG
#include "common.inc"

#include "gff.hpp"




namespace GFF_sp
{



// Cds

Cds::Cds (const string &contig_arg,
				  size_t start_arg,
				  size_t stop_arg,
				  bool strand_arg,
          size_t crossOriginSeqLen_arg)
: contig (contig_arg)
, start (start_arg)
, stop (stop_arg)
, strand (strand_arg)
, crossOriginSeqLen (crossOriginSeqLen_arg)
{ 
	ASSERT (! contig. empty ());
	if (crossOriginSeqLen)
	{
		swap (start, stop);
		start--;
		stop++;
	}
  ASSERT (start < stop); 
}
  
  
    
bool Cds::operator< (const Cds& other) const
{ 
	LESS_PART (*this, other, contig)
  LESS_PART (*this, other, start)
  LESS_PART (*this, other, stop)
  LESS_PART (*this, other, strand)
  LESS_PART (*this, other, crossOriginSeqLen);
  return false;
}




// Gff

Gff::Gff (const string &fName,
	        bool locus_tag)
{
	if (fName. empty ())
		throw runtime_error ("Empty GFF file name");
	
  LineInput f (fName /*, 100 * 1024, 1*/);
  while (f. nextLine ())
  {
    trim (f. line);
    if (   f. line. empty () 
        || f. line [0] == '#'
       )
      continue;

    replace (f. line, ' ', '_');  // to use '\t' as delimiter

   	const string errorS ("File " + fName + ", line " + toString (f. lineNum) + ": ");

    string seqid, source, type, startS, stopS, score /*real number*/, strand, phase /*frame*/, attributes;
    {
	    istringstream iss (f. line);
	    iss >> seqid >> source >> type >> startS >> stopS >> score >> strand >> phase >> attributes;
	  }

    if (attributes. empty ())
    	throw runtime_error (errorS + "9 fields are expected in each line");

    if (contains (seqid, ":"))
      findSplit (seqid, ':');  // = project_id
    if (seqid. empty ())
    	throw runtime_error (errorS + "empty sequence indentifier");
	  for (const char c : seqid)
	  	if (! printable (c))
	  		throw runtime_error (errorS + "Non-printable character in the sequence identifier: " + c);

    if (   type != "CDS"
        && type != "gene"
        && type != "pseudogene"
       )
      continue;
      
    long start = -1;
    try { start = str2<long> (startS); }
      catch (...) 
      {
	    	throw runtime_error (errorS + "Cannot read start");
      }
    if (start <= 0)
    	throw runtime_error (errorS + "start should be >= 1");
      
    long stop = -1;
    try { stop = str2<long> (stopS); }
      catch (...) 
      {
	    	throw runtime_error (errorS + "Cannot read stop");
      }
    if (stop <= 0)
    	throw runtime_error (errorS + "stop should be >= 1");

    if (start > stop)
    	throw runtime_error (errorS + "start cannot be greater than stop");

    start--;    	
    	
    if (   strand != "+" 
        && strand != "-"
       )
    	throw runtime_error (errorS + "strand should be '+' or '-'");
           
    const bool pseudo =    contains (attributes, "pseudo=true")
                        || contains (attributes, "gene_biotype=pseudogene")
                        || type == "pseudogene";
  //if (pseudo && type == "CDS")  
    //continue;

    string locusTag;
    const string locusTagName (locus_tag || pseudo ? "locus_tag=" : "Name=");
    while (! attributes. empty ())
    {
	    locusTag = findSplit (attributes, ';');
	    while (trimPrefix (locusTag, "_"));  // trim leading spaces
	    if (isLeft (locusTag, locusTagName))
	      break;
	  }
    if (! isLeft (locusTag, locusTagName))
    	continue;
    //throw runtime_error (errorS + "No attribute '" + locusTagName + "': " + f. line);
	  if (contains (locusTag, ":"))
	    { EXEC_ASSERT (isLeft (findSplit (locusTag, ':'), locusTagName)); }
	  else
	    findSplit (locusTag, '='); 
	  trimPrefix (locusTag, "\"");
	  trimSuffix (locusTag, "\"");
	  
    seqid2cdss [locusTag] << Cds (seqid, (size_t) start, (size_t) stop, strand == "+", 0);
  }
}
  
  

}
