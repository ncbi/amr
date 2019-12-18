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



// Locus

Locus::Locus (size_t lineNum_arg,
              const string &contig_arg,
    				  size_t start_arg,
    				  size_t stop_arg,
    				  bool strand_arg,
              bool partial_arg,
              size_t crossOriginSeqLen_arg)
: lineNum (lineNum_arg)
, contig (contig_arg)
, start (start_arg)
, stop (stop_arg)
, strand (strand_arg)
, partial (partial_arg)
, contigLen (crossOriginSeqLen_arg)
, crossOrigin (bool (crossOriginSeqLen_arg))
{ 
//QC_ASSERT (lineNum >= 1);
	trim (contig, '_');
	if (contig. empty ())
		throw runtime_error ("Empty contig name");
	if (crossOrigin)
	{
		swap (start, stop);
		start--;
		stop++;
	  QC_ASSERT (contigLen);
		QC_ASSERT (stop <= contigLen);
	}
  QC_ASSERT (start < stop); 
}
  
  
    
bool Locus::operator< (const Locus& other) const
{ 
	LESS_PART (*this, other, contig)
  LESS_PART (*this, other, start)
  LESS_PART (*this, other, stop)
  LESS_PART (*this, other, strand)
  LESS_PART (*this, other, contigLen);
  LESS_PART (*this, other, crossOrigin);
  return false;
}




// Gff

Annot::Annot (Gff,
	            const string &fName,
	            bool trimProject,
	            bool locus_tag,
	            bool pgap)
{
  IMPLY (pgap, ! locus_tag);
  IMPLY (pgap, ! trimProject);

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

    constexpr char tmpSpace {'_'};
    replace (f. line, ' ', tmpSpace);  // to use '\t' as delimiter

   	const string errorS ("File " + fName + ", line " + toString (f. lineNum) + ": ");

    string contig, source, type, startS, stopS, score /*real number*/, strand, phase /*frame*/, attributes;
    static Istringstream iss;
    iss. reset (f. line);
    iss >> contig >> source >> type >> startS >> stopS >> score >> strand >> phase >> attributes;
    trim (contig, tmpSpace);
    trim (source, tmpSpace);
    trim (type, tmpSpace);
    trim (startS, tmpSpace);
    trim (stopS, tmpSpace);
    trim (score, tmpSpace);
    trim (strand, tmpSpace);
    trim (phase, tmpSpace);
    trim (attributes, tmpSpace);

    if (attributes. empty ())
    	throw runtime_error (errorS + "9 fields are expected in each line");

    if (trimProject)
      if (contains (contig, ":"))
        findSplit (contig, ':');  // = project_id
    if (contig. empty ())
    	throw runtime_error (errorS + "empty sequence indentifier");
	  for (const char c : contig)
	  	if (! printable (c))
	  		throw runtime_error (errorS + "Non-printable character in the sequence identifier: " + c);

    if (   type != "CDS"
        && type != "gene"
        && type != "pseudogene"
       )
      continue;
      
    if (pgap && type != "CDS")
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
    
    const bool partial = contains (attributes, "partial=true");

    string locusTag;
    const string locusTagName (! pgap && (locus_tag || pseudo) ? "locus_tag=" : "Name=");
    while (! attributes. empty ())
    {
	    locusTag = findSplit (attributes, ';');
  	  trim (locusTag, tmpSpace);
	    if (isLeft (locusTag, locusTagName))
	      break;
	  }
    if (! isLeft (locusTag, locusTagName))
    	continue;
    //throw runtime_error (errorS + "No attribute '" + locusTagName + "': " + f. line);
	  if (! pgap && contains (locusTag, ":"))
	    { EXEC_ASSERT (isLeft (findSplit (locusTag, ':'), locusTagName)); }
	  else
	    findSplit (locusTag, '='); 
	  trimPrefix (locusTag, "\"");
	  trimSuffix (locusTag, "\"");
	  trim (locusTag, tmpSpace);
	  if (pgap)
	  {
	    size_t pos = locusTag. rfind (':');
	    if (pos != string::npos)
	    {
	      locusTag [pos] = '|';
	      locusTag = "gnl|" + locusTag;
	    }
	    pos = contig. rfind (':');
	    if (pos != string::npos)
	      contig [pos] = '|';
	  }
	  QC_ASSERT (! locusTag. empty ());
	  
	  const Locus locus (f. lineNum, contig, (size_t) start, (size_t) stop, strand == "+", partial, 0);
	#if 0
	  // DNA may be truncated
    if (type == "CDS" && ! pseudo && locus. size () % 3 != 0)
    {
      cout << "Locus tag: " << locusTag << endl;
      locus. print (cout);
      ERROR;
    }
  #endif

    prot2cdss [locusTag] << locus;
  }
}
  
  

Annot::Annot (Bed,
	            const string &fName)
{
	if (fName. empty ())
		throw runtime_error ("Empty BED file name");
	
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

    string contig, locusTag;
    size_t start, stop;
    double score;
    char strand = ' ';
    static Istringstream iss;
    iss. reset (f. line);
    iss >> contig >> start >> stop >> locusTag >> score >> strand;

    if (strand == ' ')
    	throw runtime_error (errorS + "at least 5 fields are expected in each line");

	  for (const char c : contig)
	  	if (! printable (c))
	  		throw runtime_error (errorS + "Non-printable character in the sequence identifier: " + c);

    if (start >= stop)
    	throw runtime_error (errorS + "start should be less than stop");

    if (   strand != '+'
        && strand != '-'
       )
    	throw runtime_error (errorS + "strand should be '+' or '-'");
           
	  trim (locusTag, '_');
	  ASSERT (! locusTag. empty ());
    prot2cdss [locusTag] << Locus (f. lineNum, contig, start, stop, strand == '+', false/*partial*/, 0);
  }
}
  
  

}
