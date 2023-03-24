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
              size_t crossOriginSeqLen_arg,
              string gene_arg,
              string product_arg)
: lineNum (lineNum_arg)
, contig (contig_arg)
, start (start_arg)
, stop (stop_arg)
, strand (strand_arg)
, partial (partial_arg)
, contigLen (crossOriginSeqLen_arg)
, crossOrigin (bool (crossOriginSeqLen_arg))
, gene (gene_arg)
, product (product_arg)
{ 
//QC_ASSERT (lineNum >= 1);
  trim (contig);
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
//LESS_PART (*this, other, contigLen);
  LESS_PART (*this, other, crossOrigin);
  return false;
}




// Gff

const StringVector Gff::names {"bakta", "genbank", "microscope", "patric", "pgap", "prokka", "pseudomonasdb", "rast"};


Gff::Type Gff::name2type (const string &name)
{
  if (name == "bakta")          return bakta;
  if (name == "genbank")        return genbank;
  if (name == "microscope")     return microscope;
  if (name == "patric")         return patric;
  if (name == "pgap")           return pgap;
  if (name == "prokka")         return prokka;
  if (name == "pseudomonasdb")  return pseudomonasdb;
  if (name == "rast")           return rast;
  throw runtime_error ("Unknown GFF type: " + strQuote (name));
}




// Annot

namespace
{
  
void pgap_accession (string &accession,
                     bool lcl)
// Update: accession
{
  static const string gnlPrefix ("gnl|");
  static const string lclPrefix ("lcl|");
  
  size_t pos = accession. rfind (':');
  if (pos == string::npos)
  {
    if (lcl)
      accession = lclPrefix + accession;
  }
  else 
  {
    if (lcl)
      throw runtime_error ("Accession " + strQuote (accession) + " cannot have " + strQuote (gnlPrefix) + " and " + strQuote (lclPrefix) + " at the same time");
    accession [pos] = '|';
    accession = gnlPrefix + accession;
  }
}
    
}



Annot::Annot (const string &fName,
              Gff::Type gffType,
	            bool protMatch,
	            bool lcl)
{
  IMPLY (protMatch, gffType == Gff::genbank || gffType == Gff::microscope);
  IMPLY (gffType == Gff::microscope, protMatch);
//IMPLY (trimProject, gffType == Gff::genbank);
  IMPLY (lcl, gffType == Gff::pgap);

	if (fName. empty ())
		throw runtime_error ("Empty GFF file name");
	
  LineInput f (fName /*, 100 * 1024, 1*/);
  while (f. nextLine ())
  {
    trim (f. line);
    
    if (   (   gffType == Gff::prokka
            || gffType == Gff::bakta
           )
        && f. line == "##FASTA"
       )
      break;
    
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

  #if 0
    if (trimProject)
      if (contains (contig, ":"))
        findSplit (contig, ':');  // = project_id
  #endif
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
      
    if (gffType == Gff::pgap && type != "CDS")
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
    
    string protAttr = "Name";
    switch (gffType)
    {
      case Gff::bakta:         protAttr = "ID"; break;
      case Gff::genbank:       protAttr = (protMatch || pseudo) ? "locus_tag" : "Name"; break;
      case Gff::microscope:    protAttr = "ID"; break;
      case Gff::patric:        protAttr = "ID"; break;
      case Gff::prokka:        protAttr = "ID"; break;
      case Gff::pseudomonasdb: protAttr = "Alias"; break;  // for type = "gene", "locus" for type = "CDS"
      case Gff::rast:          protAttr = "ID"; break;
      default: break;
    }
    ASSERT (! protAttr. empty ());
    protAttr += "=";
        
    string prot;
    string gene;
    string product;
    {
      string locusTag;
      while (! attributes. empty ())
      {
  	    string s (findSplit (attributes, ';'));
    	  trim (s, tmpSpace);
  	    if (isLeft (s, protAttr))
  	    {
  	      prot = s;
          findSplit (prot, '='); 
  	    }
  	    else if (isLeft (s, "gene="))
  	    {
  	      gene = s;
  	      findSplit (gene, '=');
  	    }
  	    else if (isLeft (s, "product="))
  	    {
  	      product = s;
  	      findSplit (product, '=');
  	      replace (product, tmpSpace, ' ');
  	    }
  	    else if (gffType == Gff::patric && isLeft (s, "locus_tag="))
  	    {
  	      locusTag = s;
  	      findSplit (locusTag, '=');
  	    }
  	  }
  	  trimPrefix (prot, "\"");
  	  trimSuffix (prot, "\"");
      if (prot. empty ())
      	continue;
      //throw runtime_error (errorS + "no attribute '" + protAttr + "': " + f. line);
      switch (gffType)
      {
        case Gff::genbank:       if (contains (prot, ":"))  findSplit (prot, ':');  break;
        case Gff::patric:        if (! locusTag. empty ())  prot += "|" + locusTag;  
                                 if (isLeft (contig, "accn|"))  contig. erase (0, 5);
                                 break;
        default:  break;
      }
    }
	  trim (prot, tmpSpace);
	  if (gffType == Gff::pgap)
	  {
	    pgap_accession (prot, false);
	    pgap_accession (contig, lcl);
	  }
	  QC_ASSERT (! prot. empty ());
	  
	  Locus locus (f. lineNum, contig, (size_t) start, (size_t) stop, strand == "+", partial, 0, gene, product);
	#if 0
	  // DNA may be truncated
    if (type == "CDS" && ! pseudo && locus. size () % 3 != 0)
    {
      cout << "Locus tag: " << prot << endl;
      locus. print (cout);
      ERROR;
    }
  #endif

    prot2loci [prot] << move (locus);
  }
}
  
  

Annot::Annot (const string &fName)
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

    string contig, prot;
    size_t start, stop;
    double score;
    char strand = ' ';
    static Istringstream iss;
    iss. reset (f. line);
    iss >> contig >> start >> stop >> prot >> score >> strand;

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
           
	  trim (prot, '_');
	  ASSERT (! prot. empty ());
    prot2loci [prot] << Locus (f. lineNum, contig, start, stop, strand == '+', false/*partial*/, 0, string (), string ());
  }
}
  
  

void Annot::load_fasta2gff_prot (const string &fName)
{
  ASSERT (fasta2gff_prot. empty ());

  LineInput f (fName);
	Istringstream iss;
	string fasta_prot, gff_prot;
  while (f. nextLine ())
  {
  	iss. reset (f. line);
  	fasta_prot. clear ();
  	gff_prot. clear ();
  	iss >> fasta_prot >> gff_prot;
  	QC_ASSERT (! gff_prot. empty ());
  	fasta2gff_prot [fasta_prot] = gff_prot;
  }
  if (fasta2gff_prot. empty ())
  	throw runtime_error ("File " + fName + " is empty");
}



void Annot::load_fasta2gff_dna (const string &fName)
{
  map<string/*DNA GFF id*/,string/*DNA FASTA id*/> gff2fasta;  
  {
    LineInput f (fName);
  	Istringstream iss;
  	string fasta_dna, gff_dna;
    while (f. nextLine ())
    {
    	iss. reset (f. line);
    	fasta_dna. clear ();
    	gff_dna. clear ();
    	iss >> fasta_dna >> gff_dna;
    	QC_ASSERT (! gff_dna. empty ());
    	gff2fasta [gff_dna] = fasta_dna;
    }
  }
  if (gff2fasta. empty ())
  	throw runtime_error ("File " + fName + " is empty");
  
  for (auto& it : prot2loci)
  {
    Set<Locus>& loci = it. second;
    for (const Locus& locus : loci)
    {
    	string s;
    	if (! find (gff2fasta, locus. contig, s))
    	  throw runtime_error ("FASTA DNA contig " + strQuote (locus. contig) + " is not found in GFF-DNA match file " + strQuote (fName));
    	var_cast (locus). contig = move (s);
    }
  }
}



const Set<Locus>& Annot::findLoci (const string &fasta_prot) const
{
  ASSERT (! fasta_prot. empty ());

  string gff_prot (fasta_prot);
  if (! fasta2gff_prot. empty ())
  {
  	string s;
  	if (! find (fasta2gff_prot, gff_prot, s))
  	  throw runtime_error ("FASTA protein " + strQuote (fasta_prot) + " is not found in GFF-protein match file");
  	gff_prot = move (s);
  }
  ASSERT (! gff_prot. empty ());
  
  const Set<Locus>* loci = findPtr (prot2loci, gff_prot);
  if (! loci)
    throw runtime_error ("FASTA protein " + fasta_prot + (fasta_prot == gff_prot ? "" : " (converted to GFF protein " + gff_prot +")") + " is misssing in .gff-file");
  ASSERT (! loci->empty ());

  return *loci;
}



}
