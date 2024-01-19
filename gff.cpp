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

#include "gff.hpp"

#include "common.inc"



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

const StringVector Gff::names {"bakta", "genbank", "microscope", "patric", "pgap", "prodigal", "prokka", "pseudomonasdb", "rast", "standard"};


Gff::Type Gff::name2type (const string &name)
{
  if (name == "bakta")          return bakta;
  if (name == "genbank")        return genbank;
  if (name == "microscope")     return microscope;
  if (name == "patric")         return patric;
  if (name == "pgap")           return pgap;
  if (name == "prodigal")       return prodigal;
  if (name == "prokka")         return prokka;
  if (name == "pseudomonasdb")  return pseudomonasdb;
  if (name == "rast")           return rast;
  if (name == "standard")       return standard;
  throw runtime_error ("Unknown GFF type: " + strQuote (name));
}




// Annot

namespace
{
  
string unescape (const string &s)
{
  string r (unpercent (s));
  trim (r);
  return r;
}

  
  
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

  QC_ASSERT (! accession. empty ());
}
    
}



Annot::Annot (const string &fName,
              Gff::Type gffType,
	            bool protMatch,
	            bool lcl)
{
  IMPLY (protMatch,    gffType == Gff::genbank 
                    || gffType == Gff::microscope 
                    || gffType == Gff::prodigal
        );
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

    try
    {
      /*1*/       string contig     (unescape (findSplit (f. line, '\t')));
      /*2*/ const string source     (unescape (findSplit (f. line, '\t')));
      /*3*/ const string type       (unescape (findSplit (f. line, '\t')));
      /*4*/ const string startS     (unescape (findSplit (f. line, '\t')));
      /*5*/ const string stopS      (unescape (findSplit (f. line, '\t')));
      /*6*/ const string score      (unescape (findSplit (f. line, '\t')));  // real number
      /*7*/ const string strand     (unescape (findSplit (f. line, '\t')));
      /*8*/ const string phase      (unescape (findSplit (f. line, '\t')));  // frame
      /*9*/ string attributes (f. line);  
      
      trim (attributes);
      if (attributes. empty ())
      	throw runtime_error ("9 fields are expected in each line");

    #if 0
      if (trimProject)
        if (contains (contig, ":"))
          findSplit (contig, ':');  // = project_id
    #endif
      if (contig. empty ())
      	throw runtime_error ("empty sequence indentifier");
  	  for (const char c : contig)
  	  	if (! printable (c))
  	  		throw runtime_error ("Non-printable character in the sequence identifier: " + to_string (c));

      if (   type != "CDS"
          && type != "gene"
          && type != "pseudogene"
         )
        continue;
        
      if (gffType == Gff::pgap && type != "CDS")
        continue;
        
      long start = -1;
      if (! str2<long> (startS, start))
	    	throw runtime_error ("Cannot read start");
      if (start <= 0)
      	throw runtime_error ("start should be >= 1");
        
      long stop = -1;
      if (! str2<long> (stopS, stop))
 	    	throw runtime_error ("Cannot read stop");
      if (stop <= 0)
      	throw runtime_error ("stop should be >= 1");

      if (start > stop)
      	throw runtime_error ("start cannot be greater than stop");

      start--;    	
      	
      if (   strand != "+" 
          && strand != "-"
         )
      	throw runtime_error ("strand should be '+' or '-'");
             
      const bool pseudo =    contains (attributes, "pseudo=true")
                          || contains (attributes, "gene_biotype=pseudogene")
                          || type == "pseudogene";
    //if (pseudo && type == "CDS")  
      //continue;
      
      const bool partial =    contains (attributes, "partial=true")
                              // Gff::prodigal
                           || contains (attributes, "partial=01")
                           || contains (attributes, "partial=10")
                           || contains (attributes, "partial=11");
      
      string protAttr = "Name";
      switch (gffType)
      {
        case Gff::bakta:         protAttr = "ID"; break;
        case Gff::genbank:       protAttr = (protMatch || pseudo) ? "locus_tag" : "Name"; break;
        case Gff::microscope:    protAttr = "ID"; break;
        case Gff::patric:        protAttr = "ID"; break;
        case Gff::prodigal:      protAttr = "ID"; break;
        case Gff::prokka:        protAttr = "ID"; break;
        case Gff::pseudomonasdb: protAttr = "Alias"; break;  // for type = "gene", "locus" for type = "CDS"
        case Gff::rast:          protAttr = "ID"; break;
        default: break;
      }
      ASSERT (! protAttr. empty ());
      protAttr += "=";
          
      string prot_;
      string gene_;
      string product_;
      {
        string locusTag;
        while (! attributes. empty ())
        {
    	    string s (findSplit (attributes, ';'));
      	  trim (s);
    	    if (isLeft (s, protAttr))
    	    {
    	      prot_ = s;
            findSplit (prot_, '='); 
    	    }
    	    else if (isLeft (s, "gene="))
    	    {
    	      gene_ = s;
    	      findSplit (gene_, '=');
    	    }
    	    else if (isLeft (s, "product="))
    	    {
    	      product_ = s;
    	      findSplit (product_, '=');
    	    //replace (product, tmpSpace, ' ');
    	    }
    	    else if (gffType == Gff::patric && isLeft (s, "locus_tag="))
    	    {
    	      locusTag = s;
    	      findSplit (locusTag, '=');
    	    }
    	  }
    	  trimPrefix (prot_, "\"");
    	  trimSuffix (prot_, "\"");
        if (prot_. empty ())
        	continue;
        //throw runtime_error ("no attribute '" + protAttr + "': " + f. line);
        
        switch (gffType)
        {
          case Gff::genbank:       if (contains (prot_, ":"))  findSplit (prot_, ':');  break;
          case Gff::patric:        if (! locusTag. empty ())  prot_ += "|" + locusTag;  
                                   if (isLeft (contig, "accn|"))  contig. erase (0, 5);
                                   break;
          default:  break;
        }
      }
  	  QC_ASSERT (! prot_. empty ());
  	
  	        string prot    (unescape (prot_));
  	  const string gene    (unescape (gene_));
  	  const string product (unescape (product_));
  	
  	  if (gffType == Gff::pgap)
  	  {
  	    pgap_accession (prot, false);
  	    pgap_accession (contig, lcl);
  	  }
  	  QC_ASSERT (! prot. empty ());
  	  
  	  Locus locus ((size_t) f. lineNum, contig, (size_t) start, (size_t) stop, strand == "+", partial, 0, gene, product);
  	#if 0
  	  // DNA may be truncated
      if (type == "CDS" && ! pseudo && locus. size () % 3 != 0)
      {
        cout << "Locus tag: " << prot << endl;
        locus. print (cout);
        ERROR;
      }
    #endif

      prot2loci [prot] << std::move (locus);
    }
    catch (const exception &e)
    {
      throw runtime_error ("File " + fName + ", " + f. lineStr () + ": " + e. what ());
    }
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

   	const string errorS ("File " + fName + ", " + f. lineStr () + ": ");

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
	  		throw runtime_error (errorS + "Non-printable character in the sequence identifier: " + to_string (c));

    if (start >= stop)
    	throw runtime_error (errorS + "start should be less than stop");

    if (   strand != '+'
        && strand != '-'
       )
    	throw runtime_error (errorS + "strand should be '+' or '-'");
           
	  trim (prot, '_');
	  ASSERT (! prot. empty ());
    prot2loci [prot] << Locus ((size_t) f. lineNum, contig, start, stop, strand == '+', false/*partial*/, 0, noString, noString);
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
    	var_cast (locus). contig = std::move (s);
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
  	gff_prot = std::move (s);
  }
  ASSERT (! gff_prot. empty ());
  
  const Set<Locus>* loci = findPtr (prot2loci, gff_prot);
  if (! loci)
    throw runtime_error ("FASTA protein " + fasta_prot + (fasta_prot == gff_prot ? "" : " (converted to GFF protein " + gff_prot +")") + " is misssing in .gff-file");
  ASSERT (! loci->empty ());

  return *loci;
}



}
