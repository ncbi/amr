// gff_check.cpp

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
*   Check the correctness of a .gff-file
*
*/
   
   
#undef NDEBUG 
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;
#include "gff.hpp"
using namespace GFF_sp;



namespace 
{


const string locus_tagS ("[locus_tag=");



const string noFile ("emptystring");



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Check the correctness of a .gff-file. Exit with an error if it is incorrect.")
    {
      // Input
      addPositional ("gff", ".gff-file, if " + strQuote (noFile) + " then exit 0");
      addKey ("prot", "Protein FASTA file");
      addKey ("dna", "DNA FASTA file");
      addFlag ("gpipe", "Protein identifiers in the protein FASTA file have format 'gnl|<project>|<accession>'");
      // Output
      addKey ("locus_tag", "Output file with matches: \"<FASTA id> <GFF id>\", where <id> is from " + strQuote (locus_tagS + "<id>") + " in the FASTA comment and from the .gff-file");
    }



  void body () const final
  {
    const string gffName        = getArg ("gff");
    const string protFName      = getArg ("prot");
    const string dnaFName       = getArg ("dna");
    const string locus_tagFName = getArg ("locus_tag");
    const bool   gpipe          = getFlag ("gpipe");
    
    
    if (isRight (gffName, noFile))
    	return;
    

    Annot::Gff gff;
    const Annot annot (gff, gffName, ! locus_tagFName. empty ());
    
    
    if (! protFName. empty ())
    {
	    StringVector gffIds;  gffIds. reserve (10000);  // PAR
    	{
	    	OFStream outF;
	    	if (! locus_tagFName. empty ())
	    		outF. open ("", locus_tagFName, "");
	    	StringVector fastaIds;  fastaIds. reserve (gffIds. capacity ());
			  LineInput f (protFName /*, 100 * 1024, 1*/);
    		Istringstream iss;
    		string line_orig;
			  while (f. nextLine ())
			    if (! f. line. empty ())
			    	if (f. line [0] == '>')
			    	{
			    		line_orig = f. line;
			    		string fastaId, gffId;
			    		iss. reset (f. line. substr (1));
		    		  iss >> fastaId;
		    		  // gffId
		    			gffId = fastaId;
			    		if (! locus_tagFName. empty ())
			    		{
			    		  if (gpipe)
			    		  {
			    		    if (! trimPrefix (gffId, "gnl|"))
			    		      throw runtime_error (__FILE__ ": Protein name should start with \"gnl|\" in: " + line_orig);
  			    			const size_t pos = gffId. find ("|");
  			    			if (pos == string::npos)
  			    				throw runtime_error (__FILE__ ": " + strQuote ("|") + " is not found in: " + line_orig);
  			    			gffId. erase (0, pos + 1);
			    		  }
			    		  else
			    		  {
  			    			const size_t pos = f. line. find (locus_tagS);
  			    			if (pos == string::npos)
  			    				throw runtime_error (__FILE__ ": " + strQuote (locus_tagS) + " is not found in: " + line_orig);
  			    			gffId = f. line. substr (pos + locus_tagS. size ());
  			    			const size_t end = gffId. find (']');
  			    			if (end == string::npos)
  			    				throw runtime_error (__FILE__ ": ']' is not found after " + strQuote (locus_tagS) + " in: " + line_orig);
  			    		  gffId. erase (end);
  			    		}
			    		}
			    		ASSERT (! contains (fastaId, ' '));
			    		if (contains (gffId, ' '))
			    			throw runtime_error (__FILE__ ": " + strQuote (gffId) + " contains space");
			    		if (gffId. empty ())
			    			throw runtime_error (__FILE__ ": No protein identifier in: " + line_orig);
		    			gffIds << gffId;
		    			fastaIds << fastaId;
		    			if (outF. is_open ())
		    				outF << fastaId << '\t' << gffId << endl;
			    	}
			  const size_t n = fastaIds. size ();
			  fastaIds. sort ();
			  fastaIds. uniq ();
			  if (fastaIds. size () != n)
			  	throw runtime_error (__FILE__ ": Duplicate FASTA ids");
			  gffIds. sort ();
			  {
  			  const string* s_prev = nullptr;
  			  for (const string& s : gffIds)
  			  {
  			    if (s_prev && *s_prev == s)
    			  	throw runtime_error (__FILE__ ": GFF identifier " + strQuote (s) + " is not unique");
  			    s_prev = & s;
  			  }
  			}
			//gffIds. uniq ();
			  ASSERT (gffIds. size () == fastaIds. size ());
			}
			if (verbose ())
			  cout << "# Proteins in GFF: " << annot. prot2cdss. size () << endl;
		  for (const string& seqid : gffIds)
		  	if (! contains (annot. prot2cdss, seqid))
		  		throw runtime_error (__FILE__ ": Protein id " + strQuote (seqid) + " is not in the .gff-file");
    }   


    if (! dnaFName. empty ())
    {
    	StringVector contigIds;  contigIds. reserve (10000);  // PAR
    	{
			  LineInput f (dnaFName /*, 100 * 1024, 1*/);
    		Istringstream iss;
			  while (f. nextLine ())
			    if (! f. line. empty ())
			    	if (f. line [0] == '>')
			    	{
			    		string contigId;
			    		iss. reset (f. line. substr (1));
		    		  iss >> contigId;
			    		ASSERT (! contains (contigId, ' '));
			    		if (contigId. empty ())
			    			throw runtime_error (__FILE__ ": No contig identifier in: " + f. line);
		    			contigIds << move (contigId);
			    	}
			  contigIds. sort ();
			  {
  			  const string* s_prev = nullptr;
  			  for (const string& s : contigIds)
  			  {
  			    if (s_prev && *s_prev == s)
    			  	throw runtime_error (__FILE__ ": Contig identifier " + strQuote (s) + " is not unique");
  			    s_prev = & s;
  			  }
  			}
			//contigIds. uniq ();
			}
			for (const auto& it : annot. prot2cdss)
			  for (const Locus& cds : it. second)
			    if (! contigIds. contains (cds. contig))
  		  		throw runtime_error (__FILE__ ": Contig id " + strQuote (cds. contig) + " is not in the file " + dnaFName);
    }   
  }
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



