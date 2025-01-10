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
*   Check the correctness of a GFF file
*
*/
   
   
#undef NDEBUG 

#include "common.hpp"
using namespace Common_sp;
#include "gff.hpp"
using namespace GFF_sp;

#include "common.inc"



namespace 
{


const string locus_tagS ("[locus_tag=");
const string prodigal_ID (" ID=");
const string noFile ("emptystring");



struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Check the correctness of a GFF file. Exit with an error if it is incorrect.")
    {
      // Input
      addPositional ("gff", "GFF file, if " + strQuote (noFile) + " then exit 0");
      addKey ("gfftype", "Type of GFF file: " + Gff::names. toString (", "), "genbank");
      addKey ("prot", "Protein FASTA file");
      addKey ("dna", "DNA FASTA file");
      addFlag ("lcl", "Nucleotide FASTA created by PGAP has \"lcl|\" prefix in accessions");  
      // Output
      addKey ("gff_prot_match", "Output file with pairs: \"<protein FASTA id> <protein GFF id>\", \n\
where for genbank: <protein GFF id> is from " + strQuote (locus_tagS + "<id>") + " in the protein FASTA comment, \n\
for microscope: <protein GFF id> is ID:<num> from '><acc>|ID:<num>|<gene>|', \n\
and for prodigal: <protein GFF id> is ID=<num> in the protein FASTA comment\n\
");
      addKey ("gff_dna_match",  "Output file with pairs: \"<DNA FASTA id> <DNA GFF id>\", where for pseudomonasdb: <DNA GFF id> is the suffix after '|' in the DNA FASTA identifier");
	    version = SVN_REV;
    }



  void body () const final
  {
    const string    gffName        = getArg ("gff");
    const Gff::Type type           = Gff::name2type (getArg ("gfftype"));
    const string    protFName      = getArg ("prot");
    const string    dnaFName       = getArg ("dna");
    const string    protMatchFName = getArg ("gff_prot_match");
    const string    dnaMatchFName  = getArg ("gff_dna_match");
    const bool      lcl            = getFlag ("lcl"); 
    
    if (lcl && type != Gff::pgap)
      throw runtime_error ("-lcl requires type pgap");
    
    
    if (isRight (gffName, noFile))
    	return;
    

    const Annot annot (gffName, type, ! protMatchFName. empty (), lcl);
    
    
    if (! protFName. empty ())
    {
	    StringVector gffIds;  gffIds. reserve (10000);  // PAR
    	{
	    	OFStream outF;
	    	if (! protMatchFName. empty ())
	    		outF. open ("", protMatchFName, "");
	    	StringVector fastaIds;  fastaIds. reserve (gffIds. capacity ());
			  LineInput f (protFName /*, 100 * 1024, 1*/);
    		Istringstream iss;
    		string line_orig;
     		string fastaId;
			  while (f. nextLine ())
			  {
          trimTrailing (f. line);
			    if (f. line. empty ())
			      continue;
		    	if (f. line [0] != '>')
		    	  continue;
	    		line_orig = f. line;
	    		iss. reset (f. line. substr (1));
	    		fastaId. clear ();
    		  iss >> fastaId;
    		  QC_ASSERT (! fastaId. empty ());
	    		ASSERT (! contains (fastaId, ' '));
    			fastaIds << fastaId;
    		  // gffId
    			string gffId (fastaId);
	    		if (! protMatchFName. empty ())
	    		  switch (type)
	    		  {
	    		    case Gff::genbank:
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
      	    		break;
      	    	case Gff::microscope:
      	    	  {
      	    	    string s (std::move (gffId));
      	    	    findSplit (s, '|');
      	    	    gffId = findSplit (s, '|');
      	    	    const string idS ("ID:");
      	    	    if (! isLeft (gffId, idS))
      	    	      throw runtime_error (__FILE__ ": 'ID:' is not found in: " + line_orig);
      	    	    gffId. erase (0, idS. size ());
      	    	  }      	    	  
      	    	  break;
	    		    case Gff::prodigal:
      	    		{
      	    			const size_t pos = f. line. find (prodigal_ID);
      	    			if (pos == string::npos)
      	    				throw runtime_error (__FILE__ ": " + strQuote (prodigal_ID) + " is not found in: " + line_orig);
      	    			gffId = f. line. substr (pos + prodigal_ID. size ());
      	    			const size_t end = gffId. find (';');
      	    			if (end == string::npos)
      	    				throw runtime_error (__FILE__ ": ';' is not found after " + strQuote (prodigal_ID) + " in: " + line_orig);
      	    		  gffId. erase (end);
      	    		}
      	    		break;
      	      default: break;
      	    }
	    	  //
	    		if (contains (gffId, ' '))
	    			throw runtime_error (__FILE__ ": " + strQuote (gffId) + " contains space");
	    		if (gffId. empty ())
	    			throw runtime_error (__FILE__ ": No protein identifier in: " + line_orig);
    			gffIds << gffId;
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
			  ASSERT (gffIds. size () == fastaIds. size ());
			}
			if (verbose ())
			  cout << "# Proteins in GFF: " << annot. prot2loci. size () << endl;
		  for (const string& seqid : gffIds)
		  	if (! contains (annot. prot2loci, seqid))
		  		throw runtime_error (__FILE__ ": Protein FASTA id " + strQuote (seqid) + " is not in the GFF file");
    #if 0
		  for (const auto& it : annot. prot2loci)
		    if (! gffIds. containsFast (it. first))
		  		throw runtime_error (__FILE__ ": GFF protein id " + strQuote (it. first) + " is not in the protein FASTA file");  // pseudogene ??
		#endif
    }   


    if (! dnaFName. empty ())
    {
    	StringVector contigIds;  contigIds. reserve (10000);  // PAR
    	StringVector gffIds;  gffIds. reserve (10000);  // PAR
    	{
      	OFStream outF;
      	if (! dnaMatchFName. empty ())
      		outF. open ("", dnaMatchFName, "");
			  LineInput f (dnaFName /*, 100 * 1024, 1*/);
    		Istringstream iss;
    		string contigId;
			  while (f. nextLine ())
	    	{
          trimTrailing (f. line);
			    if (f. line. empty ())
			      continue;
		    	if (f. line [0] != '>')
		    	  continue;
	    		iss. reset (f. line. substr (1));
	    		contigId. clear ();
    		  iss >> contigId;
	    		ASSERT (! contains (contigId, ' '));
    		  // gffId
    			string gffId (contigId);
	    		if (! dnaMatchFName. empty ())
	    		  switch (type)
	    		  {
	    		    case Gff::pseudomonasdb:
      	    		{
      	    			const size_t pos = gffId. rfind ('|');
      	    			if (pos != string::npos)
      	    				gffId. erase (0, pos + 1);
      	    		}
      	    		break;
      	      default: break;
      	    }
      	  //
	    		if (gffId. empty ())
	    			throw runtime_error (__FILE__ ": No contig identifier in:\n" + f. line);
	    		if (lcl && ! isLeft (gffId, "lcl|"))
	    			throw runtime_error (__FILE__ ": Contig identifier does not start with " + strQuote ("lcl|") + ":\n" + f. line);
    			gffIds << gffId;
    			contigIds << contigId;
    			if (outF. is_open ())
    				outF << contigId << '\t' << gffId << endl;
	    	}
			}
			ASSERT (contigIds. size () == gffIds. size ());
		  gffIds. sort ();
		  {
			  const string* s_prev = nullptr;
			  for (const string& s : gffIds)
			  {
			    if (s_prev && *s_prev == s)
  			  	throw runtime_error (__FILE__ ": DNA GFF identifier " + strQuote (s) + " is not unique");
			    s_prev = & s;
			  }
			}
		  contigIds. sort ();
		  {
			  const string* s_prev = nullptr;
			  for (const string& s : contigIds)
			  {
			    if (s_prev && *s_prev == s)
  			  	throw runtime_error (__FILE__ ": DNA contig identifier " + strQuote (s) + " is not unique");
			    s_prev = & s;
			  }
			}
			for (const auto& it : annot. prot2loci)
			  for (const Locus& cds : it. second)
			    if (! gffIds. contains (cds. contig))
  		  		throw runtime_error (__FILE__ ": GFF contig id " + strQuote (cds. contig) + " is not in the DNA FASTA file");
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



