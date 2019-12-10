// dna_mutation.cpp

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
*   Identification of mutations at DNA level
*
* Release changes:
*   3.2.4 11/15/2019 PD-3191  neighborhoodMismatch <= 0.04; good(): length >= min (refLen, 2 * flankingLen + 1)
*
*/
   
   
#undef NDEBUG 
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;
#include "alignment.hpp"
using namespace Alignment_sp;




// PD-3096
// PAR!
#ifdef SVN_REV
  #define SOFTWARE_VER SVN_REV
#else
  #define SOFTWARE_VER "3.4.2"
#endif



namespace
{  

map <string/*accession*/, Vector<Mutation>>  accession2mutations;
unique_ptr<OFStream> mutation_all;  // ??



struct BlastnAlignment : Alignment
// BLASTN alignment
{
	// PD-2001
 	static constexpr const size_t flankingLen = 200;  // PAR  
  string refAccessionFrag;
  string product;
  

  explicit BlastnAlignment (const string &line)
    : Alignment (line, false, false)
    {
      try
      {
	      ASSERT (! refName. empty ());
	      string s (refName);
	      refAccessionFrag = findSplit (s, '@');
	      product = findSplit (s, ':');
	      replace (product, '_', ' ');
	      QC_ASSERT (! s. empty ());
	      refAccessionFrag += ":" + s;

  	    if (const Vector<Mutation>* refMutations = findPtr (accession2mutations, refName))
  	      setSeqChanges (*refMutations, flankingLen, mutation_all. get ());
		  }
		  catch (...)
		  {
		  	cout << line << endl;
		  	throw;
		  }
    }
  void saveText (ostream& os) const 
    { const string na ("NA");
      for (const SeqChange& seqChange : seqChanges)
      {
        if (! seqChange. mutation)
          continue;
        TabDel td (2, false);
        td << na  // PD-2534
           << targetName 
           << targetStart + 1
           << targetEnd
           << (targetStrand ? '+' : '-');
        td << seqChange. mutation->geneMutation
           << seqChange. mutation->name + ifS (seqChange. empty (), " [NO_CALL]")
           << "core"  // PD-2825
           // PD-1856
           << "AMR"
           << "POINT"
           << nvl (seqChange. mutation->classS, na)
           << nvl (seqChange. mutation->subclass, na)
           //
           << "POINTN"  // PD-2088
           << targetLen;
        td << refLen
           << refCoverage () * 100  
           << pIdentity () * 100  
           << targetSeq. size ()
           << refAccessionFrag  // refName
           << product  // pm. gene
           ;
        // HMM
        td << na
           << na;
        os << td. str () << endl;
      }
    }
    

  bool good () const
    { return targetSeq. size () >= min (refLen, 2 * flankingLen + 1); }
  bool operator< (const BlastnAlignment &other) const
    { LESS_PART (*this, other, targetName);
      LESS_PART (other, *this, pIdentity ());
      LESS_PART (*this, other, targetStart);
      LESS_PART (*this, other, refName);
      return false;
    }
};




struct Batch
{
  vector<BlastnAlignment> blastAls;   
  
  
  explicit Batch (const string &mutation_tab)
	  {
	    {
        LineInput f (mutation_tab);
  	  	string accession, geneMutation, classS, subclass, name;
  			int pos;
   	  	Istringstream iss;
    	  while (f. nextLine ())
    	  {
     	  	iss. reset (f. line);
    	  	iss >> accession >> pos >> geneMutation >> classS >> subclass >> name;
    	  	QC_ASSERT (pos > 0);
   	  		accession2mutations [accession]. push_back (move (Mutation ((size_t) pos, geneMutation, classS, subclass, name)));
    	  }	    
    	}
  	  for (auto& it : accession2mutations)
  	  {
  	  	it. second. sort ();
  	    if (! it. second. isUniq ())
  	  	  throw runtime_error ("Duplicate reference mutations for " + it. first);
  	  }
	  }
	  	  

	void report (ostream &os) const
	{
    {
    	// Cf. BlastnAlignment::saveText()
	    TabDel td;
	    td << "Protein identifier"   // targetName  // PD-2534
         // Contig
         << "Contig id"
         << "Start"  // targetStart
         << "Stop"  // targetEnd
         << "Strand"   // targetStrand
	       << "Gene symbol"
	       << "Mutation name"
	       << "Scope"  // PD-2825
	       // PD-1856
	       << "Element type"
	       << "Element subtype"
	       << "class"
	       << "Subclass"
	       //
	       << "Method"
	       << "Target length" 
	       //
	       << "Reference gene length"         // refLen
	       << "% Coverage of reference gene"  // queryCoverage
	       << "% Identity to reference gene"  
	       << "Alignment length"                 // length
	       << "Accession of reference gene"    
	       << "Name of reference gene"
	       //
	       << "HMM id"
	       << "HMM description"
	       ;
	    os << td. str () << endl;
	  }

  	for (const auto& blastAl : blastAls)
  	{
   	  blastAl. qc ();
   	  blastAl. saveText (os);
    }
	}
};




// ThisApplication

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Find mutations at DNA level and report in the format of amr_report.cpp")
    {
      addPositional ("blastn", string ("blastn output in the format: ") + Alignment::format + ". sseqid is the 1st column of <mutation_tab> table");  
      addPositional ("mutation", "Mutations table");
	    version = SOFTWARE_VER;
    }



  void body () const final
  {
    const string blastnFName  = getArg ("blastn");
    const string mutation_tab = getArg ("mutation");  
    
    

    Batch batch (mutation_tab);  
  
  
    // Input 
    {
      LineInput f (blastnFName);
  	  while (f. nextLine ())
  	  {
  	    { 
  	      Unverbose unv;
  	      if (verbose ())
  	        cout << f. line << endl;  
  	    }
  	    BlastnAlignment al (f. line);
  	    al. qc ();  
  	    if (al. good ())
  	      batch. blastAls. push_back (move (al));
  	  }
  	}
  	if (verbose ())
  	  cout << "# Good Blasts: " << batch. blastAls. size () << endl;
  	
    
    // Output
    // Group by targetName and process each targetName separately for speed ??    
  //Common_sp::sort (batch. blastAls);
    if (verbose ())
    {
	    cout << "After process():" << endl;
		  for (const auto& blastAl : batch. blastAls)
		  {
		    blastAl. saveText (cout);
		    cout << ' ' << blastAl. seqChanges. size () << endl;
		  }
		}
		
    
    for (const BlastnAlignment& blastAl1 : batch. blastAls)
      for (const SeqChange& seqChange1 : blastAl1. seqChanges)
      {
        if (! seqChange1. mutation)
          continue;
        for (BlastnAlignment& blastAl2 : batch. blastAls)
          if (   blastAl2. targetName   == blastAl1. targetName
              && blastAl2. targetStrand == blastAl1. targetStrand
              && & blastAl2 != & blastAl1
             )  
            for (SeqChange& seqChange2 : blastAl2. seqChanges)
            {
              if (! seqChange2. mutation)
                continue;
              if (   seqChange1. start_target         == seqChange2. start_target 
                  && seqChange1. neighborhoodMismatch <  seqChange2. neighborhoodMismatch
                 )
                seqChange2 = SeqChange ();
            }
      }
		
		
    batch. report (cout);
  }
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



