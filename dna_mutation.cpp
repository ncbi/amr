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
* Release changes: see amrfinder.cpp
*
*/
   
   
#undef NDEBUG 
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;
#include "alignment.hpp"
using namespace Alignment_sp;




namespace
{  

map <string/*accession*/, Vector<AmrMutation>>  accession2mutations;
unique_ptr<OFStream> mutation_all;  
string input_name;



struct BlastnAlignment : Alignment
{
	// PD-2001
 	static constexpr const size_t flankingLen = 200;  // PAR  
  string organism;
  string refAccessionFrag;
  string product;
  string gene;
  

  BlastnAlignment (const string &line,
                   const string &organism_arg)
    : Alignment (line, false, false)
    , organism (organism_arg)
    {
	    replace (organism, '_', ' ');
      try
      {
        // refName =      NC_022347.1@23S_ribosomal_RNA@23S@-159:1040292-1037381
        //           accesion_version@gene_name@gene_symbol@offset:start-stop
	      ASSERT (! refName. empty ());
	      // PD-3419
	      {
  	      string s (refName);
  	      refAccessionFrag =       findSplit (s, '@');
  	      product          =       findSplit (s, '@');
  	      gene             =       findSplit (s, ':');
  	    //ref_offset       = stoi (findSplit (s, ':'));
  	      refAccessionFrag += ":" + s;
  	    }
	      replace (product, '_', ' ');
	      qc ();
  	    if (const Vector<AmrMutation>* refMutations = findPtr (accession2mutations, refName))
  	      setSeqChanges (*refMutations, flankingLen);
		  }
		  catch (...)
		  {
		  	cout << line << endl;
		  	throw;
		  }
    }
  void qc () const override
    { if (! qc_on)
        return;
      Alignment::qc ();
      QC_ASSERT (! refAccessionFrag. empty ());
      QC_ASSERT (! product. empty ());
      QC_ASSERT (! gene. empty ());
      QC_ASSERT (! organism. empty ());
    }
  void saveText (ostream& os) const override
    { const string na ("NA");
      for (const SeqChange& seqChange : seqChanges)
      {
        VectorPtr<AmrMutation> mutations (seqChange. mutations);
        if (mutations. empty ())
          mutations << nullptr;
        for (const AmrMutation* mutation : mutations)
        {
          ASSERT (! (seqChange. empty () && ! mutation));
          TabDel td (2, false);
    	    if (! input_name. empty ())
    	      td << input_name;;
          td << na  // PD-2534
             << nvl (targetName, na)
             << (empty () ? 0 : targetStart + 1)
             << (empty () ? 0 : targetEnd)
             << (empty () ? na : (targetStrand ? "+" : "-"))
             << (mutation
                   ? seqChange. empty ()
                     ? mutation->wildtype ()
                     : mutation->geneMutation
                   : gene + "_" + seqChange. getMutationStr ()
                )
             << (mutation
                   ? seqChange. empty ()
                       ? organism + " " + product + " [WILDTYPE]"
                       : mutation->name
                   : organism + " " + product + " [UNKNOWN]"
                )
             << "core"  // PD-2825
             // PD-1856
             << "AMR"
             << "POINT"
             << (mutation ? nvl (mutation->classS,   na) : na)
             << (mutation ? nvl (mutation->subclass, na) : na);
           if (empty ())
             td << na
                << na
                << na
                << na
                << na
                << na
                << na
                << na;
           else
             td << "POINTN"  // PD-2088
                << targetEnd - targetStart  // was: targetLen  // PD-3796
                << refLen
                << refCoverage () * 100  
                << pIdentity () * 100  
                << targetSeq. size ()
                << refAccessionFrag  // refName
                << product;  // pm.gene
          // HMM
          td << na
             << na;
          if (! seqChange. empty () && mutation && ! seqChange. replacement)  // resistant mutation
            os << td. str () << endl;
          if (mutation_all. get ())
  	        *mutation_all << td. str () << endl;
  	    }
      }
    }
    

  bool good () const
    { return targetSeq. size () >= min (refLen, 2 * flankingLen + 1); }
#if 0
  bool operator< (const BlastnAlignment &other) const
    { LESS_PART (*this, other, targetName);
      LESS_PART (other, *this, pIdentity ());
      LESS_PART (*this, other, targetStart);
      LESS_PART (*this, other, refName);
      return false;
    }
#endif
};




struct Batch
{
  VectorOwn<BlastnAlignment> blastAls;   
  
  
  explicit Batch (const string &mutation_tab)
	  {
	    {
        LineInput f (mutation_tab);
  	  	string accession, geneMutation, classS, subclass, name;
  			int pos;
   	  	Istringstream iss;
    	  while (f. nextLine ())
    	  {
	  	    if (isLeft (f. line, "#"))
	  	      continue;
     	  	iss. reset (f. line);
    	  	iss >> accession >> pos >> geneMutation >> classS >> subclass >> name;
    	  	QC_ASSERT (pos > 0);
   	  		accession2mutations [accession]. push_back (move (AmrMutation ((size_t) pos, geneMutation, classS, subclass, name)));
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
	    if (! input_name. empty ())
	      td << "Name";
	    td << "Protein identifier"   // targetName  // PD-2534
         // Contig
         << "Contig id"
         << "Start"  // targetStart
         << "Stop"  // targetEnd
         << "Strand"   // targetStrand
	       << "Gene symbol"
	       << "AmrMutation name"
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
      if (mutation_all. get ())
        *mutation_all << td. str () << endl;
	  }

  	for (const BlastnAlignment* blastAl : blastAls)
  	{
  	  ASSERT (blastAl);
   	  blastAl->qc ();
   	  blastAl->saveText (os);
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
      addPositional ("organism", "Organism name");
      addKey ("mutation_all", "File to report all mutations");
      addKey ("name", "Text to be added as the first column \"name\" to all rows of the report");
	    version = SVN_REV;
    }



  void body () const final
  {
    const string blastnFName        = getArg ("blastn");
    const string mutation_tab       = getArg ("mutation");  
    const string organism           = getArg ("organism");  
    const string mutation_all_FName = getArg ("mutation_all");
                 input_name         = getArg ("name");
    
    
    if (! mutation_all_FName. empty ())
      mutation_all. reset (new OFStream (mutation_all_FName));


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
  	    auto al = new BlastnAlignment (f. line, organism);
  	    al->qc ();  
  	    if (al->good ())
  	      batch. blastAls << al;
  	    else
  	      delete al;
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
		  for (const BlastnAlignment* blastAl : batch. blastAls)
		  {
		    ASSERT (blastAl);
		    blastAl->saveText (cout);
		    cout << ' ' << blastAl->seqChanges. size () << endl;
		  }
		}
		
    
    for (const BlastnAlignment* blastAl1 : batch. blastAls)
      for (const SeqChange& seqChange1 : blastAl1->seqChanges)
      {
        ASSERT (seqChange1. al == blastAl1);
      //ASSERT (seqChange1. mutation);
        for (const BlastnAlignment* blastAl2 : batch. blastAls)
          if (   blastAl2->targetName   == blastAl1->targetName
              && blastAl2->targetStrand == blastAl1->targetStrand
              && blastAl2 != blastAl1
             )  
          //for (Iter<Vector<SeqChange>> iter (var_cast (blastAl2) -> seqChanges); iter. next (); )
            for (SeqChange& seqChange2 : var_cast (blastAl2) -> seqChanges)
            {
            //SeqChange& seqChange2 = *iter;
              ASSERT (seqChange2. al == blastAl2);
            //ASSERT (seqChange2. mutation);
              if (   seqChange1. start_target == seqChange2. start_target 
                  && seqChange1. better (seqChange2)                
                 )
              //iter. erase ();
                seqChange2. replacement = & seqChange1;
            }
      }
		
		
  #if 0
  	// [UNKNOWN]
  	{
    	map<AmrMutation, const AmrMutation*> mutation2ptr;
    	for (const auto& it : accession2mutations)
    	  for (const AmrMutation& mut : it. second)
    	    mutation2ptr [mut] = & mut;
    	for (const BlastnAlignment* al : batch. blastAls)
    	  for (const SeqChange& seqChange : al->seqChanges)
    	    if (const AmrMutation* mut = seqChange. mutation)
    	      mutation2ptr. erase (*mut);
    	for (const auto& it : mutation2ptr)
    	{
    	  const auto al = new BlastnAlignment (* it. second);
    	  batch. blastAls << al;
    	}
    }
  #endif


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



