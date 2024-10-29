// mutate.cpp

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
*   Mutate a FASTA file
*
*/


#undef NDEBUG

#include "common.hpp"
using namespace Common_sp;
#include "alignment.hpp"
using namespace Alignment_sp;
#include "seq.hpp"
using namespace Seq_sp;

#include "common.inc"



namespace 
{



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Mutate a FASTA file")
    {
  	  addPositional ("in", "Input FASTA file");
  	  addPositional ("mut", "AmrMutation table: <seq_id> <1-based pos> <mutation_std> <mutation_report>");
  	  addFlag ("aa", "Protein/DNA");
  	  addFlag ("orig", "Add the original, non-mutated sequences");

 	    version = SVN_REV;  
    }


	
	void body () const final
  {
	  const string inFName  = getArg ("in");
	  const string mutFName = getArg ("mut");
	  const bool aa         = getFlag ("aa");
	  const bool orig       = getFlag ("orig");
  
  
    map <string/*seqId*/,Vector<AmrMutation>> id2mutation;
    {
      LineInput in (mutFName);
      Istringstream iss;
      while (in. nextLine ())
      {
        iss. reset (in. line);
        string seqId;
        size_t pos;
        string mutation_std;
        string mutation_report;
        iss >> seqId >> pos >> mutation_std >> mutation_report;
        QC_ASSERT (! mutation_report. empty ());
        AmrMutation mut (pos, mutation_std, mutation_report, "X", "X", "X");
        mut. qc ();
        id2mutation [seqId] << std::move (mut);
      }
    }
  
  
	  Multifasta fIn (inFName, aa);
	  while (fIn. next ())
	  {
 	    unique_ptr<const Seq> seq;
	    try
  	  {
  	    if (aa)
  	    {
  	      auto pep = new Peptide (fIn, 1000, false);  // PAR
  	      pep->pseudo = true;
  	      seq. reset (pep);  
  	    }
  	    else
  	      seq. reset (new Dna     (fIn, 100000, false));  // PAR
  	    seq->qc ();
  	    if (orig)
  	      seq->saveText (cout);
  	    if (const Vector<AmrMutation>* muts = findPtr (id2mutation, seq->getId ()))
  	      for (const AmrMutation& mut : *muts)
    	    {
    	      unique_ptr<Seq> seq1 (seq->copy ());
    	      mut. apply (seq1->seq);
    	      if (! aa)
    	        strLower (seq1->seq);
      	    seq1->name += ":" + to_string (mut. pos_real + 1) + ":" + mut. geneMutation;
      	    seq1->qc ();
      	    seq1->saveText (cout);
    	    }
  	  }
  	  catch (const exception &e)
  	  {
  	    if (seq)
  	      throw runtime_error (seq->name + "\n" + e. what ());
  	    throw;
  	  }
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



