// alignment.hpp

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
*   Protein or DNA mutations library.
*
*/
   
   
#include "common.hpp"
using namespace Common_sp;
#include "seq.hpp"
using namespace Seq_sp;


   

namespace Alignment_sp
{



static constexpr char pm_delimiter = '_';



struct AmrMutation final : Root
// Database
{
	size_t pos_real {0};  
	  // In whole reference sequence
	  // = start of reference
	
	string geneMutation_std;
  // Function of geneMutation_std
  // Upper-case
  string reference;
	string allele;
	string gene;
	int pos_std {0};  
	size_t frameshift {no_index};
	  // Position of '*' after getStop()
	int frameshift_insertion {0};

	// To be reported
	// !empty()
	string geneMutation;
	string classS;
	string subclass;
	string name;
	  // Species binomial + resistance
	
	
  // Input: pos_arg: 1-based
	AmrMutation (size_t pos_real_arg,
  		 				 const string &geneMutation_std_arg,
  		 				 const string &geneMutation_arg,
  			 			 const string &class_arg,    
  				 		 const string &subclass_arg, 
  					 	 const string &name_arg);    
	AmrMutation (size_t pos_arg,
  		 				 const string &geneMutation_std_arg)
    : AmrMutation (pos_arg, geneMutation_std_arg, geneMutation_std_arg, "X", "X", "X")
    {}
	AmrMutation () = default;
private:
	static void parse (const string &geneMutation_std,
	                   string &reference,
	                   string &allele,
                     string &gene,
                     int &pos_std,
                     size_t &frameshift,
                     int &frameshift_insertion);
public:
  void qc () const override;
  void saveText (ostream &os) const override
    { if (empty ())
        os << "empty";
      else
        os        << pos_real + 1 
           << ' ' << geneMutation 
           << ' ' << frameshift_insertion 
           << ' ' << name; 
    }
  bool empty () const override
    { return geneMutation_std. empty (); }


  size_t getStop () const
    { return pos_real + reference. size (); }
  string wildtype () const
    { return gene + "_" + reference + to_string (pos_std + 1) + reference; }
  bool operator< (const AmrMutation &other) const;
  bool operator== (const AmrMutation &other) const
    { return geneMutation_std == other. geneMutation_std; }
  void apply (string &seq) const
    { if (pos_real >= seq. size ())
        throw runtime_error ("AmrMutation position " + to_string (pos_real) + " is outside the sequence: " + seq);
      if (frameshift != no_index)
        throw runtime_error ("AmrMutation is a frameshift");
      if (verbose ())
        cerr         << seq. substr (0, pos_real) 
             << endl << allele 
             << endl << seq. substr (pos_real + reference. size ())
             << endl;
      seq = seq. substr (0, pos_real) + allele + seq. substr (pos_real + reference. size ());
    }
};



struct Alignment;



struct SeqChange final : Root
// Observation
{
  const Alignment* al {nullptr};
    // !nullptr    
//bool fromAllele {false};
  
  // In alignment
  size_t start {0};
  size_t len {0};
  
  // No '-'
  string reference;
    // Insertion => start is artifically decremented and len is incremented => !empty()
  string allele;
  // empty() <=> frame shift

  size_t start_ref {0};
  size_t stop_ref {0};
  size_t start_target {0};
	double neighborhoodMismatch {0.0};
	  // 0..1
	  
	VectorPtr<AmrMutation> mutations;
	  // !nullptr
	  // Matching AmrMutation's
    		  
	Disruption disr;

	const SeqChange* replacement {nullptr}; 
	  // !nullptr => *this is replaced by *replacement
  
  
  SeqChange () = default;
  explicit SeqChange (const Alignment* al_arg/*,
             bool fromAllele_arg*/)
    : al (al_arg)
  //, fromAllele (fromAllele_arg)
    {}
  SeqChange (const Alignment* al_arg,
             const AmrMutation* mutation_arg)
    : al (al_arg)
    { mutations << checkPtr (mutation_arg); }
  SeqChange (const Alignment* al_arg,
             const Disruption &disr_arg)
    : al (al_arg)
    , disr (disr_arg)
    {}
  void qc () const override;
  void saveText (ostream &os) const override
    { os        << start + 1 
         << ' ' << len 
         << ' ' << strQuote (reference) << " -> " << strQuote (allele)
         << ' ' << start_ref + 1 << ".." << stop_ref
         << ' ' << start_target + 1 
         << ' ' << neighborhoodMismatch;
      if (! disr. empty ())
        disr. saveText (os);
      for (const AmrMutation* mutation : mutations)
      { os << ' ' ;
        mutation->saveText (os);
      }
      os << endl; 
    }
  bool empty () const override
    { return ! len && disr. empty (); }
    
    
  bool hasMutation () const
    { return ! empty () && ! mutations. empty () && ! replacement; }
  bool hasFrameshift () const
    { return hasMutation () && mutations [0] -> frameshift != no_index; }
  bool isFrameshift () const
    { return reference. empty (); }
  string getMutationStr () const;
  size_t getStop () const
    { return start + len; }
  bool operator< (const SeqChange &other) const;
  bool better (const SeqChange &other) const;
  bool finish (const string &refSeq,
               size_t flankingLen);
    // Return: good match
    // Invokes: finishPos()
  bool finishPos (size_t flankingLen);
    // Return: good match
private:
  void setSeq ();
  void setStartStopRef ();
  void setStartTarget ();
  void setNeighborhoodMismatch (size_t flankingLen);
public:
  bool matchesMutation (const AmrMutation& mut) const;
};



struct Alignment : Hsp
{
  AmrMutation refMutation;
    // !empty() => qseq contains AmrMutation::allele
//int ref_offset {0};

  Vector<SeqChange> seqChanges;

  
  Alignment (const string &line,
             bool qProt_arg,  
             bool sProt_arg)
    : Hsp (line, qProt_arg, sProt_arg, qProt_arg || sProt_arg /*aProt*/, /*false*/ qProt_arg /*qStopCodon*/, true/*bacterialStartCodon*/)
    {}
  Alignment () = default;  
protected:
  void setSeqChanges (const Vector<AmrMutation> &refMutations,
                      size_t flankingLen/*,
                      bool allMutationsP*/);
    // Input: flankingLen: valid if > 0
private:
  size_t refMutation2refSeq_pos ();
    // Return: no_index <=> refMutation is not detected
public:
  void qc () const override;
  void saveText (ostream &os) const override
    { Hsp::saveText (os);
      if (! refMutation. empty ())
        os << ' ' << refMutation;
	    os << " #seqChanges:" << seqChanges. size ();
    }


  bool hasMutation () const
    { for (const SeqChange& seqChange : seqChanges)
        if (seqChange. hasMutation ())
          return true;
      return false;
    }
  bool hasDeclarativeFrameshift () const
    { return seqChanges. size () == 1 && seqChanges [0]. hasFrameshift (); }
};




}  // namespace


