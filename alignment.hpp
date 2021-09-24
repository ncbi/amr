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


   

namespace Alignment_sp
{



static constexpr char pm_delimiter = '_';



struct AmrMutation : Root
// Database
{
	size_t pos {0};
	  // In whole reference sequence
	  // = start of reference
	// !empty()
	string geneMutation;
	string classS;
	string subclass;
	string name;
	  // Species binomial + resistance
	
  // Function of geneMutation
  // Upper-case
  string reference;
	string allele;
	string gene;
	int ref_pos {0};

	
	AmrMutation (size_t pos_arg,
						const string &geneMutation_arg,
						const string &class_arg = "X",
						const string &subclass_arg = "X",
						const string &name_arg = "X");
		// Input: pos_arg: 1-based
	AmrMutation () = default;
	AmrMutation (AmrMutation &&other) = default;
	AmrMutation& operator= (AmrMutation &&other) = default;
private:
	static void parse (const string &geneMutation,
	                   string &reference,
	                   string &allele,
                     string &gene,
                     int &ref_pos);
public:
  void saveText (ostream &os) const override
    { if (empty ())
        os << "empty";
      else
        os << pos + 1 << ' ' << geneMutation << ' ' << name; 
    }
  bool empty () const override
    { return geneMutation. empty (); }


  size_t getStop () const
    { return pos + reference. size (); }
  string wildtype () const
    { return gene + "_" + reference + to_string (ref_pos + 1) + reference; }
  bool operator< (const AmrMutation &other) const;
  bool operator== (const AmrMutation &other) const
    { return geneMutation == other. geneMutation; }
  void apply (string &seq) const
    { if (pos >= seq. size ())
        throw runtime_error ("AmrMutation position " + to_string (pos) + " is outside the sequence: " + seq);
      if (verbose ())
        cerr         << seq. substr (0, pos) 
             << endl << allele 
             << endl << seq. substr (pos + reference. size ())
             << endl;
      seq = seq. substr (0, pos) + allele + seq. substr (pos + reference. size ());
    }
};



struct Alignment;



struct SeqChange : Root
// Observation
{
  const Alignment* al {nullptr};
    // !nullptr
    
  
  // In alignment
  size_t start {0};
  size_t len {0};
  
  // No '-'
  string reference;
    // Insertion => !empty() by artifically decrementing start and incrementing len
  string allele;

  size_t start_ref {0};
  size_t stop_ref {0};
  size_t start_target {0};
	double neighborhoodMismatch {0.0};
	  // 0..1
	  
	char prev {'\0'};
	VectorPtr<AmrMutation> mutations;
	  // !nullptr
	
	const SeqChange* replacement {nullptr};
	  // !nullptr => *this is replaced by *replacement
  
  
  SeqChange () = default;
  explicit SeqChange (const Alignment* al_arg)
    : al (al_arg)
    {}
  SeqChange (const Alignment* al_arg,
             const AmrMutation* mutation_arg)
    : al (al_arg)
  //, mutation (checkPtr (mutation_arg))
    { mutations << checkPtr (mutation_arg); }
#if 0
  SeqChange (const Alignment* al_arg,
             size_t targetStopPos);    
#endif
  void qc () const override;
  void saveText (ostream &os) const override
    { os        << start + 1 
         << ' ' << len 
         << ' ' << strQuote (reference) << " -> " << strQuote (allele)
         << ' ' << start_ref + 1 << ".." << stop_ref
         << ' ' << start_target + 1 
         << ' ' << neighborhoodMismatch;
    //if (mutation)
      for (const AmrMutation* mutation : mutations)
      { os << ' ' ;
        mutation->saveText (os);
      }
      os << endl; 
    }
  bool empty () const override
    { return ! len; }
    
    
  bool hasMutation () const
    { return ! mutations. empty () && ! empty () && ! replacement; }
  string getMutationStr () const;
  size_t getStop () const
    { return start + len; }
  bool operator< (const SeqChange &other) const;
  bool better (const SeqChange &other) const;
  bool finish (const string &refSeq,
               size_t flankingLen);
    // Return: good match
    // Input: flankingLen: valid if > 0
private:
  void setSeq ();
  void setStartStopRef ();
  void setStartTarget ();
  void setNeighborhoodMismatch (size_t flankingLen);
public:
  bool matchesMutation (const AmrMutation& mut) const;
};



void normalizeSeqs (string &seq1,
                    string &seq2);



struct Alignment : Root
// No TBLASTX
{
  // Positions are 0-based
  // start < end <= len
  
  // BLAST query
  bool targetProt {false};
    // false <=> DNA
  string targetName; 
  string targetSeq;
    // Uppercase
  size_t targetStart {0};
  size_t targetEnd {0};
  size_t targetLen {0};
  bool targetStrand {true}; 
    // false <=> negative  
  
  // BLAST subject
  bool refProt {false};
    // false <=> DNA
  string refName; 
  string refSeq;
    // Uppercase
  size_t refStart {0};
  size_t refEnd {0};
  size_t refLen {0};  
  AmrMutation refMutation;
    // !empty() => original refSeq is the result of refMutation.apply() to the original refSeq
//int ref_offset {0};
  
  // Alignment
  bool alProt {false};
  size_t nident {0};
  size_t al2ref_len {1};
  size_t al2target_len {1};

  Vector<SeqChange> seqChanges;

  
  Alignment (const string &line,
             bool targetProt_arg,
             bool refProt_arg);
  static constexpr const char* format {"qseqid sseqid qstart qend qlen sstart send slen sseq"};
    // 1-based
  Alignment () = default;
protected:
  void set_nident ();
    // Output: nident
  void refMutation2refSeq ();
    // Update: refSeq
  void setSeqChanges (const Vector<AmrMutation> &refMutations,
                      size_t flankingLen/*,
                      bool allMutationsP*/);
    // Input: flankingLen: valid if > 0
public:
  bool empty () const override
    { return targetName. empty (); }
  void qc () const override;
  void saveText (ostream &os) const override
    { os         << targetProt
         << '\t' << targetName 
         << '\t' << targetStart
         << '\t' << targetEnd 
         << '\t' << targetLen
         << '\t' << targetStrand
         << '\t' << refName
         << '\t' << refStart
         << '\t' << refEnd
         << '\t' << refLen
         << '\t' << nident
         << '\t' << targetSeq
         << '\t' << refSeq
         << endl;
    }


  bool hasMutation () const
    { for (const SeqChange& seqChange : seqChanges)
        if (seqChange. hasMutation ())
          return true;
      return false;
    }
  double pIdentity () const
    { return (double) nident / (double) targetSeq. size (); }
  double refCoverage () const
    { return (double) (refEnd - refStart) / (double) refLen; }
  double targetCoverage () const
    { return targetProt ? (double) (targetEnd - targetStart) / (double) targetLen : NaN; }
  size_t targetTail (bool upstream) const
    { return targetStrand == upstream ? targetStart : (targetLen - targetEnd); }
  bool refExactlyMatched () const
    { return    refProt
             && refLen   
             && nident == refLen 
             && nident == targetSeq. size ();
	  }
	bool getFrameShift (const Alignment &other,
	                    size_t diff_max) const
	  // Return: success
	  // Input: diff_max: in bp
	  // Requires: !targetProt, refProt, rightPart.refProt
	  { return    nident >= other. nident
	           && (          getFrameShift_right (other, diff_max) 
	               || other. getFrameShift_right (*this, diff_max) 
	              );
	  }
private:
	bool getFrameShift_right (const Alignment &rightPart,
	                          size_t diff_max) const;
};




}  // namespace


