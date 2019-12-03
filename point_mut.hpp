// point_mut.hpp

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
*   Point mutations library. Protein or DNA point mutations.
*
*/
   
   
#include "common.hpp"
using namespace Common_sp;


   

namespace PointMut_sp
{



static constexpr char pm_delimiter = '_';



void normalizeAlignment (string &seq1,
                         string &seq2);



struct PointMut
{
	size_t pos {0};
	  // In whole reference sequence
	  // = start of reference
	  // >= 0
	// !empty()
	string geneMutation;
	  // Depends on the above
	string classS;
	string subclass;
	string name;
	  // Species binomial + resistance
	bool additional {false};
	
	// Replacement
  // Upper-case
  string reference;
	string allele;

  // Target information
	double neighborhoodMismatch {0.0};
	size_t targetStart {0};

	
	PointMut (size_t pos_arg,
						const string &geneMutation_arg,
						const string &class_arg,
						const string &subclass_arg,
						const string &name_arg);
	PointMut (const string &gene,
	          size_t pos_arg,
	          const string &reference_arg,
	          const string &allele_arg);
	PointMut () = default;
private:
	static void parse (const string &geneMutation,
	                   string &reference,
	                   string &allele);
public:


  bool empty () const
    { return name. empty (); }
  void print (ostream &os) const
    { os << pos + 1 << ' ' << geneMutation << ' ' << name << endl; }
  string getResistance () const;
  bool better (const PointMut &other) const;
  bool operator< (const PointMut &other) const;
  bool operator== (const PointMut &other) const
    { return geneMutation == other. geneMutation; }
};



struct SeqChange
{
  size_t start {0};
  size_t len {0};
  string reference;
  string allele;

  size_t refStart {0};
  size_t targetStart {0};
	double neighborhoodMismatch {0.0};
  
  
  SeqChange () = default;
  void saveText (ostream &os) const
    { os        << start + 1 
         << ' ' << len 
         << ' ' << reference << " -> " << allele 
         << ' ' << refStart + 1 
         << ' ' << targetStart + 1 
         << ' ' << neighborhoodMismatch
         << endl; 
    }
  
  
  size_t getStop () const
    { return start + len; }
  void setSeq (const string &targetSeq,
               const string &refSeq);
  void setRefStart (const string &refSeq, 
                    size_t refStart_arg);
  void setTargetStart (const string &targetSeq, 
                       size_t targetStart_arg, 
                       bool targetStrand);
};




}  // namespace


