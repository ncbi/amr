// gff.hpp

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
   

#ifndef GFF_HPP
#define GFF_HPP


#include "common.hpp"
using namespace Common_sp;



namespace GFF_sp
{



struct Locus 
{
  static constexpr size_t end_delta = 3;  // PAR
  size_t lineNum {0};
    // >= 1
    // 0 - unknown
  string contig;
  size_t start {0};
  size_t stop {0};
    // start <= stop
  bool strand {false};
  bool partial {false};
  size_t contigLen {0};
    // 0 <=> unknown
  bool crossOrigin {false};

  
  Locus (size_t lineNum_arg,
         const string &contig_arg,
         size_t start_arg,
         size_t stop_arg,
         bool strand_arg,
         bool partial_arg,
         size_t crossOriginSeqLen);
  Locus () = default;
    

  bool empty () const
    { return contig. empty (); }
  void print (ostream &os) const
    { os << contig << ' ' << start << ' ' << stop << ' ' << strand << ' ' << contigLen << ' ' << crossOrigin << endl; }
  bool operator< (const Locus& other) const;
  size_t size () const
    { return crossOrigin
               ? contigLen - stop + start
               : stop - start;            
    }
  bool atContigStart () const
    { return start <= end_delta; }
  bool atContigStop () const
    { return contigLen && contigLen - stop <= end_delta;}
};



struct Annot : Root
{	
  map<string/*protein accession*/, Set<Locus>> prot2cdss; 


  class Gff {};
  Annot (Gff,
         const string &fName,
         bool trimProject,
         bool locus_tag,
         bool pgap);
		// https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
		// Requirement: the protein id should be in the attribute "Name=<id>" (9th field) of the rows with type "CDS" or "gene" (3rd field)
    // Input: fName may be empty
  class Bed {};
  Annot (Bed,
         const string &fName);
		// https://genome.ucsc.edu/FAQ/FAQformat.html#format1
};




}



#endif
