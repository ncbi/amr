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



struct Cds 
{
  string contig;
  size_t start {0};
  size_t stop {0};
  bool strand {false};
  size_t crossOriginSeqLen {0};
    // !0 => cross origin
    // To be set externally
  
  Cds (const string &contig_arg,
       size_t start_arg,
       size_t stop_arg,
       bool strand_arg,
       size_t crossOriginSeqLen);
  Cds ()
    {} 
    
  void print (ostream &os) const
    { os << contig << ' ' << start << ' ' << stop << ' ' << strand << ' ' << crossOriginSeqLen << endl; }
  bool operator< (const Cds& other) const;
};



struct Gff : Root
// https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
// Requirement: the protein id should be in the attribute "Name=<id>" (9th field) of the rows with type "CDS" or "gene" (3rd field)
{	
  map<string, Set<Cds> > seqid2cdss; 

  Gff (const string &fName,
       bool locus_tag);
    // Input: fName may be empty
};



}



#endif
