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
    // DNA FASTA id
  size_t start {0};
  size_t stop {0};
    // start <= stop
  bool strand {false};
  bool partial {false};
  size_t contigLen {0};
    // 0 <=> unknown
  bool crossOrigin {false};
  string gene;
  string product;

  
  Locus (size_t lineNum_arg,
         const string &contig_arg,
         size_t start_arg,
         size_t stop_arg,
         bool strand_arg,
         bool partial_arg,
         size_t crossOriginSeqLen,
         string gene_arg,
         string product_arg);
  Locus () = default;
    

  bool empty () const
    { return contig. empty (); }
  void print (ostream &os) const
    { os        << contig 
         << ' ' << start 
         << ' ' << stop 
         << ' ' << strand 
         << ' ' << contigLen 
         << ' ' << crossOrigin 
         << ' ' << gene 
         << ' ' << product 
         << endl; 
    }
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



struct Gff
{
  enum Type  {bakta, genbank, microscope, patric, pgap, prodigal, prokka, pseudomonasdb, rast, standard/*PD-4548*/};
    // Alphabetic order
  static const StringVector names;
  static Type name2type (const string &name);
};



struct Annot : Root
{	
  // Protein GFF id is a function of attributes (column in GFF)
  map<string/*protein GFF id*/,Set<Locus>> prot2loci; 
  map<string/*protein FASTA id*/,string/*protein GFF id*/> fasta2gff_prot;  
    // empty() => protein FASTA id = protein GFF id


  Annot (const string &fName,
         Gff::Type gffType,
         bool protMatch,
         bool lcl);
    // GFF
		// https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
		// https://github.com/ncbi/amr/issues/91
		// Input: protMatch: property of protein FASTA:  
		//                     genbank: "[locus_tag=...]" in comment
		//                     microscope: "><acc>|ID:<num>|<gene>|
		//                     prodigal: "ID=" in comment
		//        lcl: property of DNA FASTA: >lcl|...
    /*
      gffType        protein GFF id
      -------        --------------
      bakta          ID=       
      genbank        locus_tag=[project:]acc  // if pseudo or protMatch
                     Name=[project:]acc       // else                    
      microscope     ID=
      patric         ID=...;locus_tag=...
      pgap           Name=    
      prodigal       ID=            
      prokka         ID=      
      pseudomonasdb  Alias= (or locus=)
      rast           ID=       
      standard       Name=         
        
                     ["]acc["]
    */
  explicit Annot (const string &fName);
    // Bed
		// https://genome.ucsc.edu/FAQ/FAQformat.html#format1
		
		
  void load_fasta2gff_prot (const string &fName);
    // Input: fName: file is created by gff_check.cpp -gff_prot_match
    // Output: fasta2gff_prot
  void load_fasta2gff_dna (const string &fName);
    // Input: fName: file is created by gff_check.cpp -gff_dna_match
    // Output: Locus::contig
  const Set<Locus>& findLoci (const string &fasta_prot) const;
    // Return: !empty()
    // throw if not found
};




}



#endif
