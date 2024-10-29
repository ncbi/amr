// columns.hpp

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
*   AMRFinderPlus column names
*
*/



// PD-5085

constexpr const char* prot_colName = "Protein id";   // PD-2534  
constexpr const char* contig_colName = "Contig id";
// Target 
constexpr const char* start_colName = "Start"; 
constexpr const char* stop_colName = "Stop";
constexpr const char* strand_colName = "Strand";
//        
constexpr const char* genesymbol_colName = "Element symbol";  // PD-4924 
constexpr const char* elemName_colName = "Element name";  // PD-4910
constexpr const char* scope_colName = "Scope";  // PD-2825 
// PD-1856
constexpr const char* type_colName = "Type";
constexpr const char* subtype_colName = "Subtype";
constexpr const char* class_colName = "Class";
constexpr const char* subclass_colName = "Subclass";
//        
constexpr const char* method_colName = "Method";
constexpr const char* targetLen_colName = "Target length";  // was: "Element length" (temporarily)
constexpr const char* refLen_colName = "Reference sequence length";
constexpr const char* refCov_colName = "% Coverage of reference";
constexpr const char* refIdent_colName = "% Identity to reference";
constexpr const char* alignLen_colName = "Alignment length";
constexpr const char* closestRefAccession_colName = "Closest reference accession";
constexpr const char* closestRefName_colName = "Closest reference name";
constexpr const char* hmmAccession_colName = "HMM accession";
constexpr const char* hmmDescr_colName = "HMM description";
constexpr const char* hierarchyNode_colName = "Hierarchy node";


// PD-5155
constexpr const char* fusion_infix = "::";  // was: "/"
  