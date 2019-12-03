// point_mut.cpp

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
*   Point mutations library
*
*/
   
   

#undef NDEBUG 
#include "common.inc"

#include "point_mut.hpp"




namespace PointMut_sp
{



namespace
{

void normalizeAlignment_ (const string &seq1,
                          string &seq2)
{
  ASSERT (seq1. size () == seq2. size ());
  
  size_t start = NO_INDEX;
  FFOR (size_t, i, seq1. size ())
    if (start == NO_INDEX)
    {
      if (seq2 [i] == '-')
        start = i;
    }
    else
      if (seq2 [i] != '-')
      {
        if (   seq1 [i] == seq2 [i]
            && seq2 [i] == seq1 [start]
           )
          swap (seq2 [i], seq2 [start]);
        else
          start = NO_INDEX;
      }

  ASSERT (seq1. size () == seq2. size ());
}

}



void normalizeAlignment (string &seq1,
                         string &seq2)
{
  normalizeAlignment_ (seq1, seq2);
  normalizeAlignment_ (seq2, seq1);
}




// PointMut


PointMut::PointMut (size_t pos_arg,
            				const string &geneMutation_arg,
            				const string &class_arg,
            				const string &subclass_arg,
            				const string &name_arg)
: pos (pos_arg)
, geneMutation (geneMutation_arg)
, classS (class_arg)
, subclass (subclass_arg)
, name (name_arg)
{ 
  QC_ASSERT (pos > 0);
	pos--;
	
	QC_ASSERT (! name. empty ());
  QC_ASSERT (! contains (name, '\t'));
  replace (name, '_', ' ');
  QC_ASSERT (! contains (name, "  "));
  
  // reference, allele
	parse (geneMutation, reference, allele);
	if (allele == "STOP")
	  allele = "*";
	else if (allele == "del")
	  allele. clear ();
	QC_ASSERT (isUpper (reference));
	QC_ASSERT (isUpper (allele));
	QC_ASSERT (reference != allele);
}



PointMut::PointMut (const string &gene,
                    size_t pos_arg,
                    const string &reference_arg,
                    const string &allele_arg)
: pos (pos_arg)
, geneMutation (gene + pm_delimiter + reference_arg + to_string (pos + 1) + allele_arg)
, name (string (reference_arg == allele_arg ? "wildtype" : "mutation") + " " + strUpper1 (gene))
, additional (true)
, reference (reference_arg)
, allele (allele_arg)
{ 
  ASSERT (! gene. empty ());
	ASSERT (isUpper (reference));
	ASSERT (isUpper (allele));
}



void PointMut::parse (const string &geneMutation,
                      string &reference,
                      string &allele)
{ 
  QC_ASSERT (! geneMutation. empty ());
  QC_ASSERT (reference. empty ());
  QC_ASSERT (allele. empty ());
  enum Type {inAllele, inPos, inRef};
  Type type = inAllele;
  FOR_REV (size_t, i, geneMutation. size ())
  {
    const char c = geneMutation [i];
    switch (type)
    {
      case inAllele:
        if (isAlpha (c))
          allele += c;
        else
          type = inPos;
        break;
      case inPos:
        if (isAlpha (c))
        {
          type = inRef;
          reference = string (1, c);
        }
        break;
      case inRef:
        if (isAlpha (c))
          reference += c;
        else
        {
          reverse (allele);
          reverse (reference);
          return;
        }
        break;
    }
  }
  NEVER_CALL;
}



string PointMut::getResistance () const
{ 
  size_t p = name. find (' ');
	ASSERT (p != string::npos);
	p = name. find (' ', p + 1);
	ASSERT (p != string::npos);
	return name. substr (p + 1);
}



bool PointMut::better (const PointMut &other) const  
{ 
  if (geneMutation != other. geneMutation)
    return false;
/*if (pos != other. pos)
    return false; */
  LESS_PART (*this, other, neighborhoodMismatch);
  // Tie
  LESS_PART (*this, other, name); 
  LESS_PART (*this, other, geneMutation);  
  return false;
}



bool PointMut::operator< (const PointMut &other) const
{ 
  LESS_PART (*this, other, pos);
  LESS_PART (*this, other, geneMutation);
  return false;
}




// SeqChange


void SeqChange::setSeq (const string &targetSeq,
                        const string &refSeq)
{ 
  ASSERT (targetSeq. size () == refSeq. size ());
  ASSERT (start + len < targetSeq. size ());
  
  reference = refSeq.    substr (start, len);
  allele    = targetSeq. substr (start, len);
  
  replaceStr (reference, "-", "");
  replaceStr (allele,    "-", "");
  
  QC_ASSERT (reference != allele);
  ASSERT (isUpper (reference));
  ASSERT (isUpper (allele));
}



void SeqChange::setRefStart (const string &refSeq, 
                             size_t refStart_arg)
{
  ASSERT (start < refSeq. size ());
  
  refStart = refStart_arg;
  FOR (size_t, i, start)
    if (refSeq [i] != '-')
      refStart++;
}



void SeqChange::setTargetStart (const string &targetSeq, 
                                size_t targetStart_arg, 
                                bool targetStrand)
{
  ASSERT (start < targetSeq. size ());
  
  targetStart = targetStart_arg;
  if (targetStrand)
  {
    FOR (size_t, i, start)
      if (targetSeq [i] != '-')
        targetStart++;
  }
  else
  {
    FOR_REV_END (size_t, i, targetSeq. size (), start + 1)
      if (targetSeq [i] != '-')
        targetStart++;
  }  
}



}  // namespace


