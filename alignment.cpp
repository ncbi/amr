// alignment.cpp

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
   
   

#undef NDEBUG 
#include "common.inc"

#include "alignment.hpp"




namespace Alignment_sp
{



// Mutation


Mutation::Mutation (size_t pos_arg,
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



void Mutation::parse (const string &geneMutation,
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



bool Mutation::operator< (const Mutation &other) const
{ 
  LESS_PART (*this, other, pos);
  LESS_PART (*this, other, geneMutation);
  return false;
}




// SeqChange


bool SeqChange::operator< (const SeqChange &other) const
{ 
  LESS_PART (*this, other, start_target);  
  return false;
}



void SeqChange::setSeq ()
{ 
  ASSERT (al);
  ASSERT (al->targetSeq. size () == al->refSeq. size ());
  ASSERT (start + len < al->targetSeq. size ());
  
  reference = al->refSeq.    substr (start, len);
  allele    = al->targetSeq. substr (start, len);
  
  replaceStr (reference, "-", "");
  replaceStr (allele,    "-", "");
  
  QC_ASSERT (reference != allele);
  ASSERT (isUpper (reference));
  ASSERT (isUpper (allele));
}



void SeqChange::setStartStopRef ()
{
  ASSERT (al);
  ASSERT (start < al->refSeq. size ());
  
  start_ref = al->refStart;
  FOR (size_t, i, start)
    if (al->refSeq [i] != '-')
      start_ref++;
      
  stop_ref = start_ref;
  FOR_START (size_t, i, start, start + len)
    if (al->refSeq [i] != '-')
      stop_ref++;
}



void SeqChange::setStartTarget ()
{
  ASSERT (al);
  ASSERT (start < al->targetSeq. size ());
  
  start_target = al->targetStart;
  if (al->targetStrand)
  {
    FOR (size_t, i, start)
      if (al->targetSeq [i] != '-')
        start_target++;
  }
  else
  {
    FOR_REV_END (size_t, i, start + len, al->targetSeq. size ())
      if (al->targetSeq [i] != '-')
        start_target++;
  }  
}



void SeqChange::setNeighborhoodMismatch (size_t flankingLen) 
{ 
  ASSERT (al);
  ASSERT (start < al->targetSeq. size ())
  
    
  // PD-2001
  size_t span = 0;
  size_t mismatches = 0; 
  
  size_t j = 0;  // auxiliary

  // Left flank
  if (start)
  {
    j = start - 1; 
    while (start - j <= flankingLen)
    { 
      span++;
      if (al->targetSeq [j] != al->refSeq [j])
        mismatches++;
      if (j == 0)
        break;
      j--;
    }
  }
  ASSERT (start >= j);
  const size_t missed_left = min (min (al->targetTail (true), al->refStart), flankingLen + 1 - (start - j)); 
  span       += missed_left;
  mismatches += missed_left; 
  if (verbose ())
    cout << "missed_left: " << start  << ' ' << j << ' ' << missed_left << endl;

  // Right flank
  const size_t stop = getStop ();
  for (j = stop + 1; j < al->targetSeq. size () && j - stop <= flankingLen; j++)
  { 
    span++;
    if (al->targetSeq [j] != al->refSeq [j])
      mismatches++;
  }
  ASSERT (j >= stop);
  const size_t missed_right = min (min (al->targetTail (false), al->refLen - al->refEnd), flankingLen + 1 - (j - stop));
  span       += missed_right;
  mismatches += missed_right;     
  if (verbose ())
    cout << "missed_right: " << j << ' ' << missed_right << ' ' << mismatches << ' ' << span << endl;  
    
  ASSERT (mismatches <= span);
  
  neighborhoodMismatch = (double) mismatches / (double) span;
}



bool SeqChange::matchesMutation (const Mutation& mut) const
{ 
  if (   mut. pos        < start_ref
      || mut. getStop () > stop_ref
     )
    return false;
    
  const size_t head = mut. pos - start_ref;
  const size_t tail = stop_ref - mut. getStop ();
  ASSERT (head + mut. reference. size () + tail == reference. size ());
  ASSERT (reference. substr (head, mut. reference. size ()) == mut. reference);
  
  // PD-3267 
  // ??
  if (head + mut. allele. size () + tail != allele. size ())
    return false;   
  return allele. substr (head, mut. allele. size ()) == mut. allele;
}




// Alignment

namespace
{

bool normalizeSeq (const string &seq1,
                   string &seq2)
// Return: seq2 is changed
{
  ASSERT (seq1. size () == seq2. size ());
  
  bool changed = false;
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
        {
          swap (seq2 [i], seq2 [start]);
          changed = true;
        }
        else
          start = NO_INDEX;
      }
  ASSERT (seq1. size () == seq2. size ());
  
  return changed;
}

}



void normalizeSeqs (string &seq1,
                    string &seq2)
{
  while (  normalizeSeq (seq1, seq2)
         | normalizeSeq (seq2, seq1)
        );
  ASSERT (seq1. size () == seq2. size ());
  
  while (   seq1. back () == '-'
         && seq2. back () == '-'
        )
  {
    seq1. erase (seq1. size () - 1, 1);
    seq2. erase (seq2. size () - 1, 1);
  }
  ASSERT (seq1. size () == seq2. size ());
}




namespace 
{

char complementaryNucleotide (char wildNucleotide)
{
  char r = ' ';
  switch (toLower (wildNucleotide))
  {
    case 'a': r = 't'; break;
    case 'c': r = 'g'; break;
    case 'g': r = 'c'; break;
    case 't': r = 'a'; break;
    case 'm': r = 'k'; break;
    case 'r': r = 'y'; break;
    case 'w': r = 'w'; break;
    case 's': r = 's'; break;
    case 'y': r = 'r'; break;
    case 'k': r = 'm'; break;
    case 'v': r = 'b'; break;
    case 'h': r = 'd'; break;
    case 'd': r = 'h'; break;
    case 'b': r = 'v'; break;
    case 'n': r = 'n'; break;
    default: 
    	throw runtime_error ("Bad wild nucleotide " + toString (wildNucleotide));
  }
  if (isupper (wildNucleotide))
    r = toUpper (r);

  return r;
}



void reverseDna (string &seq)
{
  const size_t len = seq. size ();
  if (! len)
  	return;

  FFOR (size_t, i, len / 2)
  {
  	const size_t j = len - 1 - i;
  	if (seq [i] != '-')
  		seq [i] = complementaryNucleotide (seq [i]);
  	if (seq [j] != '-')
  		seq [j] = complementaryNucleotide (seq [j]);
    swap (seq [i], seq [j]);
  }

  if (len % 2)
  {
  	const size_t i = len / 2;
  	if (seq [i] != '-')
  		seq [i] = complementaryNucleotide (seq [i]);
  }
}

}



Alignment::Alignment (const string &line,
                      bool targetProt_arg,
                      bool refProt_arg)
: targetProt (targetProt_arg)
, refProt (refProt_arg)
{
  try
  {
    istringstream iss (line);
    iss >> targetName >> refName >> targetStart >> targetEnd >> targetLen >> refStart >> refEnd >> refLen >> targetSeq >> refSeq;
    strUpper (targetSeq);
    strUpper (refSeq);
    
    QC_ASSERT (! refSeq. empty ());	
    QC_ASSERT (targetSeq. size () == refSeq. size ());    

    if (! targetProt && ! refProt)
    {
      QC_ASSERT (refStart != refEnd);
      const bool refStrand = refStart < refEnd;  
      if (! refStrand)
      {
        swap (refStart, refEnd);
      	reverseDna (targetSeq);
      	reverseDna (refSeq);
      	toggle (targetStrand);
      }
    }

    if (! targetProt && refProt)
    {
	    QC_ASSERT (targetStart != targetEnd);
	    targetStrand = targetStart < targetEnd;  
	    if (! targetStrand)
	      swap (targetStart, targetEnd);
	  }
    
    QC_ASSERT (refStart >= 1);
    QC_ASSERT (targetStart >= 1);
    refStart--;
    targetStart--;

    normalizeSeqs (targetSeq, refSeq);
    
    FFOR (size_t, i, targetSeq. size ())
     if (targetSeq [i] == refSeq [i])
       nident++;
  }
  catch (...)
  {
  	cout << line << endl;
  	throw;
  }
}



void Alignment::setSeqChanges (const Vector<Mutation> &refMutations,
                               size_t flankingLen,
                               bool allMutationsP)
{
  ASSERT (seqChanges. empty ());
  ASSERT (! refMutations. empty ());  	      
  
  // pmGene
  string s (refMutations. front (). geneMutation);
  rfindSplit (s, pm_delimiter);  // Only gene symbol
  const string pmGene (s);

  {
	  SeqChange seqChange (this);
	  bool inSeqChange = false;
		FFOR (size_t, i, refSeq. size ())
		  if (inSeqChange)
		  {
        if (targetSeq [i] == refSeq [i])
        {
          seqChange. setSeq ();
          if (refSeq [seqChange. start] == '-')
          {
            QC_ASSERT (seqChange. start);
            seqChange. start--;
            seqChange. len++;
	          seqChange. setSeq ();
          }
          ASSERT (! seqChange. reference. empty ());
          seqChange. setStartStopRef ();
          seqChange. setStartTarget ();
          if (flankingLen)
    	      seqChange. setNeighborhoodMismatch (flankingLen);
    	    if (verbose ())
    	      seqChange. saveText (cout);
    			if (seqChange. neighborhoodMismatch <= 0.04)   // PAR  // PD-3191
            seqChanges << move (seqChange);
          seqChange = SeqChange (this);
          inSeqChange = false;
        }
        else
          seqChange. len++;
      }
      else
        if (targetSeq [i] != refSeq [i])
        {
          seqChange. start = i;
          seqChange. len = 1;
          inSeqChange = true;
        }
  }
  if (verbose ())
    PRINT (seqChanges. size ());
  
	size_t j = 0;
  while (j < refMutations. size () && refMutations [j]. pos < refStart)
  {
  	if (allMutationsP)
  	  seqChanges << SeqChange (this, & refMutations [j]);
    j++;
  }
  
	size_t start_ref_prev = NO_INDEX;
  for (SeqChange& seqChange : seqChanges)
  {
    if (verbose ())
      seqChange. saveText (cout);
    QC_IMPLY (start_ref_prev != NO_INDEX, start_ref_prev < seqChange. start_ref);
    while (j < refMutations. size ())
    {
		  const Mutation& mut = refMutations [j];
		  if (verbose ())
		  {
		    cout << "Processing Mutation: ";
		    mut. print (cout);
		  }
		  if (mut. pos < seqChange. start_ref)
	    	if (allMutationsP)
      	  seqChanges << SeqChange (this, & mut);
		  if (mut. getStop () > seqChange. stop_ref)
		    break;
		  if (seqChange. matchesMutation (mut))
		  {
		    seqChange. mutation = & mut;
		  	if (verbose ())
		  	{
		  	  cout << "Found: ";
		  	  mut. print (cout);
		  	}
		  }
	    j++;
		}
	  start_ref_prev = seqChange. start_ref;
	}
    
	if (allMutationsP)
	  while (j < refMutations. size ())
	  {
  	  seqChanges << SeqChange (this, & refMutations [j]);
	    j++;
	  }
	
  for (Iter<Vector<SeqChange>> it (seqChanges); it. next ();)
    if (! it->mutation)
      it. erase ();	
	
	seqChanges. sort ();
}



void Alignment::qc () const
{
  if (! qc_on)
    return;
  
  QC_ASSERT (targetSeq. size () == refSeq. size ());  

  QC_ASSERT (targetStart <= targetEnd);
  QC_ASSERT (targetEnd <= targetLen);
  QC_IMPLY (! targetSeq. empty (), targetEnd - targetStart <= targetSeq. size ());
  QC_IMPLY (targetProt, targetStrand);

  QC_ASSERT (refStart <= refEnd);
  QC_ASSERT (refEnd <= refLen);
  QC_IMPLY (! refSeq. empty (), refEnd - refStart <= refSeq. size ());

  QC_ASSERT (nident <= refEnd    - refStart);
  QC_ASSERT (nident <= targetEnd - targetStart);
}




}  // namespace


