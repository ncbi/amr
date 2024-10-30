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

#include "alignment.hpp"

#include "common.inc"



namespace Alignment_sp
{



// AmrMutation

AmrMutation::AmrMutation (size_t pos_real_arg,
            		 				  const string &geneMutation_std_arg,
            		 				  const string &geneMutation_arg,
                  				const string &class_arg,
                  				const string &subclass_arg,
                  				const string &name_arg)
: pos_real (pos_real_arg)
, geneMutation_std (geneMutation_std_arg)
, geneMutation     (geneMutation_arg)
, classS (class_arg)
, subclass (subclass_arg)
, name (name_arg)
{ 
  QC_ASSERT (pos_real > 0);
	pos_real--;
	
	QC_ASSERT (! name. empty ());
  QC_ASSERT (! contains (name, '\t'));
  replace (name, '_', ' ');
  QC_ASSERT (! contains (name, "  "));
  
  // reference, allele
	parse (geneMutation_std, reference, allele, gene, pos_std, frameshift, frameshift_insertion);
	if (allele == "STOP")
	  allele = "*";
	else if (allele == "del")
	  allele. clear ();
	  
	qc ();
}



void AmrMutation::parse (const string &geneMutation_std,
                         string &reference,
                         string &allele,
                         string &gene,
                         int &pos_std,
                         size_t &frameshift,
                         int &frameshift_insertion)
{ 
  QC_ASSERT (! geneMutation_std. empty ());
  QC_ASSERT (reference. empty ());
  QC_ASSERT (allele. empty ());
  QC_ASSERT (frameshift == no_index);
  QC_ASSERT (! frameshift_insertion);
  
  
  string pureMutation (geneMutation_std);

  static const string fsInfix ("fsTer");  // ...fsTer<N>{ins|del}<M>
  const size_t fsInfix_pos = geneMutation_std. find (fsInfix);
  if (fsInfix_pos != no_index)
  {
    pureMutation. erase (fsInfix_pos);
    string suffix = geneMutation_std. substr (fsInfix_pos + fsInfix. size ());  // <N>{ins|del}<M>
    size_t indel_pos = suffix. find ("ins");
    bool ins = true;  // {ins|del}
    if (indel_pos == string::npos)
    {
      indel_pos = suffix. find ("del");
      ins = false;
    }
    QC_ASSERT (indel_pos != string::npos);      
    frameshift_insertion = str2<int> (suffix. substr (indel_pos + 3));  // <M>
    if (! ins)
      frameshift_insertion *= -1;
    QC_ASSERT (frameshift_insertion % 3);
    suffix. erase (indel_pos);  // <N>
    frameshift = str2<size_t> (suffix);
    QC_ASSERT (frameshift != no_index);
  }

  enum Type {inAllele, inPos, inRef};
  Type type = inAllele;
  string pos_stdS;
  FOR_REV (size_t, i, pureMutation. size ())
  {
    const char c = pureMutation [i];
    switch (type)
    {
      case inAllele:
        if (isAlpha (c))
          allele += c;
        else
        {
          pos_stdS = string (1, c);
          type = inPos;
        }
        break;
      case inPos:
        if (isAlpha (c))
        {
          type = inRef;
          reference = string (1, c);
        }
        else
          pos_stdS += c;
        break;
      case inRef:
        if (isAlpha (c))
          reference += c;
        else
        {
          QC_ASSERT (c == '_');
          QC_ASSERT (i);
          Common_sp::reverse (allele);
          Common_sp::reverse (reference);
          Common_sp::reverse (pos_stdS);
          gene = pureMutation. substr (0, i);
          pos_std = stoi (pos_stdS) - 1;
          return;
        }
        break;
    }
  }
  NEVER_CALL;
}



void AmrMutation::qc () const
{
  if (! qc_on)
    return;
    
  if (empty ())
    return;
    
  QC_ASSERT (! geneMutation. empty ());
	QC_ASSERT (isUpper (reference));
	QC_ASSERT (isUpper (allele));
	QC_ASSERT (reference != allele);
	QC_ASSERT (! contains (reference, '-'));
	QC_ASSERT (! contains (allele,    '-'));
	QC_ASSERT (! gene. empty ());
	QC_ASSERT (! reference. empty ());  
	QC_IMPLY (frameshift != no_index, ! allele. empty ());
	QC_IMPLY (frameshift != no_index, pos_std >= 0);
	QC_ASSERT ((frameshift == no_index) == (frameshift_insertion == 0));
	QC_IMPLY (frameshift_insertion, frameshift_insertion % 3);
}



bool AmrMutation::operator< (const AmrMutation &other) const
{ 
  LESS_PART (*this, other, gene);
  LESS_PART (*this, other, pos_real);
  LESS_PART (*this, other, geneMutation_std);
  return false;
}




// SeqChange

void SeqChange::qc () const
{
  if (! qc_on)
    return;
         
  QC_ASSERT (al);
  
  if (empty ())
  {
    QC_ASSERT (! mutations. empty () /*mutation*/);
    return;
  }
  
  QC_ASSERT (len);
  QC_ASSERT (! reference. empty ());  // For an insertion the preceding character is artifically added
  QC_ASSERT (start_ref < stop_ref);
  QC_ASSERT (stop_ref <= al->refEnd);
  QC_ASSERT (reference. size () * al->al2ref_len <= stop_ref - start_ref);
  QC_ASSERT (start_target <= al->targetEnd);
  QC_ASSERT (! contains (reference, '-'));
  QC_ASSERT (! contains (allele, '-'));
  QC_ASSERT (reference != allele);
  QC_ASSERT (isUpper (reference));
  QC_ASSERT (isUpper (allele));
  for (const AmrMutation* mutation : mutations)
  {
    QC_ASSERT (mutation);
    QC_ASSERT (between (mutation->pos_real, al->refStart, al->refEnd));
    QC_IMPLY (mutation->frameshift != no_index, al->refProt);
    QC_IMPLY (mutation->frameshift != no_index, ! al->targetProt);
  }
}



string SeqChange::getMutationStr () const
{ 
  const string allele_ (allele. empty () 
                          ? "DEL" 
                          : allele == "*"
                              ? "STOP"
                              : allele
                       );
  return reference + to_string ((int) start_ref + 1 /* - al->ref_offset*/) + allele_; 
}



bool SeqChange::operator< (const SeqChange &other) const
{ 
  LESS_PART (*this, other, start_target);  
  return false;
}



bool SeqChange::better (const SeqChange &other) const
{ 
  LESS_PART (other, *this, hasMutation ());
  LESS_PART (*this, other, neighborhoodMismatch);  
  LESS_PART (other, *this, al->pIdentity ());  
  return false;
}



bool SeqChange::finish (const string &refSeq,
                        size_t flankingLen)
{
  setSeq ();
  if (refSeq [start] == '-')
  {
    QC_ASSERT (start);
    start--;
    len++;
    ASSERT (refSeq [start] != '-');
    setSeq ();
  }
  ASSERT (! reference. empty ());
  return finishPos (flankingLen);
}  



bool SeqChange::finishPos (size_t flankingLen)
{
  setStartStopRef ();
  setStartTarget ();
  setNeighborhoodMismatch (flankingLen);
  ASSERT (! reference. empty ());
  if (verbose ())
    saveText (cout);
	return neighborhoodMismatch <= 0.04;   // PAR  // PD-3191
}  



void SeqChange::setSeq ()
{ 
  ASSERT (al);
  ASSERT (al->targetSeq. size () == al->refSeq. size ());
  ASSERT (start + len <= al->targetSeq. size ());
  
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
      start_ref += al->al2ref_len;
      
  stop_ref = start_ref;
  FOR_START (size_t, i, start, start + len)
    if (al->refSeq [i] != '-')
      stop_ref += al->al2ref_len;
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
        start_target += al->al2target_len;
  }
  else
  {
    FOR_REV_END (size_t, i, start + len, al->targetSeq. size ())
      if (al->targetSeq [i] != '-')
        start_target += al->al2target_len;
  }  
}



void SeqChange::setNeighborhoodMismatch (size_t flankingLen) 
{ 
  ASSERT (al);
  ASSERT (start < al->targetSeq. size ())
  

  neighborhoodMismatch = 0.0;
  if (! flankingLen)
    return;
  
    
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



bool SeqChange::matchesMutation (const AmrMutation& mut) const
{ 
  if (empty ())
    return false;
    
  if (   mut. pos_real        < start_ref
      || mut. getStop () > stop_ref
     )
    return false;
    
  const size_t head = mut. pos_real - start_ref;
  const size_t tail = stop_ref - mut. getStop ();
  ASSERT (head + mut. reference. size () + tail == reference. size ());
  const string ref_seg (reference. substr (head, mut. reference. size ()));
  if (ref_seg != mut. reference)
    throw runtime_error ("Reference sequence has " + strQuote (ref_seg) + ", but mutation is: " + mut. geneMutation_std);
  
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
  size_t start = no_index;
  // -AAA --> AAA-
  FFOR (size_t, i, seq1. size ())
    if (start == no_index)
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
          start = no_index;
      }
  ASSERT (seq1. size () == seq2. size ());
  
  return changed;
}

}



void normalizeSeqs (string &seq1,
                    string &seq2)
{
  for (;;)
  {
    const bool b1 = normalizeSeq (seq1, seq2);
    const bool b2 = normalizeSeq (seq2, seq1);
    if (   ! b1 
        && ! b2
       )
      break;
  }
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
, alProt (targetProt || refProt)
, al2ref_len    (alProt && ! refProt    ? 3 : 1)
, al2target_len (alProt && ! targetProt ? 3 : 1)
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

    if (refProt)
    {
      if (   refEnd == refLen 
          && refSeq. back () == '*'
         )
      {
        targetStopCodon = toEbool (targetSeq. back () == '*');
        if (targetStopCodon == etrue)
        {
          nident--;        
        //if (targetProt)
          //targetLen -= al2target_len;
        }
        do
        {
          // Target
          targetSeq. erase (targetSeq. size () - 1);
          if (targetStrand)
            targetEnd -= al2target_len;
          else
            targetStart += al2target_len;
          // Reference
          refSeq. erase (refSeq. size () - 1);
          refEnd--;
          QC_ASSERT (! refSeq. empty ());
        }
        while (targetSeq. back () != refSeq. back ());
      //seqChanges ??
      }
      else if (   refEnd == refLen - 1 
               && targetTail (false) >= al2target_len
              )
        targetStopCodon = efalse;
      refLen--;
    }

    normalizeSeqs (targetSeq, refSeq);    
    set_nident ();
  }
  catch (...)
  {
  	cout << line << endl;
  	throw;
  }
}



void Alignment::set_nident ()
{
  ASSERT (targetSeq. size () == refSeq. size ());
  nident = 0;
  FFOR (size_t, i, targetSeq. size ())
    if (targetSeq [i] == refSeq [i])
      nident++;
}



void Alignment::setSeqChanges (const Vector<AmrMutation> &refMutations,
                               size_t flankingLen/*,
                               bool allMutationsP*/)
{
  ASSERT (seqChanges. empty ());
  ASSERT (! refMutations. empty ());  	
  
  
  if (! refMutation. empty ())
  {
    if (refMutation. frameshift != no_index && ! (refProt && ! targetProt))
      return;
    {
      const size_t index = refMutations. indexOf (refMutation);
      QC_ASSERT (index != no_index);
      refMutation = refMutations [index];
    }
    const size_t start = refMutation2refSeq_pos ();
    if (start != no_index)
    {
      SeqChange seqChange (this, true);
      seqChange. start = start;
      seqChange. len = refMutation. allele. size ();
      seqChange. reference = refMutation. reference;
      seqChange. allele    = refMutation. allele;
      if (seqChange. finishPos (flankingLen))
      {
  		//ASSERT (seqChange. matchesMutation (refMutation))
  		  seqChange. mutations << & refMutation;
        seqChanges << std::move (seqChange);
      }
    }
  	if (verbose ())
  	{
  	  cout << endl;
  	  cout << targetName << " vs. " << refName << ":" << endl;
  	  cout << "Declarative mutation:" << endl;
  	  for (const SeqChange& seqChange : seqChanges)
  	    seqChange. saveText (cout);
  	}
    return;
  }

  
  {
	  SeqChange seqChange (this, false);
	  bool inSeqChange = false;
		FFOR (size_t, i, refSeq. size () + 1)
		  if (inSeqChange)
		  {
        if (targetSeq [i] == refSeq [i])
        {
          if (seqChange. finish (refSeq, flankingLen))            
            seqChanges << std::move (seqChange);
          seqChange = SeqChange (this, false);
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
#if 0
  if (   targetProt 
      && refProt
      && targetEnd == targetLen
      && refEnd < refLen
     )
  {
	  SeqChange seqChange (this, targetLen);
    seqChanges << std::move (seqChange);
  }
#endif
  if (verbose ())
  {
    PRINT (seqChanges. size ());
    cout << endl;
  }

  
//Vector<SeqChange> seqChanges_add;
	size_t j = 0;
	
  while (j < refMutations. size () && refMutations [j]. pos_real < refStart)
  {
  //if (allMutationsP)
  	//seqChanges_add << SeqChange (this, & refMutations [j]);
    j++;
  }
  
  // SeqChange::mutations
	size_t start_ref_prev = no_index;
	for (SeqChange& seqChange : seqChanges)
  {
    seqChange. qc ();
    if (verbose ())
      seqChange. saveText (cout);
    IMPLY (start_ref_prev != no_index, start_ref_prev <= seqChange. start_ref);
    while (j < refMutations. size ())
    {
		  const AmrMutation& mut = refMutations [j];
		  if (verbose ())
		  {
		    cout << "Processing AmrMutation: ";
		    mut. saveText (cout);
		    cout << endl;
		  }
	  #if 0
		  if (mut. pos_real < seqChange. start_ref)
	     if (allMutationsP)
      	  seqChanges_add << SeqChange (this, & mut);  // "seqChanges <<" destroys seqChange
    #endif
		  if (mut. pos_real >= seqChange. stop_ref)
		    break;
		#if 0
		  cout << mut << endl;
		  PRINT (seqChange. matchesMutation (mut));  
		#endif
		  if (seqChange. matchesMutation (mut))
		  {
		    seqChange. mutations << & mut;
		  	if (verbose ())
		  	{
		  	  cout << "Found: ";
		  	  mut. saveText (cout);
		  	  cout << endl;
		  	}
		  }
	    j++;
		}
	  start_ref_prev = seqChange. start_ref;
	}
	
//seqChanges << std::move (seqChanges_add);
    
#if 0
	if (allMutationsP)
	  while (j < refMutations. size ())
	  {
  	  seqChanges << SeqChange (this, & refMutations [j]);
	    j++;
	  }
#endif

  // WILDTYPE
  {
    size_t refPos = refStart;
    size_t i = 0;
    for (const AmrMutation& mut : refMutations)
    {
      while (refPos < mut. pos_real)
      {
        ASSERT (refSeq [i] != '-');
        i++;
        while (refSeq [i] == '-')
          i++;
        if (! refSeq [i])
          break;
        refPos++;
      }
      if (! refSeq [i])
        break;
      if (   refPos == mut. pos_real
          && isLeft (string (& refSeq    [i]), mut. reference)
          && isLeft (string (& targetSeq [i]), mut. reference)
         )
    	  seqChanges << SeqChange (this, & mut);
    }
  }
	
	seqChanges. sort ();
	if (verbose ())
	{
	  cout << endl;
	  cout << targetName << " vs. " << refName << ":" << endl;
	  cout << "All found mutations:" << endl;
	  for (const SeqChange& seqChange : seqChanges)
	    seqChange. saveText (cout);
	}
}



size_t Alignment::refMutation2refSeq_pos ()
{
  ASSERT (! refMutation. empty ());
  
  if (   refMutation. pos_real < refStart
      || refMutation. pos_real + refMutation. allele. size () > refEnd
     )
    return no_index;
    
  size_t pos = refStart;
  size_t frameshift_i = no_index;
	size_t i = 0;
	while (i < refSeq. size ())
	{
	  ASSERT (pos <= refMutation. pos_real);
	  if (refSeq [i] != '-')
	  {
      if (pos == refMutation. pos_real)
      {
        if (   refSeq.    substr (i, refMutation. allele. size ()) == refMutation. allele
            && targetSeq. substr (i, refMutation. allele. size ()) == refMutation. allele
            && (pos == 0 || (i > 0 && refSeq [i - 1] == targetSeq [i - 1]))
           )
        {
          if (refMutation. frameshift != no_index) 
          {
            frameshift_i = i;
            i += refMutation. allele. size ();
            pos = refMutation. getStop ();
            break;
          }
          const size_t alignStop = i   + refMutation. allele. size ();
          const size_t refStop   = pos + refMutation. allele. size ();
          ASSERT (alignStop <= refSeq. size ());
          if (refStop == refLen || refSeq [alignStop] == targetSeq [alignStop])
            return i;
        }
        return no_index;
      }
  	  pos++;
    }
    i++;
	}

	while (i < refSeq. size ())
	{
	  ASSERT (frameshift_i != no_index);
	  if (refSeq [i] != '-')
	  {
	  #if 0
	    PRINT (i);
	    PRINT (pos);
	    PRINT (refMutation. getStop ());
	    PRINT (refMutation. frameshift);
	    PRINT (refSeq    [i]);
	    PRINT (targetSeq [i]);	    
	  #endif
	    ASSERT (refMutation. getStop ());
	    if (pos == refMutation. getStop () - 1 + refMutation. frameshift)
	    {
	      if (   refSeq    [i] == '*'
	          && targetSeq [i] == '*'
	         )
	        return frameshift_i;
	      return no_index;
	    } 
  	  pos++;
    }
    i++;
	}

  return no_index;
}



void Alignment::qc () const
{
  if (! qc_on)
    return;
    
  if (empty ())
  {
    QC_ASSERT (targetName. empty ());    
    QC_ASSERT (targetSeq. empty ());
    QC_ASSERT (! targetStart);
    QC_ASSERT (! targetEnd);
    QC_ASSERT (! targetLen);
    QC_ASSERT (refName. empty ());
    QC_ASSERT (refSeq. empty ());
    QC_ASSERT (! refStart);
    QC_ASSERT (! refEnd);
    QC_ASSERT (! refLen);  
    QC_ASSERT (! nident);
    return;
  }
  
  QC_ASSERT (targetSeq. size () == refSeq. size ());  

  QC_ASSERT (targetStart <= targetEnd);
  QC_ASSERT (targetEnd <= targetLen);
  QC_ASSERT ((targetEnd - targetStart) % al2target_len == 0);
  QC_IMPLY (! targetSeq. empty (), (targetEnd - targetStart) / al2target_len <= targetSeq. size ());
  QC_IMPLY (targetProt, targetStrand);
  QC_IMPLY (targetStopCodon != enull, refProt);

  QC_ASSERT (refStart <= refEnd);
  QC_ASSERT (refEnd <= refLen);
  QC_ASSERT ((refEnd - refStart) % al2ref_len == 0);
  QC_IMPLY (! refSeq. empty (), (refEnd - refStart) / al2ref_len <= refSeq. size ());
//QC_IMPLY (refProt, ! contains (refSeq, "*"));
  
  refMutation. qc ();

  QC_ASSERT (alProt == targetProt || refProt);
  
  QC_ASSERT (nident <= (refEnd    - refStart) / al2ref_len);
  QC_ASSERT (nident <= (targetEnd - targetStart) / al2target_len);
  
  for (const SeqChange& seqChange : seqChanges)
  {
    seqChange. qc ();
    QC_ASSERT (seqChange. al == this);
  }
}



long Alignment::getGlobalTargetStart () const
{ 
  ASSERT (! targetProt);
  ASSERT (refProt);
  return targetStrand
           ? (long) targetStart - 3 * (long) refStart
           : (long) targetEnd   + 3 * (long) refStart; 
}



}  // namespace


