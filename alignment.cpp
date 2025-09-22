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
	if (allele == terminatorWord)
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
    QC_ASSERT (disr. empty ());
    return;
  }
  
  disr. qc ();
  if (disr. empty ())
  {
    QC_ASSERT (len);
    QC_ASSERT (start_ref < stop_ref);
    QC_ASSERT (start_target <= al->sInt. stop);
    QC_ASSERT (! contains (reference, '-'));
    QC_ASSERT (! contains (allele, '-'));
    QC_ASSERT (isUpper (reference));
    QC_ASSERT (isUpper (allele));
    if (reference. empty ())  
    {
      ASSERT (isFrameshift ());
      QC_ASSERT (allele. empty ());
      QC_ASSERT (len == 1);
      QC_ASSERT (stop_ref == start_ref + 1);
      QC_ASSERT (! neighborhoodMismatch);
      QC_ASSERT (mutations. empty ());
    }
    else
    {
      QC_ASSERT (stop_ref <= al->qInt. stop);
      QC_ASSERT (reference. size () * al->a2q <= stop_ref - start_ref);
      QC_ASSERT (reference != allele);
    }
    for (const AmrMutation* mutation : mutations)
    {
      QC_ASSERT (mutation);
      QC_ASSERT (between (mutation->pos_real, al->qInt. start, al->qInt. stop));
      QC_IMPLY (mutation->frameshift != no_index, al->qProt);
      QC_IMPLY (mutation->frameshift != no_index, ! al->sProt);
    }
  }
  else
  {    
    QC_ASSERT (mutations. empty ());
  }
}



string SeqChange::getMutationStr () const
{ 
  ASSERT (al);
  if (! disr. empty ())
    return disr. genesymbol_raw ();
  const string allele_ (allele. empty () 
                          ? "DEL" 
                          : allele == "*"
                              ? terminatorWord
                              : allele
                       );
  return reference + to_string ((int) start_ref + 1) + allele_; 
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
  LESS_PART (other, *this, al->relIdentity ());  
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
  ASSERT (al->sseq. size () == al->qseq. size ());
  ASSERT (start <= al->sseq. size ());
  ASSERT (len <= al->sseq. size () - start);
  
  reference = al->qseq. substr (start, len);
  allele    = al->sseq. substr (start, len);
  
  replaceStr (reference, "-", noString);
  replaceStr (allele,    "-", noString);
  
  QC_ASSERT (reference != allele);
  strUpper (reference);
  strUpper (allele);
}



void SeqChange::setStartStopRef ()
{
  ASSERT (al);
  ASSERT (start < al->qseq. size ());
  
  start_ref = al->qInt. start;
  FOR (size_t, i, start)
    if (al->qseq [i] != '-')
      start_ref += al->a2q;
      
  stop_ref = start_ref;
  FOR_START (size_t, i, start, start + len)
    if (al->qseq [i] != '-')
      stop_ref += al->a2q;
}



void SeqChange::setStartTarget ()
{
  ASSERT (al);
  ASSERT (start < al->sseq. size ());
  
  start_target = al->sInt. start;
  if (al->sInt. strand == 1)
  {
    FOR (size_t, i, start)
      if (al->sseq [i] != '-')
        start_target += al->a2s;
  }
  else
  {
    FOR_REV_END (size_t, i, start + len, al->sseq. size ())
      if (al->sseq [i] != '-')
        start_target += al->a2s;
  }  
}



void SeqChange::setNeighborhoodMismatch (size_t flankingLen) 
{ 
  ASSERT (al);
  ASSERT (start < al->sseq. size ())
  

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
      if (al->sseq [j] != al->qseq [j])
        mismatches++;
      if (j == 0)
        break;
      j--;
    }
  }
  ASSERT (start >= j);
  const size_t missed_left = min (min (al->sTail (true), al->qInt. start), flankingLen + 1 - (start - j)); 
  span       += missed_left;
  mismatches += missed_left; 
  if (verbose ())
    cout << "missed_left: " << start  << ' ' << j << ' ' << missed_left << endl;

  // Right flank
  const size_t stop = getStop ();
  for (j = stop + 1; j < al->sseq. size () && j - stop <= flankingLen; j++)
  { 
    span++;
    if (al->sseq [j] != al->qseq [j])
      mismatches++;
  }
  ASSERT (j >= stop);
  const size_t missed_right = min (min (al->sTail (false), al->qlen - al->qInt. stop), flankingLen + 1 - (j - stop));
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

void Alignment::setSeqChanges (const Vector<AmrMutation> &refMutations,
                               size_t flankingLen/*,
                               bool allMutationsP*/)
{
  ASSERT (seqChanges. empty ());
  ASSERT (! refMutations. empty ());  	
  
  
  if (! refMutation. empty ())
  {
    if (   refMutation. frameshift != no_index 
        && ! blastx ()
       )
      return;
    {
      const size_t index = refMutations. indexOf (refMutation);
      QC_ASSERT (index != no_index);
      refMutation = refMutations [index];
    }
    const size_t start = refMutation2refSeq_pos ();
    if (start != no_index)
    {
      SeqChange seqChange (this/*, true*/);
      seqChange. start     = start;
      seqChange. len       = refMutation. allele. size ();
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
  	  cout << sseqid << " vs. " << qseqid << ":" << endl;
  	  cout << "Declarative mutation:" << endl;
  	  for (const SeqChange& seqChange : seqChanges)
  	    seqChange. saveText (cout);
  	}
    return;
  }

  
  {
	  SeqChange seqChange (this/*, false*/);
	  bool inSeqChange = false;
		FFOR (size_t, i, qseq. size () + 1)
		  if (inSeqChange)
		  {
        if (sseq [i] == qseq [i])
        {
          if (seqChange. finish (qseq, flankingLen))            
            seqChanges << std::move (seqChange);
          seqChange = SeqChange (this/*, false*/);
          inSeqChange = false;
        }
        else
          seqChange. len++;
      }
      else
        if (sseq [i] != qseq [i])
        {
          seqChange. start = i;
          seqChange. len = 1;
          inSeqChange = true;
        }
  }
#if 0
  if (   sProt 
      && qProt
      && send == slen
      && qend < qlen
     )
  {
	  SeqChange seqChange (this, slen);
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
	
  while (j < refMutations. size () && refMutations [j]. pos_real < qInt. start)
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
    size_t refPos = qInt. start;
    size_t i = 0;
    if (verbose ())
      cout << "refMutations: " << refMutations. size () << endl;
    for (const AmrMutation& mut : refMutations)
    {
      while (refPos < mut. pos_real)
      {
        ASSERT (qseq [i] != '-');
        i++;
        while (qseq [i] == '-')
          i++;
        if (! qseq [i])
          break;
        refPos++;
      }
      if (! qseq [i])
        break;
      if (   refPos == mut. pos_real
          && isLeft (string (& qseq [i]), mut. reference)
          && isLeft (string (& sseq [i]), mut. reference)
         )
    	  seqChanges << SeqChange (this, & mut);
    }
  }
	
	seqChanges. sort ();
	if (verbose ())
	{
	  cout << endl;
	  cout << sseqid << " vs. " << qseqid << ":" << endl;
	  cout << "All found mutations:" << endl;
	  for (const SeqChange& seqChange : seqChanges)
	    seqChange. saveText (cout);
	}
}



size_t Alignment::refMutation2refSeq_pos ()
{
  ASSERT (! refMutation. empty ());
  
  if (   refMutation. pos_real < qInt. start
      || refMutation. pos_real + refMutation. allele. size () > qInt. stop
     )
    return no_index;
    
  size_t pos = qInt. start;
  size_t frameshift_i = no_index;
	size_t i = 0;
	while (i < qseq. size ())
	{
	  ASSERT (pos <= refMutation. pos_real);
	  if (qseq [i] != '-')
	  {
      if (pos == refMutation. pos_real)
      {
        if (   qseq.    substr (i, refMutation. allele. size ()) == refMutation. allele
            && sseq. substr (i, refMutation. allele. size ()) == refMutation. allele
            && (pos == 0 || (i > 0 && qseq [i - 1] == sseq [i - 1]))
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
          ASSERT (alignStop <= qseq. size ());
          if (refStop == qlen || qseq [alignStop] == sseq [alignStop])
            return i;
        }
        return no_index;
      }
  	  pos++;
    }
    i++;
	}

	while (i < qseq. size ())
	{
	  ASSERT (frameshift_i != no_index);
	  if (qseq [i] != '-')
	  {
	  #if 0
	    PRINT (i);
	    PRINT (pos);
	    PRINT (refMutation. getStop ());
	    PRINT (refMutation. frameshift);
	    PRINT (qseq    [i]);
	    PRINT (sseq [i]);	    
	  #endif
	    ASSERT (refMutation. getStop ());
	    if (pos == refMutation. getStop () - 1 + refMutation. frameshift)
	    {
	      if (   qseq    [i] == '*'
	          && sseq [i] == '*'
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
    
  Hsp::qc ();
    
  refMutation. qc ();
  for (const SeqChange& seqChange : seqChanges)
  {
    seqChange. qc (); 
    QC_ASSERT (seqChange. al == this);
  }
}



}  // namespace


