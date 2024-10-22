// seq.cpp

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
*   Sequence utilities
*
*/


#undef NDEBUG

#include "seq.hpp"

#include "common.inc"



namespace Seq_sp
{


////////////////////////////////////////////////////////////////////////////

// Alphabets

const char* dnaAlphabet    = "acgt";
const char* extSparseDnaAlphabet = "-acgtbdhkmnrsvwyACGTBDHKMNRSVWY";
const char* extDnaAlphabet = & extSparseDnaAlphabet [1];
const char* dnaWildcards   = & extDnaAlphabet [4];


const char* peptideAlphabet              =  "ACDEFGHIKLMNPQRSTVWY";
const char* extPeptideAlphabet           =  "ACDEFGHIKLMNPQRSTVWYXBZJUO";
const char* extSparseTermPeptideAlphabet = "-ACDEFGHIKLMNPQRSTVWYXBZJUO*";
const char* terminator             = & extSparseTermPeptideAlphabet [27];
const char* extTermPeptideAlphabet = & extSparseTermPeptideAlphabet [1];
const char* peptideWildcards       = & extPeptideAlphabet [20];



size_t alphabet2Pos (const char* alphabet,
                     char        c)
{ 
  ASSERT (alphabet);
  
  if (const char* s = strchr (alphabet, c))
    return (size_t) (s - alphabet);
  cout << c << "  " << alphabet << endl;
  ERROR;
  return no_index;
}



//////////////////////////////// Seq ////////////////////////////////////

#if 0
int SparseSeq2SeqPos (const char* SparseSeq,
              		      int         SparseSeqPos)
{
  ASSERT (SparseSeq != nullptr);


  if (SparseSeqPos < 0)
    return SparseSeqPos;


  if (SparseSeqPos <= (int) strlen (SparseSeq))
    while (SparseSeq [SparseSeqPos] == '-')
      SparseSeqPos++;

  size_t N = SparseSeqPos;
  For (i, (size_t) SparseSeqPos)
    {
      if (SparseSeq [i] == '\0')
        break;
      if (SparseSeq [i] == '-')
        N--;
    }


  return N;
}



int Seq2SparseSeqPos (const char* SparseSeq,
             		       int         SeqPos)
{
  ASSERT (SparseSeq != nullptr);


  if (SeqPos < 0)
    return SeqPos;


  size_t N = SeqPos;
  For (i, N + 1)
    {
      if (SparseSeq [i] == '\0')
        break;
      if (SparseSeq [i] == '-')
        N++;
    }
  ASSERT (N > strlen (SparseSeq) || SparseSeq [N] != '-');


  return N;
}
#endif




constexpr size_t fastaLineLen = 80;




// Seq

Seq::Seq (LineInput &fasta,
          size_t reserveLen,
          bool sparse_arg,
     	    bool makeUpper)
: sparse (sparse_arg)
{
	QC_ASSERT (! fasta. eof);
	QC_ASSERT (! fasta. line. empty ());
  QC_ASSERT (fasta. line [0] == '>');
  name = fasta. line. substr (1);
  replace (name, '\t', ' ');
  qcName ();
  
  seq. reserve (reserveLen);

  while (   fasta. nextLine () 
         && ! (   strBlank (fasta. line)
               || fasta. line [0] == '>'
       	      )
       	)
  {
    trimTrailing (fasta. line);
    seq += fasta. line;
  }

  while (   ! fasta. eof        
         && strBlank (fasta. line)
        )
    fasta. nextLine ();

  if (makeUpper)
    strUpper (seq);
  else
    strLower (seq);

	if (! sparse)
	  unSparse ();
}



void Seq::qcName () const
{
	QC_ASSERT (! strBlank (name));
	QC_ASSERT (! contains (name, "\t"));
}



void Seq::qcAlphabet () const
{
  const string::const_iterator it = stringInSet (seq, getSeqAlphabet ());
  if (it != seq. end ())	
    throw runtime_error ("Bad sequence character in " + strQuote (getId ()) + ": " + to_string (int (*it)) + " " + seq. substr ((size_t) (it - seq. begin ())));
}



void Seq::printCase (ostream &os,
              	     bool makeUpper) const
{
  os << ">" << name;

  CONST_ITER (string, it, seq)
  {
    if ((size_t) (it - seq. begin ()) % fastaLineLen == 0)
    	os << endl;

    char c = *it;
    if (makeUpper)
     	c = toUpper (c);
    os << c;
  }
  os << endl;
  
  ASSERT (os. good ());
}



size_t Seq::getTaxonStart (const string &s) 
{
  size_t brackets = 0;
  FOR_REV (size_t, i, s. size ())
  {
    if (s [i] == ']')
      brackets++;
    else if (s [i] == '[')
    {
      if (brackets == 1)
        return i;
      else if (brackets == 0)
    	  return string::npos;  // Non-balanced brackets: too many '['
      else
        brackets--;
    }
  //cout << s. substr (0, i + 1) << ' ' << i << ' ' << brackets << endl;  
    if (brackets == 0)
    	return string::npos;  // No brackets
  }
 	return string::npos;  // Non-balanced brackets: too many ']'
 
#if 0
  if (s [s. size () - 1] != ']')
  	return string::npos;
  const size_t pos = s. rfind ('[');
  if (pos == string::npos)
  	return string::npos;
  ASSERT (pos < s. size () - 1);
  return pos;
#endif
}



long /*CSeq_id::TGi*/ Seq::getGi () const
{
	const char* prefix = "gi|";
	size_t start = name. find (prefix, 0);
	size_t end;
	if (start == string::npos)
	{
		start = 0;
		end = name. find (' ');
	}
	else
	{
		ASSERT (start < getIdSize ());
		start += strlen (prefix);
	  end = name. find ('|', start);
	  ASSERT (end != string::npos);
  }
  if (end == string::npos)
  	end = name. size ();
  return str2<long> (name. substr (start, end - start));
}



string Seq::getDescription (bool trimTaxon) const
{
	const size_t idSize = getIdSize ();
	if (idSize == name. size ())
		return noString;
	ASSERT (name [idSize] == ' ');
	string desc (name. substr (idSize + 1));
	if (trimTaxon)
	{
		const size_t taxonStart = getTaxonStart ();
		if (taxonStart != string::npos)
			desc. erase (taxonStart - (idSize + 1));
  }
	trim (desc);
	return desc;
}



size_t Seq::getXs () const
{
	size_t len = 0;
	for (const char c : seq)
	  if (isAmbiguous (c))
	  	len++;
 	return len;
}



size_t Seq::getContiguousXs () const
{
	size_t maxLen = 0;
	size_t len = 0;
	for (const char c : seq)
	  if (isAmbiguous (c))
	  	len++;
	  else
	  {
	  	maximize (maxLen, len);
	  	len = 0;
	  }
 	maximize (maxLen, len);
 	return maxLen;
}



map<char,size_t> Seq::getCharCount () const
{ 
  map<char,size_t> m;
  for (const char c : seq)
    m [c] ++;
  return m;
}



#if 0
size_t* Seq::GetAlphabetCount () const
{
  const size_t Len = strlen (SeqAlphabet);

  size_t* AlphabetCount = NewUintArray (Len);
  For (i, Len)
    AlphabetCount [i] = 0;

  ForString (j, seq)
    {
      bool Found = false;
      For (k, Len)
        if (seq [j] == SeqAlphabet [k])
          {
            AlphabetCount [k] ++;
            Found = true;
            break;
          }
      ASSERT (Found);
    }


  return AlphabetCount;
}



void SaveSeq (const Seq* seq,
              const char* DirName)
{
  ASSERT (GoodObject (seq));


  if (StrIsEmpty (DirName))
    seq->SaveFile (stdout, false);
  else
  {
    ASSERT (GoodDirName (DirName));
    char FName [1024];
    DirFile2PathName (FName, DirName, seq->Name);
    StrReplaceSet (FName, " ", '\0');
    seq->Save (FName, false);
  }
}



// _SEQ_COLLECTION

void _SEQ_COLLECTION::SaveFile (FILE*   F,
              		                bool    UpperCase) const
{
  ASSERT (F != nullptr);

  ForCollection (i, *this)
    GetSeq (i) -> SaveFile (F, UpperCase);
}



void _SEQ_COLLECTION::Print (bool UpperCase) const
{
  printf ("# sequences = %u\n", Count ());
  printf ("\n");

  SaveFile (stdout, UpperCase);
}



void _SEQ_COLLECTION::Save (const char* FName,
                    		 	    bool        UpperCase) const
{
  FILE* F = fopen (FName, "w");
  ASSERT (F != nullptr);

  SaveFile (F, UpperCase);

  fclose (F);
}



void _SEQ_COLLECTION::SaveDir (const char* DirName) const
{
  ForCollection (i, *this)
    SaveSeq (GetSeq (i), DirName);
}



byte* String2Qual (char* S,
                   size_t  Len)
{
  ASSERT (S != nullptr);


  if (StrIsEmpty (S))
    return nullptr;


  byte* Qual = NewByteArray (Len);
  ASSERT (Qual != nullptr);


  char* S1 = S;
  For (i, Len)
    if (StrIsBlank (S1))
      Qual [i] = 0;  // ??
    else
      {
        const char* S2 = StrChop (S1, " ");
        uint Score;
        if (sscanf (S2, "%u", & Score) != 1 ||
            Score > 99)
          {
            printf ("%s\n", S2);
            ERROR;
          }
        Qual [i] = Score;
      }


  return Qual;
}
#endif




///////////////////////////// MultiFasta /////////////////////////////////

void Multifasta::qcNewSeq () const
{
	if (   ! in. eof 
	    && ! (   ! in. line. empty () 
	          && in. line [0] == '>'
	         )
	   )
	  throw runtime_error ("Error in Multifasta, " + in. lineStr ());   
}

	  


///////////////////////////////// Dna ////////////////////////////////////

size_t nuc2num (char wildNucleotide)
{
  const char c = toLower (wildNucleotide);
  switch (c)
  {
    case 'a': return 0; break;
    case 'c': return 1; break;
    case 'g': return 2; break;
    case 't': return 3; break;
  }  
  if (isAmbigNucl (c))
    return 4;
  NEVER_CALL;
}



#if 0
void PrintNTProb (const NT_PROBABILITY Prob)
{
  For (i, 5)
    printf ("%c: %0.3f\n", "acgtb" [i], Prob [i]);
}



static void CompressNucleotide (char  WildNucleotide,
                                char* Bits,
                                uint  &BitNum)
// Input: WildNucleotide: in ExtDnaAlphabet + '-' + uppercase
// Output: Bits []
// Update: BitNum
{
  ASSERT (Bits != nullptr);


  if (CharInSet (WildNucleotide, "acgt"))
    {
      PutBit0 (Bits, BitNum);
      switch (WildNucleotide)
        {
          case 'a':
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
          Case 'c':
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
          Case 'g':
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
          Case 't':
            PutBit1 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
          Default: ERROR;
        }
    }
  else
    {
      PutBit1 (Bits, BitNum);
      switch (WildNucleotide)
        {
          case 'm':
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
          Case 'r':
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
          Case 'w':
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
          Case 's':
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
          Case 'y':
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
          Case 'k':
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
          Case 'v':
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
          Case 'h':
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
          Case 'd':
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
          Case 'b':
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
          Case 'n':
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
          Case 'A':
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
          Case 'C':
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
          Case 'G':
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
          Case 'T':
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
          Case 'M':
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
          Case 'R':
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
          Case 'W':
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
          Case 'S':
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
          Case 'Y':
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
          Case 'K':
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
          Case 'V':
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
          Case 'H':
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
          Case 'D':
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
          Case 'B':
            PutBit1 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
          Case 'N':
            PutBit1 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
          Case '-':
            PutBit1 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
            PutBit1 (Bits, BitNum);
            PutBit0 (Bits, BitNum);
          Default: ERROR;
        }
    }
}



static char UncompressNucleotide (const char* Bits,
                                  uint        &BitNum)
// Update: BitNum
// Return: Wild nucleotide in ExtDnaAlphabet + '-' + uppercase
{
  ASSERT (Bits != nullptr);


  char C;
  if (GetBit0 (Bits, BitNum))
    if (GetBit0 (Bits, BitNum))
      if (GetBit0 (Bits, BitNum))
        C = 'a';
      else
        C = 'c';
    else
      if (GetBit0 (Bits, BitNum))
        C = 'g';
      else
        C = 't';
  else
    if (GetBit0 (Bits, BitNum))
      if (GetBit0 (Bits, BitNum))
        if (GetBit0 (Bits, BitNum))
          if (GetBit0 (Bits, BitNum))
            if (GetBit0 (Bits, BitNum))
              C = 'm';
            else
              C = 'r';
          else
            if (GetBit0 (Bits, BitNum))
              C = 'w';
            else
              C = 's';
        else
          if (GetBit0 (Bits, BitNum))
            if (GetBit0 (Bits, BitNum))
              C = 'y';
            else
              C = 'k';
          else
            if (GetBit0 (Bits, BitNum))
              C = 'v';
            else
              C = 'h';
      else
        if (GetBit0 (Bits, BitNum))
          if (GetBit0 (Bits, BitNum))
            if (GetBit0 (Bits, BitNum))
              C = 'd';
            else
              C = 'b';
          else
            if (GetBit0 (Bits, BitNum))
              C = 'n';
            else
              C = 'A';
        else
          if (GetBit0 (Bits, BitNum))
            if (GetBit0 (Bits, BitNum))
              C = 'C';
            else
              C = 'G';
          else
            if (GetBit0 (Bits, BitNum))
              C = 'T';
            else
              C = 'M';
    else
      if (GetBit0 (Bits, BitNum))
        if (GetBit0 (Bits, BitNum))
          if (GetBit0 (Bits, BitNum))
            if (GetBit0 (Bits, BitNum))
              C = 'R';
            else
              C = 'W';
          else
            if (GetBit0 (Bits, BitNum))
              C = 'S';
            else
              C = 'Y';
        else
          if (GetBit0 (Bits, BitNum))
            if (GetBit0 (Bits, BitNum))
              C = 'K';
            else
              C = 'V';
          else
            if (GetBit0 (Bits, BitNum))
              C = 'H';
            else
              C = 'D';
      else
        if (GetBit0 (Bits, BitNum))
          if (GetBit0 (Bits, BitNum))
            if (GetBit0 (Bits, BitNum))
              C = 'B';
            else
              C = 'N';
          else
            if (GetBit0 (Bits, BitNum))
              C = '-';
            else
              ERROR
        else
          ERROR;


  return C;
}



void CompressDna (const char* seq,
                  char*       CompressedSeq,
                  size_t        &CompressedSeqLen)
{
  ASSERT (seq != nullptr);
  ASSERT (CompressedSeq != nullptr);


  // CompressedSeq []
  uint BitNum = 0;
  ForString (i, seq)
    CompressNucleotide (seq [i], CompressedSeq, BitNum);

  const uint RestBits = 8 - (BitNum % 8);
  if (RestBits < 8)
    For (i, RestBits)
      PutBit0 (CompressedSeq, BitNum);


  // CompressedSeqLen
  ASSERT (BitNum % 8 == 0);
  CompressedSeqLen = BitNum / 8;
}



void UncompressDna (size_t        SeqLen,
                    const char* CompressedSeq,
                    char*       &seq,
                    size_t        &CompressedSeqLen)
{
  ASSERT (CompressedSeq != nullptr);


  // seq
  seq = (char*) malloc (SeqLen + 1);
  ASSERT (seq != nullptr);

  uint BitNum = 0;
  For (i, SeqLen)
    seq [i] = UncompressNucleotide (CompressedSeq, BitNum);
  seq [SeqLen] = '\0';

  // CompressedSeqLen
  CompressedSeqLen = Bit2ByteNum (BitNum);
}
#endif



uchar wild2nucleotides (char wildNucleotide,
              		      bool acgtb [5])
{
  bool& a = acgtb [0];
  bool& c = acgtb [1];
  bool& g = acgtb [2];
  bool& t = acgtb [3];
  a = false;
  c = false;
  g = false;
  t = false;


  if (wildNucleotide == ' ')
  {
    acgtb [4] = false;
    return 0;
  }


  if (wildNucleotide == '-')
  {
    acgtb [4] = true;
    return 1;
  }


  acgtb [4] = isUpper (wildNucleotide);
  switch (toLower (wildNucleotide))
    {
      case 'a': a = true;
      Case 'c':           c = true;
      Case 'g':                     g = true;
      Case 't':                               t = true;
      Case 'm': a = true; c = true;
      Case 'r': a = true;           g = true;
      Case 'w': a = true;                     t = true;
      Case 's':           c = true; g = true;
      Case 'y':           c = true;           t = true;
      Case 'k':                     g = true; t = true;
      Case 'v': a = true; c = true; g = true;
      Case 'h': a = true; c = true;           t = true;
      Case 'd': a = true;           g = true; t = true;
      Case 'b':           c = true; g = true; t = true;
      Case 'n': a = true; c = true; g = true; t = true;
      Default: printf ("!%c! %d\n", wildNucleotide, wildNucleotide);
               ERROR;
    }


  uchar n = 0;
  FOR (uchar, i, 5)
    if (acgtb [i])
      n++;
  ASSERT (n > 0);


  return n;
}



#if 0
void Wild2NucleotideFreq (char           WildNucleotide,
                          NT_PROBABILITY acgtb)
{
  bool ACGTB_P [5];
  const uint N = Wild2Nucleotides (WildNucleotide, ACGTB_P);

  const PROBABILITY P = 1 / (Float) N;
  For (i, 5)
    acgtb [i] = ACGTB_P [i] ? P : 0;
}
#endif



char nucleotides2wild (const bool acgtb [5])
{
  const bool& a = acgtb [0];
  const bool& c = acgtb [1];
  const bool& g = acgtb [2];
  const bool& t = acgtb [3];

  char res = '\0';
  if (a)
    if (c)
      if (g)
        if (t)
          res = 'n';
        else
          res = 'v';
      else
        if (t)
          res = 'h';
        else
          res = 'm';
    else
      if (g)
        if (t)
          res = 'd';
        else
          res = 'r';
      else
        if (t)
          res = 'w';
        else
          res = 'a';
  else
    if (c)
      if (g)
        if (t)
          res = 'b';
        else
          res = 's';
      else
        if (t)
          res = 'y';
        else
          res = 'c';
    else
      if (g)
        if (t)
          res = 'k';
        else
          res = 'g';
      else
        if (t)
          res = 't';
        else
          res = '-';


  if (acgtb [4])
    return toUpper (res);
  else
    if (res == '-')
      return ' ';
    else
      return res;
}



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
    	throw runtime_error ("Bad wild nucleotide " + to_string (wildNucleotide));
  }
  if (isupper (wildNucleotide))
    r = toUpper (r);
  ASSERT (charInSet (r, extDnaAlphabet));

  return r;
}



char getUnionNucleotide (const string& charSet)
{
  if (charSet. empty ())
    return ' ';


  bool acgtb [5] = {false, false, false, false, false};
  for (const char c : charSet)
  {
    bool _acgtb [5];
    wild2nucleotides (c, _acgtb);
    FOR (size_t, j, 5)
      maximize (acgtb [j], _acgtb [j]);
  }


  return nucleotides2wild (acgtb);
}



char getIntersectNucleotide (const string& charSet)
{
  if (charSet. empty ())
    return ' ';


  bool acgtb [5] = {true, true, true, true, true};
  for (const char c : charSet)
  {
    bool _acgtb [5];
    wild2nucleotides (c, _acgtb);
    FOR (size_t, j, 5)
      minimize (acgtb [j], _acgtb [j]);
  }


  return nucleotides2wild (acgtb);
}



bool nucleotideMatch (char wildNucleotide1,
               			  char wildNucleotide2)
{
  ASSERT (wildNucleotide1 != '\0');
  ASSERT (wildNucleotide2 != '\0');


  if (wildNucleotide1 == wildNucleotide2)
    return true;

  string s;
  s. resize (2);
  s [0] = wildNucleotide1;
  s [1] = wildNucleotide2;
  return (getIntersectNucleotide (s) != ' ');
}



#if 0
bool NucleotideSeqMatch (const char* Seq1,
                         const char* Seq2)
{
  ASSERT (Seq1 != nullptr);
  ASSERT (Seq2 != nullptr);


  ForString (i, Seq1)
    if (Seq2 [i] == '\0')
      return true;
    else if (! NucleotideMatch (Seq1 [i], Seq2 [i]))
      return false;


  return true;
}



bool MoreGeneralNucleotide (char WildNucleotide1,
                  			       char WildNucleotide2)
{
  if (WildNucleotide1 == WildNucleotide2)
    return true;


  bool ACGTB1 [5], ACGTB2 [5];
  Wild2Nucleotides (WildNucleotide1, ACGTB1);
  Wild2Nucleotides (WildNucleotide2, ACGTB2);
  For (i, 5)
    if (! ACGTB1 [i] && ACGTB2 [i])
      return false;


  return true;
}



void SelectSpecificNucleotides (char* CharSet)
{
  ASSERT (CharSet != nullptr);


  ForDown (i, strlen (CharSet))
    ForString (j, CharSet)
      if (i != (int) j &&
          MoreGeneralNucleotide (CharSet [i], CharSet [j]))
        {
          StrDelete (& CharSet [i], 1);
          break;
        }
}



size_t CountAmbiguousNucleotides (const char* seq)
{
  size_t N = 0;
  ForString (i, seq)
    if (seq [i] != '-' &&
        ! CharInSet (seq [i], DnaAlphabet))
      N++;


  return N;
}
}
#endif



string& reverseDna (string &seq)
{
  const size_t len = seq. size ();
  if (! len)
  	return seq;

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
 
  return seq;
}



char codon2aa (const char codon [3],
               Gencode gencode,
               bool lowercasePossibleStartCodon)
{
  char c [3];
  FOR (size_t, i, 3)
    c [i] = toLower (codon [i]);


  char aa = 'X';
  switch (c [0])
  {
    case 't':
     	switch (c [1])
   	  {
   	    case 't':
   	      switch (c [2])
      		{
      		  case 't': case 'c': case 'y': aa = 'F';
      		  Case 'a': case 'g': case 'r': aa = 'L';
      		}
   	    Case 'c':                         aa = 'S'; 
   	    Case 'a':
   	      switch (c [2])
      		{
      		  case 't': case 'c': case 'y': aa = 'Y';
      		  Case 'a': case 'g': case 'r': aa = *terminator;
      		}
   	    Case 'g':
   	      switch (c [2])
      		{
      		  case 't': case 'c': case 'y': aa = 'C';
      		  Case 'a': 
      		    switch (gencode)
      		    {
      		      case 4: case 25:          aa = 'W';
      		      Default:                  aa = *terminator;
      		    }
      		  Case 'g':                     aa = 'W';
      		}
   	    Case 'r':
   	      if (c [2] == 'a')          aa = *terminator;
   	  }
    Case 'c':
     	switch (c [1])
   	  {
   	    case 't':                         aa = 'L';
   	    Case 'c':                         aa = 'P';
   	    Case 'a':
   	      switch (c [2])
      		{
      		  case 't': case 'c': case 'y': aa = 'H';
      		  Case 'a': case 'g': case 'r': aa = 'Q';
      		}
   	    Case 'g':                         aa = 'R';
   	  }
    Case 'a':
     	switch (c [1])
   	  {
   	    case 't':
   	      switch (c [2])
      		{
      		  case 't': case 'c': case 'a':
      		  case 'y': case 'w': case 'm':
      		  case 'h':                     aa = 'I';
      		  Case 'g':                     aa = 'M'; 
      		}
   	    Case 'c':                         aa = 'T';
   	    Case 'a':
   	      switch (c [2])
      		{
      		  case 't': case 'c': case 'y': aa = 'N';
      		  Case 'a': case 'g': case 'r': aa = 'K';
      		}
   	    Case 'g':
   	      switch (c [2])
      		{
      		  case 't': case 'c': case 'y': aa = 'S';
      		  Case 'a': case 'g': case 'r': aa = 'R';
      		}
   	  }
    Case 'g':
     	switch (c [1])
   	  {
   	    case 't':                         aa = 'V';
   	    Case 'c':                         aa = 'A';
   	    Case 'a':
   	      switch (c [2])
      		{
      		  case 't': case 'c': case 'y': aa = 'D';
      		  Case 'a': case 'g': case 'r': aa = 'E';
      		}
   	    Case 'g':                         aa = 'G';
   	  }
    Case 'm':
      if (c [1] == 'g' &&
          charInSet (c [2], "agr"))  aa = 'R';
    Case 'r':
      if (c [1] == 'a' &&
          charInSet (c [2], "tcy"))  aa = 'B';
    Case 's':
      if (c [1] == 'a' &&
          charInSet (c [2], "agr"))  aa = 'Z';
    Case 'y':
      if (c [1] == 't' &&
          charInSet (c [2], "agr"))  aa = 'L';
  }
  ASSERT (charInSet (aa, extTermPeptideAlphabet));
  
  
  if (lowercasePossibleStartCodon)
  {
    bool startcodon = false;
    switch (gencode)
    {
      case  1: if (    c[0] == 'a' && c[1] == 't' && c[2] == 'g' )  startcodon = true;
      Case  4: if (   (               c[1] == 't' && c[2] == 'g')
                   || (c[0] == 'a' && c[1] == 't'               )
                   || (c[0] == 't' && c[1] == 't' && c[2] == 'a'))  startcodon = true;
      Case 11: if (                  (c[1] == 't' && c[2] == 'g')
                   || (c[0] == 'a' && c[1] == 't')               )  startcodon = true;  
      Case 25: if (    c[0] != 'c' && c[1] == 't' && c[2] == 'g' )  startcodon = true;
      Default: throw runtime_error (FUNC "Genetic code " + to_string (gencode) + " is not implemented");
    }
    if (startcodon)
      aa = toLower (aa);
  }


  return aa;
}



#if 0
bool MayBeTerminatorCodon (const char Codon [3])
{
  return
    CodonMatch (Codon, "taa") ||
    CodonMatch (Codon, "tag") ||
    CodonMatch (Codon, "tga");
}



size_t Dna2PeptidePos (size_t  DnaPos,
              		   Frame frame,
              		   size_t  DnaLen)
{
  ASSERT (isFrame (frame));
  ASSERT (DnaPos < DnaLen);


  if (frame < 0)
    DnaPos = (DnaLen - 1) - DnaPos;


  return (DnaPos - (AbsInt (frame) - 1)) / 3;
}



size_t Peptide2DnaPos (size_t  PeptidePos,
              		   Frame frame,
              		   size_t  DnaLen)
{
  ASSERT (isFrame (frame));


  const int DnaPos = (AbsInt (frame) - 1) + PeptidePos * 3;
  ASSERT (DnaPos < (int) DnaLen);

  if (frame > 0)
    return DnaPos;
  else
    return (DnaLen - 1) - DnaPos;
}




// phred scores

PROBABILITY PhredScore2Prob (PHRED_SCORE Score)
{
  ASSERT (Score <= AccuratePhredScore);


  static PROBABILITY Res [AccuratePhredScore + 1];
  if (Res [Score] == 0)
    Res [Score] = MinimumFloat (pow (10.0, - (Float) Score / 10), 0.8); 


  return  Res [Score];
}



PHRED_SCORE Prob2PhredScore (PROBABILITY P)
{
  ASSERT (IsProbability (P));


  if (NullFloat (P))
    return MaxAutoPhredScore;
  else
    return MinimumInt (Round (10.0 * log10 (1.0 / P)), MaxAutoPhredScore);
}



const PROBABILITY MinNonAmbigProb = 0.78;  // 0.6 before 2/18/04 
// See Contig::GetConsensusNucleotide()
const PROBABILITY AmbigThreshold = (1.0 - MinNonAmbigProb) / 4.0 - 1e-3;



PROBABILITY Nucleotide_MonoErrorProb2WildNucleotide (uint        ACGTBNum,
                                                     PROBABILITY ErrorProb,
                                                     char        WildNucleotide)
// Return: Probability of WildNucleotide given ErrorProb, |WildNucleotide| and correct nucleotide ACGTBNum (including blank)
// Model: Only wild nucleotides of size K = |WildNucleotide| are considered, K \in [1,5]
//        P1 := Probability (wild nucleotide     containing the nucleotide)
//        P2 := Probability (wild nucleotide not containing the nucleotide)
//        0 <= P1, P2 <= 1
//        P1 * choose(4,K-1)           + P2 * choose(4,K) = 1
//        P1 * choose(4,K-1) * (K-1)/K + P2 * choose(4,K) = ErrorProb
//        CorrProb := 1 - ErrorProb
//        CorrProb = P1 * choose(4,K-1) / K
//        P1 = CorrProb * K / choose(4,K-1) <= 1 / choose(4,K-1)
//        CorrProb * K <= 1
//        Quality score --> raw ErrorProb
//        For Contig:GetConsensusNucleotide():
//          P2 / (K * P1 + (5-K) * P2) <= AmbigThreshold - 1e-3
{
  ASSERT (ACGTBNum <= 5);
  ASSERT (IsProbability (ErrorProb));
//ASSERT (LessEqualFloat (ErrorProb, 4.0/5.0));  Old 
  ASSERT (1.0 / (AmbigThreshold - 1e-3) - 5 > 0);


  bool acgtb [5];
  const uint K = Wild2Nucleotides (WildNucleotide, acgtb);
  ASSERT (BetweenEqual (K, 1, 5));

  const PROBABILITY CorrProb = (K == 1 ? 1.0 - ErrorProb : 0.0);
  const int C  = Choose (4, K);
  const int C1 = Choose (4, K - 1);
  const PROBABILITY P1 = 
    MinimumFloat (MaximumFloat (CorrProb * (double) K / C1, 
                                1.0 / (K * C / (1.0 / (AmbigThreshold - 1e-3) - 5 + K) + C1)
                               ),
                  1.0 / C1);
  
  if (acgtb [ACGTBNum])  // WildNucleotide is correct 
    return P1;
  else                   // WildNucleotide is incorrect
    return (1.0 - P1 * C1) / C;

/* Old
  if (K == 5)
    ErrorProb = 1.0 - 1.0 / 5.0;
  else
    maximize (ErrorProb, 1.0 - 1.0 / K);
  
  if (acgtb [ACGTBNum])  // Correct 
    return (1.0 - ErrorProb) / Choose (5 - 1, (int) K - 1);
  else
    return ErrorProb / Choose (5 - 1, K);
*/
}
#endif



// Dna

//const double Dna::stdMinComplexity = 2.0;  // PAR



Dna::Dna (Multifasta &fasta,
          size_t reserveLen,
          bool sparse_arg)
: Seq (fasta. in, reserveLen, sparse_arg, false) 
{ 
	ASSERT (! fasta. aa);
	fasta. prog (getId ());
}



double Dna::getComplexityInt (size_t start,
                              size_t end) const
{
  if (seq. empty ())
    return 0;


  // DinucFreq[][], n
  size_t dinucFreq [4] [4];
  FOR (size_t, i, 4)
    FOR (size_t, j, 4)
      dinucFreq [i] [j] = 0;
  size_t n = 0;
  minimize (end, seq. size ());
  FOR_START (size_t, i, start + 1, end)
  {  	
    const size_t a = nuc2num (seq [i - 1]);
    const size_t b = nuc2num (seq [i]);
    if (   a < 4 
    	  && b < 4
    	 )
    {
      dinucFreq [a] [b] ++;
      n++;
    }
  }

  // Entropy
  double s = 0;
  FOR (size_t, i, 4)
    FOR (size_t, j, 4)
      if (dinucFreq [i] [j])
      {
       	const double p = (double) dinucFreq [i] [j] / (double) n;
       	s -= p * log (p);
      }
  ASSERT (s >= 0);

  return s;
}



#if 0
Dna::Dna (const Dna &dna)
: Seq (dna)
, Qual      (nullptr)
{
  // Qual
  if (dna. Qual != nullptr)
    {
      CreateQual ();
      ForString (i, seq)
        Qual [i] = dna. Qual [i];
    }
}



void Dna::PrintHTML (bool        UpperCase,
                     PHRED_SCORE MinGoodQual) const
{
  printf (">%s", Name);
  bool GoodQual = true;
  ForString (i, seq)
    {
      if (i % FastaLineLen == 0)
       	printf ("\n");

      char C = seq [i];
      if (UpperCase)
       	C = ToUpper (C);
       	
      if (Qual != nullptr) 
        if (Qual [i] < MinGoodQual)
          {
            if (GoodQual)
              printf ("<font color=gray>");
            GoodQual = false;
          }
        else
          {
            if (! GoodQual)
              printf ("</font>");
            GoodQual = true;
          }
       	
       printf ("%c", C);
     }
  if (! GoodQual)
    printf ("</font>");
  printf ("\n");
}
  


SEQ_ANNOT_LIST* Dna::GetAnnotList (int GoodStart) const
{
  SEQ_ANNOT_LIST* AnnotList = new SEQ_ANNOT_LIST (strlen (seq), GoodStart);
  ASSERT (GoodObject (AnnotList));

  SEQ_ANNOT* SeqAnnot = AnnotList->AddAnnot ("seq");
  ASSERT (GoodObject (SeqAnnot));
  SeqAnnot->AddCharColor (DnaWildcards, ColorBad);
  SeqAnnot->PutStr (0, seq);


  // Qual
  if (Qual != nullptr)
    {
      SEQ_ANNOT* QualAnnot1 = AnnotList->AddAnnot ("Score 10");
      ASSERT (GoodObject (QualAnnot1));
      QualAnnot1->AddCharColor ("01", ColorBad);
      SEQ_ANNOT* QualAnnot2 = AnnotList->AddAnnot ("Score  1");
      ASSERT (GoodObject (QualAnnot2));

      ForString (i, seq)
        {
          QualAnnot1->PutNumber (i, Qual [i] / 10);
          QualAnnot2->PutNumber (i, Qual [i] % 10);
        }
    }


  return AnnotList;
}



// Qual

void Dna::CreateQual ()
{
  ASSERT (Qual == nullptr);


  Qual = NewByteArray (strlen (seq));
  ASSERT (Qual != nullptr);
}



bool Dna::GoodQual () const
{
  ASSERT (Qual != nullptr);


  ForString (i, seq)
    {
      if (Qual [i] > 99)
        return false;
      if (seq [i] == 'n' && Qual [i] != 0)
      {
      //printf ("%u %c %d\n", i + 1, seq [i], Qual [i]);  
        return false;
      }
    }


  return true;
}



void Dna::CopyQual (const PHRED_SCORE* SourceQual)
{
  ASSERT (SourceQual != nullptr);


  CreateQual ();
  ForString (i, seq)
    Qual [i] = SourceQual [i];
  if (! GoodQual ())
    {          
      printf ("%s\n", Name);
      ERROR;
    }
}



void Dna::ReadQual (const char* FName)
{
  ASSERT (Qual == nullptr);
  ASSERT (! StrIsEmpty (FName));


  // F
  FILE* F = fopen (FName, "r");
  ASSERT (F != nullptr);

  char C;
  if (fscanf (F, "%c", & C) != 1)
    ERROR;
  if (C != '>')
    ERROR
  SkipLine (F);


  // Qual
  CreateQual ();

  uint Score;
  ForString (i, seq)
    if (fscanf (F, "%u", & Score) == 1)
      Qual [i] = Score;
    else
      ERROR;

  ASSERT (GoodQual ());


  fclose (F);
}



void Dna::SetQual (PHRED_SCORE DefaultScore)
{
  ASSERT (Qual != nullptr);


  ForString (i, seq)
    Qual [i] = CharInSet (seq [i], "Nn") ? 0 : DefaultScore;
  ASSERT (GoodQual ());
}



PHRED_SCORE Dna::GetMinQual () const
{
  ASSERT (Qual != nullptr);


  PHRED_SCORE MinQual = 255;
  ForString (i, seq)
    minimize (MinQual, Qual [i]);

  return MinQual;
}



PHRED_SCORE Dna::GetMaxQual () const
{
  ASSERT (Qual != nullptr);


  PHRED_SCORE MaxQual = 0;
  ForString (i, seq)
    maximize (MaxQual, Qual [i]);

  return MaxQual;
}



void Dna::Qual2MeanVar (size_t  Start,
                        size_t  End,
                        Float &Mean,
                        Float &Var) const
{
  ASSERT (Qual != nullptr);
  ASSERT (Start <= End);
  
  
  const size_t N = End - Start;

  Float S = 0.0;
  ForStart (i, Start, End)
    S += Qual [i];
  Mean = S / N;

  S = 0.0;  
  if (N == 0)
    Var = NaN;
  else
    {
      ForStart (i, Start, End)
        S += Sqr (Qual [i] - Mean);
      Var = S / (N - 1);
    }
}



void Dna::Qual2N (PHRED_SCORE MinQual,
                  bool        GapCoded)
{
  ASSERT (Qual != nullptr);


  if (MinQual > 0)
    {
      const char NChar = GapCoded ? 'N' : 'n';
      ForString (i, seq)
        if (Qual [i] < MinQual)
          seq [i] = NChar;
    }

  delete [] Qual;
  Qual = nullptr;
}



void Dna::MinScore2Pos (PHRED_SCORE MinScore,
                        size_t        &Start,
                        size_t        &End) const
{
  ASSERT (Qual != nullptr);
  
  
  Start = 0;
  End   = 0;

  
  int L           = 0;  // Max. good quality segment length
  int L_Last      = 0;  // Max. good quality segment length ending at i
  size_t Start_Last = 0;  // Start of L_Last segment
  ForString (i, seq)
    {
      if (Qual [i] >= MinScore)
        L_Last++;
      else
        {
          L_Last = 0;
          Start_Last = i + 1;
        }
      if (maximize (L, L_Last))
        Start = Start_Last;
    }

    
  End = Start + L;
}



void Dna::QualSaveFile (FILE*       F,
                     		 PHRED_SCORE DefaultScore) const
{
  ASSERT (F != nullptr);


  fprintf (F, ">%s", Name);

  ForString (i, seq)
    {
      if (i % FastaLineLen == 0)
       	fprintf (F, "\n");

      byte B = DefaultScore;
      if (Qual != nullptr)
        B = Qual [i];

      fprintf (F, "%u ", B);
     }
  fprintf (F, "\n");
}



char* Dna::Qual2String () const
{
  ASSERT (! StrIsEmpty (seq));
  ASSERT (Qual != nullptr);
  
 
  const int SLen = strlen (seq) * 3;  // 2-digit PHRED scores
  char* S = NewCharArray (SLen + 1);
  int Pos = 0;
  ForString (i, seq)
    {
      char Score [16] = " ";
      sprintf (& Score [Pos == 0 ? 0 : 1], "%d", Qual [i]);
      ASSERT (strlen (Score) <= 3);
      strcpy (& S [Pos], Score);
      Pos += strlen (Score);
    } 
  S [Pos] = '\0';
  ASSERT (S [0] != ' ');
  
  
  return S;
 
}



void Dna::DeleteStart (size_t NewStart)
{
  const size_t OldLen = strlen (seq);
  if (NewStart > OldLen)
    {
      printf ("%s %s\n", Name, seq);
      printf ("NewStart = %u, OldLen = %u\n", NewStart, OldLen);
      ERROR;
    }

  StrDelete (seq, NewStart);
  if (Qual != nullptr)
    memmove (Qual, & Qual [NewStart], OldLen - NewStart);
}
#endif



Dna* Dna::makeComplementary () const
{
  Dna* dna = new Dna (getId () + ".rev", seq. size (), sparse);
  FFOR (size_t, i, seq. size ())
    dna->seq [i] = complementaryNucleotide (seq [seq. size () - 1 - i]);

#if 0
  if (Qual)
  {
    dna->CreateQual ();
    For (i, SeqLen)
      dna->Qual [i] = Qual [SeqLen - 1 - i];
  }
#endif

  return dna;
}



void Dna::reverse ()
{
  unique_ptr<Dna> dnaRev (makeComplementary ());
  ASSERT (seq. size () == dnaRev->seq. size ());
  seq = dnaRev->seq;

#if 0
  if (Qual)
    ForString (i, seq)
      Qual [i] = dnaRev->Qual [i];
#endif
}



Peptide Dna::makePeptide (Frame frame,
                          Gencode gencode,
                          bool lowercasePossibleStartCodon,
                          bool firstStartCodon2M,
                          size_t &translationStart) const
{
  ASSERT (isFrame (frame));
  IMPLY (firstStartCodon2M, lowercasePossibleStartCodon);

  const string peptideName (getId () + ".fr" + to_string ((int) frame));
  const size_t frameOffset = (size_t) abs (frame) - 1;
  ASSERT (frameOffset <= 2);
  const size_t aaSeqLen = seq. size () >= 3 ? (seq. size () - frameOffset) / 3 : 0;

  Peptide peptide (peptideName, aaSeqLen, false);
  peptide. pseudo = true;

  // dna
  unique_ptr<Dna> compDna;
  if (frame < 0)
    compDna. reset (makeComplementary ());
  const Dna* dna = (frame < 0) ? compDna. get () : this;
  ASSERT (dna->seq. size () == seq. size ());
  
  translationStart = frame > 0 ? frameOffset : dna->seq. size () - frameOffset;

  const char* s = dna->seq. c_str () + frameOffset;
  FOR (size_t, j, aaSeqLen)
  {
    peptide. seq [j] = codon2aa (s, gencode, lowercasePossibleStartCodon);
    s += 3;
  }
  ASSERT (peptide. seq [aaSeqLen] == '\0');
  
  if (! peptide. seq. empty () && firstStartCodon2M && isLower (peptide. seq [0]))
  	peptide. seq [0] = 'm';

  return peptide;
}



Peptide Dna::cds2prot (Gencode gencode,
                       bool trunc5,
                       bool trunc3,
	                     bool hasStopCodon,
	                     bool allowExtraStopCodon) const
{
  size_t translationStart = 0;
  Peptide pep (makePeptide (1 /*frame*/, gencode, true, true, translationStart));
  
  EXEC_ASSERT (trimSuffix (pep. name, ".fr1"));
  
  if (pep. seq. empty ())
    throw runtime_error (FUNC "empty sequence");
  if (! trunc5 && ! trunc3 && pep. seq. size () * 3 != seq. size ())
    throw runtime_error (FUNC "incomplete codon");  
  if (! trunc5 && pep. seq [0] != 'm')
    throw runtime_error (FUNC "no start codon");
    
  const size_t starPos = pep. seq. find ("*");
  if (hasStopCodon)
  {
  	if (starPos == string::npos)
  	{
  	  if (! trunc3)
        throw runtime_error (FUNC "no stop codon");
    }
    else
    {
  	  if (starPos == pep. seq. size () - 1)
  	  {
    	  EXEC_ASSERT (trimSuffix (pep. seq, "*"));
    	}
    	else
    	{
        if (allowExtraStopCodon)
          pep. seq. erase (starPos);
        else
          throw runtime_error (FUNC "stop codon at peptide position " + to_string (starPos + 1) + "\n" + seq. substr (starPos * 3));
      }
  	}
  }
  else
	  if (starPos != string::npos)
	    throw runtime_error (FUNC "has a stop codon");
  
  return pep;
}



Vector<Peptide> Dna::getOrfs (Frame frame,
                              Gencode gencode,
                              size_t len_min) const
{
  ASSERT (isFrame (frame));
  ASSERT (len_min);
  
  Vector<Peptide> peps;  peps. reserve (seq. size () / 300 + 1);  // PAR
  
  string dnaSeq (seq);
  if (frame < 0)
    reverseDna (dnaSeq);
    
  size_t start = (size_t) abs (frame) - 1;
  size_t stop = 0;
  string pepSeq;
  size_t len = 0;
  for (size_t i = start; i + 2 < dnaSeq. size (); i += 3)
  {
    const char aa = codon2aa (& dnaSeq [i], gencode, false);
    stop = i;
    if (aa == *terminator)
    {
      if (len >= len_min)
      {
        if (frame < 0)
        {
          const size_t start_ = start;
          start = seq. size () - stop;
          stop  = seq. size () - start_;
        }
        peps << std::move (Peptide (getId () + ":" + to_string (start + 1) + ".." + to_string (stop), pepSeq, false));
      }
      start = i + 3;
      pepSeq. clear ();
      len = 0;
    }
    else
    {
      pepSeq += aa;
      if (aa != 'X')
        len ++;
    }
  }
  if (len >= len_min)
  {
    if (frame < 0)
    {
      const size_t start_ = start;
      start = seq. size () - stop;
      stop  = seq. size () - start_;
    }
    peps << std::move (Peptide (getId () + ":" + to_string (start + 1) + ".." + to_string (stop), pepSeq, false));
  }
  
  return peps;
}



#if 0
bool Dna::ExistsPeptide () const
{
  bool Found = false;
  for (Frame frame = 1; frame <= 3; frame++)
    {
      Peptide* peptide = MakePeptide (frame);
      ASSERT (GoodObject (peptide));

      if (! StrIsEmpty (peptide->seq))
      {
        peptide->seq [strlen (peptide->seq) - 1] = '\0';
        if (strchr (peptide->seq, *Terminator) == nullptr)
          Found = true;
      }

      delete peptide;
    }


  return Found;
}



bool Dna::LongestCompleteCDS (size_t MinProteinLen,
                              size_t &Start,
                              size_t &End) const
{
  Start = 0;
  End   = 0;


  // OFList
  COLLECTION PeptideList;
  Orfs OFList;
  for (Frame frame = 1; frame <= 3; frame++)
  {
    Peptide* peptide = MakePeptide (frame);
    PeptideList. Append (peptide);

    peptide->AppendOpenFrameList (OFList, true, frame);
  }
  if (OFList. Empty ())
    return false;


  OFList. DeleteRepeated ();
  OFList. Sort ();
  ForCollection (i, OFList)
    {
      const Orf* OF = OFList. Get (i);
     	if (OF->Length () < MinProteinLen)
     	  return false;

     	if (! (OF->StartM && OF->EndTerminator))
     	  continue;

     	Start = OF->frame - 1 + 3 * OF->Start;
     	End   = OF->frame - 1 + 3 * (OF->End + 1);
     	ASSERT (Start < End);
     	ASSERT ((End - Start) % 3 == 0);
     	ASSERT (End <= strlen (seq));

     	return true;
    }


  return false;
}



void Dna::RefineSparse (char* Target,
                        byte* TargetQual) const
{
  ASSERT (Target != nullptr);


  size_t j = 0;
  size_t k;
  for (k = 0; Target [k] != '\0'; k++)  // 'ForString' is not used because k is used outside the loop
    if (Target [k] == '-')
      {
        ASSERT (j > 0);
        if (TargetQual != nullptr)
          TargetQual [k] = 0;
      }
    else
      {
        if (! islower (Target [k]))
          {          
            printf ("%s %u %s\n", Name, k + 1, & Target [k]);
            ERROR;
          }
        if (seq [j] == '\0' ||
            Target [k] != 'n' && ! MoreGeneralNucleotide (seq [j], Target [k])
           )
          { 
            printf ("\n");
            printf ("%s %u %u\n", Name, k, j);
            printf ("%s\n", & Target [k]);
            printf ("%s\n", & seq [j]);
            ERROR;
          }
        Target [k] = seq [j];
        if (TargetQual != nullptr)
          if (Qual == nullptr)
            TargetQual [k] = 0;
          else
            TargetQual [k] = Qual [j];
        j++;
      }
  if (seq [j] != '\0')
    { 
      printf ("\n");
      printf ("%s %u %u\n", Name, k, j);
      printf ("%s\n", & seq [j]);
      ERROR;
    }
}



bool Dna::ContainsAmbiguity () const
{
  ForString (i, DnaWildcards)
    if (strchr (seq, DnaWildcards [i]) != nullptr)
      return true;

  return false;
}



size_t Dna::GetAmbiguousPrefixEnd () const
{
  size_t i = 0;
  while (CharInSet (seq [i], DnaWildcards))
    i++;


  return i;
}



size_t Dna::GetAmbiguousSuffixStart () const
{
  int i = strlen (seq) - 1;
  while (i >= 0 &&
       	 CharInSet (seq [i], DnaWildcards))
    i--;
  i++;


  return i;
}



bool Dna::DeleteAmbiguousPrefix ()
{ 
  const size_t AmbiguousPrefixEnd = GetAmbiguousPrefixEnd ();
  DeleteStart (AmbiguousPrefixEnd); 
  
  return AmbiguousPrefixEnd > 0;
}



bool Dna::DeleteAmbiguousSuffix ()
{ 
  const size_t AmbiguousSuffixStart = GetAmbiguousSuffixStart ();
  const bool Trimmed = (seq [AmbiguousSuffixStart] != '\0');
  seq [AmbiguousSuffixStart] = '\0'; 
  
  return Trimmed;
}



void Dna::TrimN (bool GapCoded)
{
  ForDown (i, strlen (seq))
    if (seq [i] == 'N' ||
        ! GapCoded && seq [i] == 'n')
      seq [i] = '\0';
    else
      break;

  size_t k = 0;
  ForString (i, seq)
    if (seq [i] == 'N' ||
        ! GapCoded && seq [i] == 'n')
      k++;
    else
      break;
  DeleteStart (k);
}



size_t Dna::TrimBadStart (size_t        WindowLen,
              		        PROBABILITY MaxAmbigFraction) const
{
  ASSERT (IsProbability (MaxAmbigFraction));


  size_t Start = 0;
  const size_t Len = strlen (seq);
  ForEver
    {
      size_t WindowN = 0;
      bool Truncated = false;
      for (int i = Start; i < MinimumInt (Len, Start + WindowLen); i++)
        if (CharInSet (seq [i], DnaWildcards))
          {
            WindowN++;
            if (GreaterEqualFloat ((Float) WindowN / (Float) (i - Start + 1), MaxAmbigFraction))
              {
                Start = i + 1;
                Truncated = true;
                break;
              }
          }

      if (! Truncated)
        break;
    }
    

  return Start;
}



size_t Dna::TrimBadEnd (size_t        WindowLen,
              		      PROBABILITY MaxAmbigFraction) const
{
  ASSERT (IsProbability (MaxAmbigFraction));


  const Dna* DnaRev = MakeComplementary ();
  ASSERT (GoodObject (DnaRev));
  
  const size_t End = strlen (DnaRev->seq) - DnaRev->TrimBadStart (WindowLen, MaxAmbigFraction);
  
  delete DnaRev;
  

  return End;
}



size_t Dna::TrimLowComplexityStart () const
{        
  Float ThisMinComplexity = inf; 
  size_t Start = 0;
  ForString (i, seq)
  {
    if (i < 6)  // ??
      continue;

  //printf ("%3u: %c %0.1f\n", i + 1, C, Dna1. getComplexity ());
    if (   minimize (ThisMinComplexity, getComplexityInt (0, i))
        && LessEqualFloat (ThisMinComplexity, MinComplexity))
      Start = i;
  }
    

  return Start;
}



size_t Dna::TrimLowComplexityEnd () const
{        
  Float ThisMinComplexity = inf; 
  int SeqEnd = strlen (seq);    
  size_t End = SeqEnd;
  ForDown (i, SeqEnd - 6)  // ??
  {
  //printf ("%3u: %c %0.1f\n", i + 1, C, Dna1. getComplexity ());
    if (   minimize (ThisMinComplexity, getComplexityInt (i, SeqEnd))
        && LessEqualFloat (ThisMinComplexity, MinComplexity))
      End = i;
  }
    

  return End;
}
#endif



// Poly-nucleotide segment

bool Dna::polyNucWindow (size_t start,
           	             size_t windowLen,
           		           char nucleotide) const
{
  ASSERT (start <= seq. size ());
  ASSERT (charInSet (nucleotide, dnaAlphabet));

  const size_t end = min (start + windowLen, seq. size ());
  ASSERT (start <= end);
  if (start == end)
    return false;

  size_t mismatch = 0;
  FOR_START (size_t, i, start, end)
    if (! nucleotideMatch (seq [i], nucleotide))
      mismatch++;

  return (double) mismatch / (double) (end - start) <= 0.2;  // PAR
}



size_t Dna::findPolyNucEnd (size_t windowLen,
                            char nucleotide) const
{
  ASSERT (windowLen);
  ASSERT (charInSet (nucleotide, dnaAlphabet));


  size_t lastInitStart = seq. size () > windowLen ? seq. size () - windowLen : 0;
  size_t notPolyNuc = 0;
  FOR_REV (size_t, initStart, lastInitStart)
    if (polyNucWindow (initStart, windowLen, nucleotide))
    {
    	 lastInitStart = initStart;
    	 notPolyNuc = 0;
    }
    else
    {
      notPolyNuc++;
      if ((double) notPolyNuc > (double) windowLen * 1.3)  // PAR
        break;
    }


  FFOR_START (size_t, start, lastInitStart, seq. size ())
    if (polyNucWindow (start, min (windowLen, seq. size () - start), nucleotide))
    { 
      size_t nucleotide_num = 0;
      while (start < seq. size ())
      {
        if (nucleotideMatch (seq [start], nucleotide))
        {
          nucleotide_num++;
          if (nucleotide_num >= 2)  // PAR
            return start - nucleotide_num + 1;
        }
        else
          nucleotide_num = 0;
        start++;   
      }
      return start - nucleotide_num;
    }


  return seq. size ();
}



#if 0
void Dna::FindPolyAWindow (int   &BestStart,
                           int   &BestEnd,
                           Float &MaxWeight) const
{
  BestStart = -1;
  BestEnd   = -1;
  MaxWeight = 0;
 
 
  // ??  
  static const PROBABILITY BackgroundProb = 0.25;
  static const PROBABILITY SignalProb     = 0.90;
  ASSERT (SignalProb > BackgroundProb);

  static const Float SignalCoeff     = log (     SignalProb  /      BackgroundProb);
  static const Float BackgroundCoeff = log ((1 - SignalProb) / (1 - BackgroundProb));
  ASSERT (SignalCoeff     > 0);
  ASSERT (BackgroundCoeff < 0);
    
  
  const int Len = strlen (seq);
  if (Len < 1)
    return;
 
  
  const char Nuc = 'a';
 

  // Time: O(Len^2) ??
  For (Start, Len)
    if (seq [Start] == Nuc && 
        (Start == 0 || seq [Start - 1] != Nuc)  // Local maximum
       )
    {
      int M = 1;
      ForStart (End, Start + 1, Len)
        if (seq [End] == Nuc)
        {
          M++;
          const Float V = M * SignalCoeff + (End + 1 - Start - M) * BackgroundCoeff;
          if (Negative (V))
            break;
          if ((End == Len - 1 || seq [End + 1] != Nuc) &&  // Local maximum
              maximize (MaxWeight, V)
             )
          {
            BestStart = Start;
            BestEnd   = End + 1;
          }
        }
    }


  if (BestStart < 0)
    return;
    
    
  ForDown (i, BestStart - 1)
    if (NucleotideMatch (seq [i], Nuc))
      BestStart--;
    else
      break;
      
  ForStart (i, BestEnd, Len)
    if (NucleotideMatch (seq [i], Nuc))
      BestEnd++;
    else
      break;


  ASSERT (BestStart >= 0);
  ASSERT (BestStart < BestEnd);
  ASSERT (BestEnd <= Len);
}



bool Dna::DiNucLowComplexity (int   &BestStart,
                              int   &BestEnd,
                              Float &MaxWeight,
                              int   &MaxMonoNucLen) const
{
  BestStart     = -1;
  BestEnd       = -1;
  MaxWeight     = 0;
  MaxMonoNucLen = 0;
//TailWeight    = 0;
 
 
  // ??  
  static const PROBABILITY SameDiGoodProb = 0.30;
  static const PROBABILITY SameDiBadProb  = 0.80;
  ASSERT (SameDiBadProb > SameDiGoodProb);

  static const Float SameDiCoeff = log (     SameDiBadProb  /      SameDiGoodProb);
  static const Float DiffDiCoeff = log ((1 - SameDiBadProb) / (1 - SameDiGoodProb));
  ASSERT (SameDiCoeff > 0);
  ASSERT (DiffDiCoeff < 0);
    
  
  const int Len = strlen (seq);
  if (Len < 2)
    return false;


  const int DiLen = Len - 1;
  bool* SameDi = new bool [DiLen];
  For (i, DiLen)
    SameDi [i] = seq [i] == seq [i + 1]                     ||
                 CharInSet (seq [i],     DnaWildcards)      ||
                 (Qual != nullptr && Qual [i] <= Qual2N_Level) ||
                 CharInSet (seq [i + 1], DnaWildcards)      ||
                 (Qual != nullptr && Qual [i + 1] <= Qual2N_Level);

  // Time: O(DiLen^2) ??
  For (Start, DiLen)
    if (SameDi [Start] && 
        (Start == 0 || ! SameDi [Start - 1])  // Local maximum
       )
    {
      int M = 1;
      ForStart (End, Start + 1, DiLen)
        if (SameDi [End])
        {
          M++;
          const Float V = M * SameDiCoeff + (End + 1 - Start - M) * DiffDiCoeff;
          if (Negative (V))
            break;
          if ((End == DiLen - 1 || ! SameDi [End + 1]) &&  // Local maximum
              maximize (MaxWeight, V)
             )
          {
            BestStart = Start;
            BestEnd   = End;
          }
        }
    }


  bool Res = false;
  if (BestStart != -1 &&
      BestEnd - BestStart >= 18)  // Longest mononucleotide repeat in GenBank Reference genomic sequences ??
  {     
    int Mono = 0;
    ForStart (i, BestStart, BestEnd + 1)
      if (SameDi [i])
        Mono ++;
      else
      {
        maximize (MaxMonoNucLen, Mono);
        Mono = 0;
      }
    maximize (MaxMonoNucLen, Mono);

  #if 0    
    int M = 0;
    ForStart (i, BestEnd + 1, DiLen)
      if (SameDi [i])
        M++;
    TailWeight = M * SameDiCoeff + (DiLen + 1 - (BestEnd + 1) - M) * DiffDiCoeff;
  #endif
    
    BestEnd += 2;
    Res = true;
  }


  delete [] SameDi;  
  
  
  return Res;
}



size_t Dna::NucleotideFreq (NT_PROBABILITY acgtb,
                          size_t           Start,
                          size_t           End,
                          uint           Step) const
{
  ASSERT (Start <= End);
  ASSERT (Start <= strlen (seq));
  ASSERT (Step > 0);


  For (j, 5)
    acgtb [j] = 0;


  NT_PROBABILITY ACGTB1;
  size_t N = 0;
  for (size_t i = Start; i < End; i += Step)
    {
      ASSERT (seq [i] != '\0');
      Wild2NucleotideFreq (seq [i], ACGTB1);
      ASSERT (! EqualFloat (ACGTB1 [4], 1));
      For (j, 5)
        acgtb [j] += ACGTB1 [j];

      N++;
    }


  if (N > 0)
    For (j, 5)
      acgtb [j] /= N;


  return N;
}



size_t Dna::NucleotidePairFreq (PROBABILITY acgtb [5] [5]) const
{
  {
    For (j, 5)
    For (k, 5)
      acgtb [j] [k] = 0;
  }


  NT_PROBABILITY ACGTB1;
  NT_PROBABILITY ACGTB2;
  ForString (i, seq)
    {
      Wild2NucleotideFreq (seq [i], ACGTB2);
      ASSERT (! EqualFloat (ACGTB2 [4], 1));

      if (i > 0)
        {
          For (j, 5)
          For (k, 5)
            acgtb [j] [k] += ACGTB1 [j] * ACGTB2 [k];
        }

      {
        For (j, 5)
          ACGTB1 [j] = ACGTB2 [j];
      }
    }


  const size_t Len = strlen (seq);
  if (Len >= 1)
    {
      For (j, 5)
      For (k, 5)
        acgtb [j] [k] /= (Len - 1);
    }


  return Len;
}



// MeltTemp

// entropy_value [i] [j] = entropy_value [3-j] [3-i]  
static const int entropy_value [4] [4] /* acgt */ =
  {// a    c    g    t
    {240, 173, 208, 239},  // a
    {129, 266, 278, 208},  // c
    {135, 267, 266, 173},  // g
    {169, 135, 129, 240}   // t
  };



static Float MeltTempEntropy (const char* DnaSeq)
{
  Float S = 108.0;  // For non-self-complementary sequences; 124 for self-complementary sequences
  const int Len = strlen (DnaSeq);
  for (int i = 0; i < Len - 1; i++)
    {
      NT_PROBABILITY ACGTB1, ACGTB2;
      Wild2NucleotideFreq (DnaSeq [i],     ACGTB1);
      Wild2NucleotideFreq (DnaSeq [i + 1], ACGTB2);
      For (j, 4)
        For (k, 4)
          S += entropy_value [j] [k] * ACGTB1 [j] * ACGTB2 [k];
    }

  return - S * 0.1;
}



// enhalpy_value [i] [j] = enthalpy_value [3-j] [3-i]  
static const int enthalpy_value [4] [4] /* acgt */ =
  {
    {91,  65,  78, 86},
    {58, 110, 119, 78},
    {56, 111, 110, 65},
    {60,  56,  58, 91}
  };



static Float MeltTempEnthalpy (const char* DnaSeq)
{
  Float H = 0.0;
  const int Len = strlen (DnaSeq);
  for (int i = 0; i < Len - 1; i++)
    {
      NT_PROBABILITY ACGTB1, ACGTB2;
      Wild2NucleotideFreq (DnaSeq [i],     ACGTB1);
      Wild2NucleotideFreq (DnaSeq [i + 1], ACGTB2);
      For (j, 4)
        For (k, 4)
          H += enthalpy_value [j] [k] * ACGTB1 [j] * ACGTB2 [k];
    }


  return - H * 0.1;
}



Float Dna::GetMeltTemp (Float DnaConc,
                        Float SaltConc,
                        Float MgConc) const
// T_m (reverse-complement of seq) = T_m (seq)
{
  const Float T_0       = 273.15;
  const Float R         = 1.987;  // Universal gas constant
  const Float logdna    = R * log (DnaConc / (4.0 * 1.0e9));
  // Shu-ichi Nakano et al., Nucleic Acids Research, 1999, vol. 27, No. 14 2957-2965: Data
  const Float SaltEquiv = 30.8 * MgConc;  
  const Float logsalt   = 16.6 * log10 ((SaltConc + SaltEquiv) / 1.0e3);

  return (MeltTempEnthalpy (seq) * 1000.0) / (MeltTempEntropy (seq) + logdna) - T_0 + logsalt;
}



PRIMER* Dna::MakePrimer (size_t        PrimerLen,
                         size_t        SearchStart,
                         size_t        SearchEnd,
                         PHRED_SCORE MinQuality,
                         size_t        MaxMaxRepeat,
                         PROBABILITY MaxGCShare,
                         Float       DnaConcentration,
                         Float       SaltConcentration,
                         Float       MgConcentration,
                         Float       MinMeltTemp,
                         Float       MaxMeltTemp,
                         Float       TargetMeltTemp) const
{
  ASSERT (PrimerLen > 0);
  ASSERT (SearchStart <= SearchEnd);
  ASSERT (BetweenFloat (TargetMeltTemp, MinMeltTemp, MaxMeltTemp));

  
  minimize (SearchEnd, (size_t) strlen (seq));

  if ((int) PrimerLen > (int) SearchEnd - (int) SearchStart)
    return nullptr;


  PRIMER* BestPrimer = nullptr;
  Float BestMeltTempDiff = inf;
  ForStart (i, SearchStart, SearchEnd - PrimerLen)
    {
      PRIMER Primer (this, i, PrimerLen, 
                     MinQuality, DnaConcentration, SaltConcentration, MgConcentration);
      
      if (! Primer. Ambiguity &&
          Primer. MaxRepeat <= MaxMaxRepeat &&
          LessEqualFloat (Primer. GCShare, MaxGCShare) &&
          BetweenFloat (Primer. MeltTemp, MinMeltTemp, MaxMeltTemp) &&
          minimize (BestMeltTempDiff, AbsFloat (Primer. MeltTemp - TargetMeltTemp)))
        BestPrimer = (PRIMER*) Primer. Copy ();
    }


  return BestPrimer;
}



PROBABILITY Dna::NucleotideFreqFrame (size_t        Start,
                     			              Frame       frame,
                     			              size_t        Count,
                     			              const char* Alphabet) const
{
  const size_t Len = strlen (seq);
  ASSERT (Start < Len);
  ASSERT (isFrame (frame));
  ASSERT (Alphabet != nullptr);


  int _Frame = AbsInt (frame) - 1;
  ASSERT (_Frame >= 0);
  if (frame < 0)
    _Frame = - _Frame;
  const int Step = frame > 0 ? 3 : -3;
  size_t N = 0;
  size_t C = 0;
  for (int i = Start + _Frame;
       C < Count && Between (i, 0, Len);
       i += Step, C++)
    if (CharInSet (seq [i], Alphabet))
      N++;


  ASSERT (N <= C);
  return (Float (N) / C);
}



static void Frame2Prob (Frame       frame,
                     	  PROBABILITY Prob [3] [4])
// Output: Prob [] []
{
  ASSERT (isFrame (frame));
  ASSERT (frame >= 1);


  // For all dicots -??
  const PROBABILITY initProb [3] [4] =   // frame, "acgt"
    {
      {0.29, 0.19, 0.32, 0.20},
      {0.31, 0.23, 0.19, 0.27},
      {0.24, 0.22, 0.22, 0.32}
    };

  For (i, 3)
    For (j, 4)
      Prob [i] [j] = initProb [(3 - (frame - 1) + i) % 3] [j] * 4 /*to avoid machine 0*/;
}



static Float _GetFrameProb (const uint        Y    [3] [4],
                    			  const PROBABILITY Prob [3] [4])
// Return: probability * Const
{
  PROBABILITY T = 1;
  For (i, 3)
    T *= Multinomial (4, Prob [i], Y [i], true);


  return T;
}



static void GetFrameProb (const uint  Y [3] [4],
                     			PROBABILITY FrameProb [3])
{
  PROBABILITY Prob [3] [4];

  Float S = 0.0;  // ??
  for (Frame frame = 1; frame <= 3; frame++)
    {
      Frame2Prob (frame, Prob);
      FrameProb [frame - 1] = _GetFrameProb (Y, Prob);
      S += FrameProb [frame - 1];
    }

/* ??
  For (i, 3)
    FrameProb [i] /= S;
*/
}



void Dna::GetOpenFrameProb (size_t        Start,
                     			  size_t        Count,
                     			  PROBABILITY FrameProb [3]) const
{
  // Y
  uint Y [3] [4];  // frame, "acgt"
  for (Frame frame = 1; frame <= 3; frame++)
    {
    //printf ("%d: ", frame);  
      ForString (i, DnaAlphabet)  // IUPAC -??
        {
          char S [3] = "  ";
          S [0] =          DnaAlphabet [i];
          S [1] = ToUpper (DnaAlphabet [i]);
          const PROBABILITY Prob = NucleotideFreqFrame (Start, frame, Count, S);
          Y [frame - 1] [i] = Round (Prob * Count);
        //printf (" %0.3lf (%u)", Prob, Y [frame - 1] [i]);  
        }
    //printf ("\n");  
    }


  GetFrameProb (Y, FrameProb);

  // ??
  // Assume: Start = 3 * Count < strlen (seq)
  For (i, 3)
    FrameProb [i] = log (FrameProb [i]) / (3 * Count) * 1000;
}



Float Dna::CompareOrfs (const CodStatistics* StatTable1,
                        const CodStatistics* StatTable2) const
{
  ASSERT (StatTable1 != nullptr);
  ASSERT (StatTable2 != nullptr);


  try
    {
      CrsDna dna (seq);
      dna. makeLower ();

    		CrsOrf Orf1, Orf2;
    		if (! dna. compareScores (*StatTable1, *StatTable2, & Orf1, & Orf2))
    		  return 0.0;

		    const Float Diff = Orf1. score () - Orf2. score ();
		    ASSERT (! isnan (Diff) && finite (Diff));
		    return Diff;
    }
  catch (...)
    ERROR;
}



Float Dna::Orf2Strand (const CodStatistics* StatTable) const
{
  ASSERT (StatTable != nullptr);


  try
    {
      CrsDna dna (seq);
      dna. makeLower ();

    		CrsOrf Orf, RevOrf;
    		if (! dna. compareStrandScores (*StatTable, & Orf, & RevOrf))
    		  return 0.0;

		    const Float Diff = Orf. score () - RevOrf. score ();
		    ASSERT (! isnan (Diff) && finite (Diff));
		    return Diff;
    }
  catch (...)
    ERROR;
}



bool Dna::CDS_Start_nr2Orf (int CDS_Start_nr,
                            int Neighborhood,
                            int &CDS_Start,
                            int &CDS_End) const
{
  CDS_Start = 0;
  CDS_End   = 0;


  const int DnaLen = strlen (seq);
  For (i, (2 * Neighborhood) / 3)  // aa
    for (int j = -1; j <= 1; j += 2)
      {
        const int Pos = CDS_Start_nr + (int) i * 3 * j;
        if (BetweenEqual (Pos, 0, DnaLen - 3) &&
            CodonMatch (& seq [Pos], "atg"))
          {
            CDS_Start = Pos;

            // CDS_End            
            CDS_End = CDS_Start + 3;
            while (CDS_End <= DnaLen &&
                   ! TerminatorCodon (& seq [CDS_End - 3]))
              CDS_End += 3;
            if (CDS_End > DnaLen)
              return false;
              
            return true;
          }
      }


  return false;
}

 

bool Dna::CDS_End_nr2Orf (int CDS_End_nr,
                          int Neighborhood,
                          int &CDS_Start,
                          int &CDS_End) const
{
  CDS_Start = 0;
  CDS_End   = 0;


  const int DnaLen = strlen (seq);
  For (i, (2 * Neighborhood) / 3)  // aa
    for (int j = -1; j <= 1; j += 2)
      {
        const int Pos = (CDS_End_nr - 3) + (int) i * 3 * j;
        if (BetweenEqual (Pos, 0, DnaLen - 3) &&
            MayBeTerminatorCodon (& seq [Pos]))
          {
            CDS_End = Pos + 3;

            // CDS_Start            
            CDS_Start = MAXINT;
            int k = Pos - 3;
            while (k >= 0 &&
                   ! TerminatorCodon (& seq [k]))
              {
                if (Codon2AA (& seq [k]) == 'M')
                  CDS_Start = k;
                k -= 3;
              }
            if (CDS_Start == MAXINT)
              return false;
              
            return true;
          }
      }


  return false;
}

 


////////////////////////////////////////////////////////////

Dna* ReadDna (const char* FName,
       	      bool        GapCoded)
{
  FILE_INPUT F (10 * 1024, FName);
  if (F. ErrorNum != 0)
    return nullptr;


  F. NextLine ();
  ASSERT (! F. Eof &&
	         StrIsLeft (F. Line, ">"));

  Dna* dna = new Dna (F, GapCoded);
  ASSERT (GoodObject (dna));


  return dna;
}



///////////////////////////// DNA_COLLECTION /////////////////////////////

DNA_COLLECTION::DNA_COLLECTION (const char* FName,
                                bool        GapCoded):
  inherited ()
{
  FILE_INPUT F (255, FName);
  ASSERT (F. ErrorNum == 0);


  F. NextLine ();
  while (! F. Eof &&
       	 StrIsLeft (F. Line, ">"))
    {
      Dna* dna = new Dna (F, GapCoded);
      ASSERT (GoodObject (dna));

      if (StrIsEmpty (dna->seq))
        delete dna;
      else
        Append (dna);
    }
}



void DNA_COLLECTION::SaveQual (const char* FName,
                               uint        DefaultQualScore) const
{
  ASSERT (! StrIsBlank (FName));


  FILE* F = fopen (FName, "w");
  ASSERT (F != nullptr);

  ForCollection (i, *this)
    GetDna (i) -> QualSaveFile (F, DefaultQualScore);

  fclose (F);
}
#endif




/////////////////////////////// SubstMat ///////////////////////////////

SubstMat::SubstMat (const string &fName)
{
  const string charsS (chars);
  
  for (const char c : charsS)
    ASSERT (charInSet (c, extTermPeptideAlphabet));
  
  FOR (size_t, i, sim_size)
    FOR (size_t, j, sim_size)
      sim [i] [j] = numeric_limits<AlignScore>::quiet_NaN ();
  
  LineInput f (fName);
  size_t charNum = 0;
  while (f. nextLine ())
  {
    const char c1 = f. line [0];
    if (c1 == '#')  // Comment
      continue;
    if (c1 == ' ')
    {
      QC_ASSERT (! charNum);
      replaceStr (f. line, " ", "");
      QC_ASSERT (f. line == charsS);
    }
    else 
    {
      QC_ASSERT (charNum < charsS. size ());
      QC_ASSERT (c1 == charsS [charNum]);
      findSplit (f. line);
      for (const char c : charsS)
      {
        ASSERT ((size_t) c < sim_size);
        trimLeading (f. line);
        const string numS (findSplit (f. line));
        QC_ASSERT (! numS. empty ());
        const AlignScore r = atof (numS. c_str ());
        QC_ASSERT (r == r);  // Not a NaN
        sim [(size_t) c1] [(size_t) c] = r;
      }
      charNum++;
    }
  }
  QC_ASSERT (charNum == charsS. size ());
}



void SubstMat::qc () const
{
  if (! qc_on)
    return;
    
  // sim[][]
  FOR (size_t, i, sim_size)
  {
    if (! goodChar (i))
      continue;
    QC_IMPLY (! isAmbigAa (char (i)), sim [i] [i] >= 0.0);
    FOR (size_t, j, sim_size)
    {
      if (! goodChar (j))
        continue;
      QC_ASSERT (fabs (sim [i] [j] - sim [j] [i]) < 1e-6);  // PAR
      if (sim [i] [i] < sim [i] [j])
        cout << "sim(" << char (i) << ',' << char (i) << ") < sim(" << char(i) << ',' << char(j) << ")" << endl;          
    }
  }
}  



void SubstMat::saveText (ostream& os) const 
{
  FOR (size_t, i, sim_size)
    if (goodChar (i))
      os << ' ' << char (i);
  os << endl;
  
  FOR (size_t, i, sim_size)
    if (goodChar (i))
    {
      os << char (i);
      FOR (size_t, j, sim_size)
        if (goodChar (j))
          os << ' ' << sim [i] [j];
      os << endl;
    }
}





///////////////////////////////// PeptideOrf //////////////////////////////////

PeptideOrf::PeptideOrf (size_t translationStart_arg,
              	        Strand strand_arg,
              			    const Peptide* peptide,
                        size_t start_arg)
: translationStart (translationStart_arg)
, strand (strand_arg)
, start (start_arg)
, startM (Peptide::isStartAa (peptide->seq [start]))
, stop (peptide->seq. size ())
, stopTerminator (false)
{
  ASSERT (isStrand (strand));
	ASSERT (peptide);
  ASSERT (! peptide->seq. empty ());
  ASSERT (! peptide->sparse);
  ASSERT (start < peptide->seq. size ());
	
  // stop, stopTerminator
  FOR_START (size_t, i, start, peptide->seq. size ())
    if (peptide->seq [i] == *terminator)
    {
     	stop = i;
     	stopTerminator = true;
     	break;
    }
  ASSERT (stop <= peptide->seq. size ());

  ASSERT (! empty ());
}



void PeptideOrf::qc () const
{
  if (! qc_on)
    return;

	if (empty ())
		return;

  QC_ASSERT (isStrand (strand));
  QC_IMPLY (! strand, ! translationStart);  	

  QC_ASSERT (start <= stop);
  QC_IMPLY (start == stop, ! startM && stopTerminator);
  QC_IMPLY (start, startM);
}



Peptide* PeptideOrf::toPeptide (const Peptide* peptide) const
{
	ASSERT (peptide);
  ASSERT (! peptide->sparse);
	ASSERT (! empty ());
  ASSERT (stop <= peptide->seq. size ());
  IMPLY (stop < peptide->seq. size (), stopTerminator && peptide->seq [stop] == *terminator);
  IMPLY (stop == peptide->seq. size (), ! stopTerminator);
	
  Peptide* p = new Peptide (peptide->name, peptide->seq, false);
  p->seq. erase (stop + stopTerminator);
  p->seq. erase (0, start);
  
  return p;
}




/////////////////////////////// Peptide //////////////////////////////////

size_t aa2num (char wildAminoacid)
{
  const char c = toUpper (wildAminoacid);

  size_t i = 0;
  for (const char* a = peptideAlphabet; *a; a++)
    if (c == *a)
      return i;
    else
    	i++;
  if (wildAminoacid == *terminator)
    return 20;
  if (isAmbigAa (c))
    return 21;

  NEVER_CALL;
}



bool moreGeneralAminoacid (char wildAminoacid1,
                           char wildAminoacid2)
{
  if (wildAminoacid1 == wildAminoacid2)
    return true;


  switch (wildAminoacid1)
  {
    case 'x': return true;
    case 'b': return charInSet (wildAminoacid2, "dn");
    case 'z': return charInSet (wildAminoacid2, "eq");
    case 'j': return charInSet (wildAminoacid2, "il");
    case 'u': return wildAminoacid2 == '*'; 
    case 'o': return wildAminoacid2 == '*'; 
    case '*': return charInSet (wildAminoacid2, "uo");
  }

  return false;
}




// Peptide

Peptide::Peptide (Multifasta &fasta,
                  size_t reserveLen,
                  bool sparse_arg)
: Seq (fasta. in, reserveLen, sparse_arg, true) 
{ 
	ASSERT (fasta. aa);
	fasta. prog (getId ());
}



void Peptide::qc () const
{
  if (! qc_on)
    return;
  Seq::qc ();

	QC_IMPLY (hasInsideStop (), pseudo);
}



bool Peptide::isDescriptionPartial () const
{
	string desc (getDescription (false));
	strLower (desc);
	// PAR
	return    contains (desc, "fragmented")  // "fragment" may mean a complex unit
	       || contains (desc, "fragmentary")
	       || contains (desc, "partial")
	       || contains (desc, "truncat");
}



Vector<PeptideOrf> Peptide::getPeptideOrfs (size_t translationStart,
              	                            Strand strand,
                                            bool includeInitial,
                                            bool longestOnly,
                                            size_t len_min) const
{
  Vector<PeptideOrf> orfs;

  PeptideOrf prev;
  FFOR (size_t, i, seq. size ())
    if (isStartAa (seq [i]))
    {
      const PeptideOrf orf (translationStart, strand, this, i);
      ASSERT (! orf. empty ());
      orf. qc ();
      if (! prev. empty ())
      {
        ASSERT (prev. start < orf. start);
        if (longestOnly && prev. stop == orf. stop)
          continue;
      }
      if (orf. size () >= len_min)
   	    orfs << orf;
   	  prev = orf;
   	}

  if (includeInitial && (orfs. empty () || orfs [0]. start > 0))
  {
    const PeptideOrf orf (translationStart, strand, this, 0);
    ASSERT (! orf. empty ());
    orf. qc ();
    if (orf. size () >= len_min)
  	  orfs << orf;
  }

  return orfs;
}



double Peptide::getComplexityInt (size_t start,
                                  size_t end) const
{
  if (seq. empty ())
    return 0;

  const size_t maxAANum = 21;

  // diFreq[][], n
  size_t diFreq [maxAANum] [maxAANum];
  FOR (size_t, i, maxAANum)
    FOR (size_t, j, maxAANum)
      diFreq [i] [j] = 0;
  size_t n = 0;
  minimize (end, seq. size ());
  FOR_START (size_t, i, start + 1, end)
    {
      const size_t a = aa2num (seq [i - 1]);
      const size_t b = aa2num (seq [i]);
      if (   a < maxAANum 
          && b < maxAANum
         )
      {
        diFreq [a] [b] ++;
        n++;
      }
    }

  // Entropy
  double s = 0;
  FOR (size_t, i, maxAANum)
    FOR (size_t, j, maxAANum)
      if (diFreq [i] [j])
      {
       	const double p = (double) diFreq [i] [j] / (double) n;
       	s -= p * log (p);
      }
  ASSERT (s >= 0);

  return s;
}



double Peptide::getSelfSimilarity (const SubstMat &mat,
                                   size_t start,
                                   size_t stop) const
{
  ASSERT (start <= stop);
  ASSERT (stop <= seq. size ());
  
  if (! start && ! stop)
    stop = seq. size ();
  
  double s = 0.0;
  FOR_START (size_t, i, start, stop)
  {
    const char c = seq [i];
    if (c != '-')
      s += mat. sim [(size_t) c] [(size_t) c];
  }

  return s;
}



double Peptide::getSimilarity (const Peptide &other,
                               const SubstMat &mat,
                               double gapOpenCost,
                               double gapCost) const
{
  ASSERT (gapOpenCost >= 0.0);
  ASSERT (gapCost >= 0.0);
  ASSERT (sparse);
  ASSERT (other. sparse);
  ASSERT (seq. size () == other. seq. size ());
  
  double s = 0.0;
  char c1_prev = '\0';
  char c2_prev = '\0';
  FFOR (size_t, i, seq. size ())
  {
    const char c1 =        seq [i];
    const char c2 = other. seq [i];
    if (   c1 == '-' 
        && c2 == '-'
       )
      ;
    else if (   c1 == '-' 
             || c2 == '-'
            )
    {
      s -= gapCost;
      if (   (c1 == '-' && c1_prev != '-')
          || (c2 == '-' && c2_prev != '-')
         )
        s -= gapOpenCost;
    }
    else
      s += mat. sim [(size_t) c1] [(size_t) c2];
    c1_prev = c1;
    c2_prev = c2;
  }

  return s;
}



size_t Peptide::ambig2X ()
{
  size_t n = 0;
  for (char& c : seq)
    if (charInSet (c, "BZJUO"))
    {
      c = 'X';
      n++;
    }
  
  return n;
}



void Peptide::toGBMR4 ()
{
  for (char& c : seq)
    switch (c)
    {
      case 'A':
      case 'D':
      case 'K':
      case 'E':
      case 'R':
      case 'N':
      case 'T':
      case 'S':
      case 'Q': c = 'A'; break;
      case 'Y':
      case 'F':
      case 'L':
      case 'I':
      case 'V':
      case 'M':
      case 'C':
      case 'W':
      case 'H': c = 'Y'; break;
      case 'G':
      case 'P':
      case 'X':
      case '*': break;
      default: ERROR_MSG ("Unknown aa " + string (1, c));
    }
}



#if 0
size_t Peptide::GetLeftMPos (size_t Start) const
{
  ASSERT (Start < strlen (seq));
  ForDown (i, Start)
    if (seq [i] == 'M')
      return i;

  return no_index;
}



size_t Peptide::GetClosestMPos (size_t Start) const
{
  size_t ClosestMPos = no_index;
  size_t Distance = MAXUINT;
  ForString (i, seq)
    if (seq [i] == 'M' &&
        minimize (Distance, (size_t) AbsInt ((int) Start - (int) i)))
      ClosestMPos = i;

  return ClosestMPos;
}



Peptide* ReadPeptide (const char* FName)
{
  FILE_INPUT F (255, FName);
  if (F. ErrorNum != 0)
    return nullptr;


  F. NextLine ();
  ASSERT (! F. Eof &&
	 StrIsLeft (F. Line, ">"));

  Peptide* peptide = new Peptide (F);
  ASSERT (GoodObject (peptide));


  return peptide;
}




// PEPTIDE_COLLECTION

PEPTIDE_COLLECTION::PEPTIDE_COLLECTION (const char* FName):
  inherited ()
{
  FILE_INPUT F (255, FName);
  ASSERT (F. ErrorNum == 0);


  F. NextLine ();
  while (! F. Eof &&
       	 StrIsLeft (F. Line, ">"))
    {
      Peptide* peptide = new Peptide (F);
      ASSERT (GoodObject (peptide));

      Append (peptide);
    }
}




///////////////////////////// SEQ_CLUSTER_FILE ////////////////////////////////

// SEQ_CLUSTER_ITEM

SEQ_CLUSTER_ITEM::SEQ_CLUSTER_ITEM (const char* initName,
                        		    		    int         initLeftOffset,
                        		    		    size_t        initLen,
                        		    		    size_t        initMatchNum):
  inherited  (initName),

  LeftOffset (initLeftOffset),
  Len        (initLen),
  MatchNum   (initMatchNum)
{
  ASSERT (! StrIsEmpty (Name));
  ASSERT (Len > 0);
  ASSERT ((int) MatchNum >= 0);
}



void SEQ_CLUSTER_ITEM::Print () const
{
  printf ("%s %d %u %u\n", Name, LeftOffset, Len, MatchNum);
}




// SEQ_CLUSTER_FILE

SEQ_CLUSTER_FILE::SEQ_CLUSTER_FILE (const char* FName):
  inherited  (),

  F          (512, FName),
  ClusterNum (0),
  ItemNum    (0)
{
  ASSERT (F. ErrorNum == 0);
}



void SEQ_CLUSTER_FILE::FirstCluster ()
{
  while (F. NextLine () &&
  	      StrIsLeft (F. Line, "#"));

  ASSERT (! F. Eof);
  ASSERT (F. BlankLine ());
}



bool SEQ_CLUSTER_FILE::NextCluster ()
{
  if (! F. NextLine ())
    {
      fprintf (stderr, "\n");
      return false;
    }


  uint N;
  if (sscanf (F. Line, "Cluster %u", & N) != 1)
    ERROR;

  ClusterNum++;
//ASSERT (N == ClusterNum);
  fprintf (stderr, "\r# Clusters = %u", ClusterNum);
  ItemNum = 0;


  return true;
}



SEQ_CLUSTER_ITEM* SEQ_CLUSTER_FILE::NextItem ()
{
  ASSERT (F. NextLine ());


  if (StrIsLeft (F. Line, "total = "))
    {
      ASSERT (ItemNum >= 1);

      uint M;
      if (sscanf (F. Line, "total = %u", & M) != 1)
       	ERROR;
      ASSERT (M == ItemNum);

      return nullptr;
    }


  ItemNum++;

  char Name [255];
  int  LeftOffset;
  size_t Len, MatchNum;
  if (sscanf (F. Line, "%s %d %u %u",
      	       Name, & LeftOffset, & Len, & MatchNum) != 4)
    ERROR;

  SEQ_CLUSTER_ITEM* Item = new SEQ_CLUSTER_ITEM (Name, LeftOffset, Len, MatchNum);
  ASSERT (GoodObject (Item));


  return Item;
}
#endif




// Cds

void Cds::qc () const
{
  if (! qc_on)
    return;
  QC_ASSERT (! refProt. empty ());
	QC_ASSERT (start != stop);
	QC_ASSERT (divisible ((uint) size (), 3)); 
	QC_ASSERT (sizeEffective () > 0);
	QC_ASSERT (positives > 0.5);  // PAR
	QC_ASSERT (positives <= 1);
}



size_t Cds::getOverlap_ (const Cds &other) const
{ 
  ASSERT (left () <= other. left ());
	ASSERT (other. left () <= right ());
	if (other. right () <= right ())
		return other. size ();
	return right () - other. left ();
}



bool Cds::coexists (const Cds &next) const
{ 
  ASSERT (left () <= next. left ());
  if (! strand () && next. strand ())
    return right () + promoter_min <= next. left ();
  const size_t overlap_min = min ( size_min - 3/*PAR*/
				    	                   , (size_t) ((double) min (size (), next. size ()) * 0.2/*PAR*/)
                                 );
  return getOverlap (next) <= overlap_min;
}



int Cds::getLenIncrease (const Cds &prev) const  
{ 
  ASSERT (prev. left () <= left ());
	return sizeEffective () /*- (int) getOverlap (prev) */; 
}



bool Cds::operator< (const Cds &other) const
{
  LESS_PART (*this, other, left ());
  LESS_PART (other, *this, right ());
  LESS_PART (other, *this, positives);
  return refProt < refProt;
}




// DnaAnnot

const Cds* DnaAnnot::run ()
{
  sort (cdss);
  
  
  // Dynamic programming
  Vector<const Cds*> nexts;
  FFOR (size_t, i, cdss. size ())
  {
  	// cdss[j] are optimized for all j < i
  	
  	Cds& cds = cdss [i];
  	IMPLY (i, cds. left () >= cdss [i - 1]. left ());

  	// Optimization of cdss[i]
    // cds.{sumLen,bestPrev}    	
  	ASSERT (! cds. bestPrev);
  	cds. sumLen = cds. sizeEffective ();
  	for (const Cds* prev : cds. prevs)
  	{
  	  const long sumLen_new = prev->sumLen  + cds. getLenIncrease (*prev);
  	  if (   cds. sumLen < sumLen_new
  	      || (   cds. sumLen == sumLen_new   // tie resolution
  	        #if 1
  	          && prev < cds. bestPrev
  	        #else
  	          && cds. bestPrev 
  	          && prev->refProt < cds. bestPrev->refProt
  	        #endif
  	         )
  	     )  
  	  {
  	    cds. sumLen = sumLen_new;
  	  	cds. bestPrev = prev;
  	  }
  	}
  	IMPLY (! cds. prevs. empty (), cds. bestPrev && cds. sumLen > 0);
  	

    // Cds::prevs
  	nexts. clear ();  // Cds::coexists(cds)
  	FOR_START (size_t, j, i + 1, cdss. size ())
  	{
  		Cds& next = cdss [j];
  		if (! cds. coexists (next))
  			continue;
  		
  		bool tooFar = false;
  	  const int nextLenIncrease = next. getLenIncrease (cds);
  	  for (const Cds* next_old : nexts)
  	    if (   next_old->coexists (next)
  	    	  && next. getLenIncrease (*next_old) + next_old->getLenIncrease (cds) > nextLenIncrease
  	    	 )
  	    {
  	    	tooFar = true;
  	    	break;
  	    }
  	  if (tooFar)
  	  	break;
  	  	
  	  nexts << & next;
  		next. prevs << & cds;
  	}
  }


  const Cds* cds_last = nullptr;
  long sumLen_max = 0;
  for (const Cds& cds : cdss)
    if (maximize (sumLen_max, cds. sumLen))
    	cds_last = & cds;


  return cds_last;
}




// Align

namespace 
{
  
struct Cell
{
  enum Dir {dirLeft, dirUp, dirDiag};
    // Direction to the previous Cell
  typedef  array <AlignScore, 3/*Dir*/>  Score;
  typedef  array <Dir,    3/*Dir*/>  Dir_best;
    // Direction from the previous Cell to the previous previous Cell
  Score scores;
  Dir_best dirs;

  
  Cell () 
    { constexpr AlignScore score_inf = numeric_limits<AlignScore>::infinity ();
      scores [dirLeft] = - score_inf;
      scores [dirUp]   = - score_inf;
      scores [dirDiag] = - score_inf;
    }
};
  
}



Align::Align (const Peptide &pep1,
	            const Peptide &pep2,
	            const SubstMat &substMat,
	            AlignScore gapOpenCost,
	            AlignScore gapCost,
	            size_t semiglobalMatchLen_min)
{
  ASSERT (! pep1. sparse);
  ASSERT (! pep2. sparse);
  ASSERT (gapOpenCost >= 0.0);
  ASSERT (gapCost >= 0.0);
  
  
#if 0
  // PAR
	int gap_open   = 1;  
	int gap_extent = 1;
	if (blosum62)
	{
	  gap_open   = -11;  
		gap_extent =  -2; 
	}
	else
	{
	  gap_open   = -8;  
		gap_extent = -2;  
	}
  ASSERT (gap_open <= 0);
  ASSERT (gap_extent < 0);

	CNWAligner al (pep1. seq, pep2. seq, blosum62 ? & NCBISM_Blosum62 : & NCBISM_Pam30);
  al. SetWg (gap_open);
  al. SetWs (gap_extent);
	al. SetEndSpaceFree (semiglobal, semiglobal, semiglobal, semiglobal);
	score = al. Run ();
  tr = al. GetTranscriptString ();  
  
	finish (pep1, pep2, semiglobal, match_len_min);
#endif


  const string seq1 ("-" + pep1. seq);
  const string seq2 ("-" + pep2. seq);


  Vector<Vector<Cell>> cells (seq1. size ());
  for (Vector<Cell>& rowVec : cells)
    rowVec = Vector<Cell> (seq2. size ());
  FFOR (size_t, row_, seq1. size () + seq2. size () - 1)
  {
    size_t row = row_;
    size_t col = 0;
    if (row_ >= seq1. size ())
    {
      row = seq1. size () - 1;
      col = row_ - row;
    }
    
    cells [0] [0]. scores [Cell::dirDiag] = 0.0;
    while (col < seq2. size ())
    {
      Cell::Score&    scores_cur  = cells [row] [col]. scores; 
      Cell::Dir_best& dirs_cur    = cells [row] [col]. dirs; 
      if (col)
      {
        AlignScore&    score_cur = scores_cur [Cell::dirLeft];
        Cell::Dir& dir_cur   = dirs_cur   [Cell::dirLeft];
        const Cell::Score& prev = cells [row] [col - 1]. scores;
        if (semiglobalMatchLen_min && (! row || row == seq1. size () - 1))
        {
          score_cur = prev [Cell::dirLeft];
          dir_cur = Cell::dirLeft;
          if (maximize (score_cur, prev [Cell::dirUp]))
            dir_cur = Cell::dirUp;
          if (maximize (score_cur, prev [Cell::dirDiag]))
            dir_cur = Cell::dirDiag;
        }
        else
        {
          score_cur = prev [Cell::dirUp];
          dir_cur = Cell::dirUp;
          if (maximize (score_cur, prev [Cell::dirDiag]))
            dir_cur = Cell::dirDiag;
          score_cur -= gapOpenCost;
          if (maximize (score_cur, prev [Cell::dirLeft]))
            dir_cur = Cell::dirLeft;
          score_cur -= gapCost;
        //score_cur= - gapCost + max (prev [Cell::dirLeft], max (prev [Cell::dirUp], prev [Cell::dirDiag]) - gapOpenCost);
        }
      }
      if (row)
      {
        AlignScore&    score_cur = scores_cur [Cell::dirUp];
        Cell::Dir& dir_cur   = dirs_cur   [Cell::dirUp];
        const Cell::Score& prev = cells [row - 1] [col]. scores;
        if (semiglobalMatchLen_min && (! col || col == seq2. size () - 1))
        {
          score_cur = prev [Cell::dirLeft];
          dir_cur = Cell::dirLeft;
          if (maximize (score_cur, prev [Cell::dirUp]))
            dir_cur = Cell::dirUp;
          if (maximize (score_cur, prev [Cell::dirDiag]))
            dir_cur = Cell::dirDiag;
        }
        else
        {
          score_cur = prev [Cell::dirLeft];
          dir_cur = Cell::dirLeft;
          if (maximize (score_cur, prev [Cell::dirDiag]))
            dir_cur = Cell::dirDiag;
          score_cur -= gapOpenCost;
          if (maximize (score_cur, prev [Cell::dirUp]))
            dir_cur = Cell::dirUp;
          score_cur -= gapCost;        
        //score_cur = - gapCost + max (prev [Cell::dirUp], max (prev [Cell::dirLeft], prev [Cell::dirDiag]) - gapOpenCost);
        }
      }
      if (row && col)
      {
        AlignScore&    score_cur = scores_cur [Cell::dirDiag];
        Cell::Dir& dir_cur   = dirs_cur   [Cell::dirDiag];
        const Cell::Score& prev = cells [row - 1] [col - 1]. scores;
        score_cur = prev [Cell::dirLeft];
        dir_cur = Cell::dirLeft;
        if (maximize (score_cur, prev [Cell::dirUp]))
          dir_cur = Cell::dirUp;
        if (maximize (score_cur, prev [Cell::dirDiag]))
          dir_cur = Cell::dirDiag;
        score_cur += substMat. sim [size_t (seq1 [row])] [size_t (seq2 [col])];
      #if 0
        score_cur =   substMat. sim [size_t (seq1 [row])] [size_t (seq2 [col])] 
                    + max (prev [Cell::dirUp], max (prev [Cell::dirLeft], prev [Cell::dirDiag]));
      #endif
      }
      
      if (! row)
        break;
      row--;
      col++;
    }
  }
  
  
  transformations. reserve (pep1. seq. size () + pep2. seq. size ());
  {
    size_t row = seq1. size () - 1;
    size_t col = seq2. size () - 1;
    Cell::Dir dir = Cell::dirUp;
    {
      const Cell::Score& scores_last = cells [row] [col]. scores;
      score = scores_last [Cell::dirUp];
      if (maximize (score, scores_last [Cell::dirLeft]))
        dir = Cell::dirLeft;
      if (maximize (score, scores_last [Cell::dirDiag]))
        dir = Cell::dirDiag;
    //score = max (score_last [Cell::dirUp], max (score_last [Cell::dirLeft], score_last [Cell::dirDiag]));
    }
    while (row || col)
    {
      char transformation = '\0';
      const Cell::Dir dir_new = cells [row] [col]. dirs [dir];
      switch (dir)
      {
        case Cell::dirUp:   transformation = '_'; 
                            deletions++;
                            ASSERT (row); 
                            row--; 
                            break;
        case Cell::dirLeft: transformation = '-';  
                            insertions++;
                            ASSERT (col); 
                            col--; 
                            break;
        case Cell::dirDiag: ASSERT (row); 
                            ASSERT (col);
                            {
                              const char c1 = pep1. seq [row - 1];
                              const char c2 = pep2. seq [col - 1];
                              if (   c1 == c2
                                  && ! isAmbigAa (c1)
                                  && ! isAmbigAa (c2)
                                 )
                              {
                                transformation = '|';
                                matches++;
                              }
                              else
                              {
                                transformation = ' ';
                                substitutions++;
                              }
                            }
                            row--; 
                            col--; 
                            break;
        default: ERROR;
      }
      ASSERT (transformation);
      transformations += transformation;
      dir = dir_new;
    }
  }
  reverse (transformations);  
#if 0
  PRINT (seq1);
  PRINT (seq2);
  PRINT (transformations);
  PRINT (pep2. seq. size ());
  PRINT (size2 ());
#endif
  ASSERT (pep1. seq. size () == size1 ());
  ASSERT (pep2. seq. size () == size2 ());
  
  
  stop1 = size1 ();
  stop2 = size2 ();
  if (semiglobalMatchLen_min)
  {
  	for (const char* c = & transformations [0]; *c == '-'; c++)
  	  start2++;
  	for (const char* c = & transformations [0]; *c == '_'; c++)
  	  start1++;
  	for (const char* c = & transformations [transformations. size () - 1]; stop2 && *c == '-'; c--)
  	  stop2--;
  	for (const char* c = & transformations [transformations. size () - 1]; stop1 && *c == '_'; c--)
  	  stop1--;
  }

		
	self_score1 = pep1. getSelfSimilarity (substMat, start1, stop1);
	self_score2 = pep2. getSelfSimilarity (substMat, start2, stop2);
}



#if 0
Align::Align (const Dna &dna1,
	            const Dna &dna2,
	            bool semiglobal,
	            size_t match_len_min)
{
  ASSERT (! dna1. sparse);
  ASSERT (! dna2. sparse);
  
	CNWAligner al;

#ifdef WU_BLASTN
	constexpr int match_score    =  5;
	constexpr int mismatch_score = -4;
	constexpr int gap_open       =   0;
	constexpr int gap_extent     = -10;
#else	
  // NCBI BLASTN
	constexpr int match_score    =  2;
	constexpr int mismatch_score = -3;
	constexpr int gap_open       = -5;
	constexpr int gap_extent     = -2;
#endif
  static_assert (match_score > 0, "match_score");
  static_assert (mismatch_score < 0, "mismatch_score");
  static_assert (gap_open <= 0, "gap_open");
  static_assert (gap_extent < 0, "gap_extent");	
  
  al. SetWm (match_score);
  al. SetWms (mismatch_score);
  al. SetWg (gap_open);
  al. SetWs (gap_extent);
  al. SetScoreMatrix (nullptr);
	
	{
		string seq1 (dna1. seq);
		string seq2 (dna2. seq);	
		strUpper (seq1);
		strUpper (seq2);
		al. SetSequences (seq1, seq2);
	}

	al. SetEndSpaceFree (semiglobal, semiglobal, semiglobal, semiglobal);
	score = al. Run ();		
  tr = al. GetTranscriptString();  	
	finish (dna1, dna2, semiglobal, match_len_min);
	
	self_score1 = match_score * (int) (stop1 - start1) /*dna1. seq. size ()*/;
	self_score2 = match_score * (int) (stop2 - start2) /*dna2. seq. size ()*/;
	ASSERT (self_score1 >= 0);
	ASSERT (self_score2 >= 0);
}



void Align::finish (const Seq &seq1,
	                  const Seq &seq2,
	                  bool semiglobal,
	                  size_t match_len_min)
{
  ASSERT ((bool) match_len_min == semiglobal);

  const string& s1 = seq1. seq;
  const string& s2 = seq2. seq;

  ASSERT (! s1. empty ());
  ASSERT (! s2. empty ());
  ASSERT (! tr. empty ());

  matches       = strCountSet (tr, "M");
  substitutions = strCountSet (tr, "R");
  insertions    = strCountSet (tr, "I");
  deletions     = strCountSet (tr, "D");
  ASSERT (matches + substitutions + insertions + deletions == tr. size ());
  // Global alignment
  ASSERT (s1. size () == size1 ());
  ASSERT (s2. size () == size2 ());
  
	stop1 = s1. size ();
	stop2 = s2. size ();
  if (semiglobal)
  {
  	for (const char* c = & tr [0]; *c == 'I'; c++)
  	  start2++;
  	for (const char* c = & tr [0]; *c == 'D'; c++)
  	  start1++;
  	for (const char* c = & tr [tr. size () - 1]; stop2 && *c == 'I'; c--)
  	  stop2--;
  	for (const char* c = & tr [tr. size () - 1]; stop1 && *c == 'D'; c--)
  	  stop1--;
  }
	ASSERT (start1 <= stop1);
	ASSERT (start2 <= stop2);
}
#endif



void Align::qc () const
{
  if (! qc_on)
    return;
    
  QC_ASSERT (transformations. size () == insertions + deletions + matches + substitutions);
    
	QC_ASSERT (start1 <= stop1);
	QC_ASSERT (start2 <= stop2);
  
  QC_ASSERT (self_score1 >= 0.0);
  QC_ASSERT (self_score2 >= 0.0);
  QC_ASSERT (score <= self_score1);
  QC_ASSERT (score <= self_score2);
}



AlignScore Align::getMinEditDistance () const
{	
  ASSERT (self_score1 >= 0.0);
  ASSERT (self_score2 >= 0.0);
  
  const AlignScore dist = self_score1 + self_score2 - 2 * score;
  ASSERT (dist >= 0.0);
	return dist;
}



void Align::printAlignment (const string &seq1,
	                          const string &seq2,
	                          size_t line_len) const
{
  ASSERT (! contains (seq1, '-'));
  ASSERT (! contains (seq2, '-'));
  ASSERT (line_len);

  string sparse1;    sparse1.   reserve (transformations. size ());
  string sparse2;    sparse2.   reserve (transformations. size ());
  string consensus;  consensus. reserve (transformations. size ());
  size_t i1 = 0;
  size_t i2 = 0;
  for (const char c : transformations)
  {
  	char c1 = '-';
  	char c2 = '-';
  	char cons = ' ';
  	switch (c)
  	{
  		case '|': cons = '|';      // Match
  		          c1 = seq1 [i1];  
  			        c2 = seq2 [i2];
  			        break;
  		case ' ': c1 = seq1 [i1];  // Replacement
  			        c2 = seq2 [i2];
  			        break;
  		case '-': c2 = seq2 [i2];  // Insertion
  			        break;
  		case '_': c1 = seq1 [i1];  // Deletion
  			        break;
  		default : ERROR;
  	}
    sparse1 += c1; 
  	sparse2 += c2; 
  	consensus += cons;
  	if (c1 != '-')
  		i1++;
  	if (c2 != '-')
  		i2++;
  }
  ASSERT (i1 == seq1. size ());
  ASSERT (i2 == seq2. size ());
  ASSERT (sparse1. size () == transformations. size ());
  ASSERT (sparse2. size () == transformations. size ());
  
  size_t qEnd = 0;
  size_t sEnd = 0;
  for (size_t i = 0; i < transformations. size (); i += line_len)
  {
    const string qStr (sparse1. substr (i, line_len));
    const string sStr (sparse2. substr (i, line_len));
    qEnd += (qStr. size () - strCountSet (qStr, "-"));
    sEnd += (sStr. size () - strCountSet (sStr, "-"));
    cout << qStr << ' ' << qEnd << endl;
    cout << consensus. substr (i, line_len) << endl;
    cout << sStr << ' ' << sEnd << endl;
    cout << endl;
  }
}




// Mutation

Mutation::Mutation (bool prot_arg,
                    const string& line)
: prot (prot_arg)
{ 
  string s (line);
  trim (s);
  
  const size_t dash_pos = s. find ('-');
  
  // geneName
  if (dash_pos != string::npos)
  {
    geneName = s. substr (0, dash_pos);
    QC_ASSERT (! geneName. empty ());
  }

  // ref, pos, allele
  const size_t refStart = dash_pos == string::npos ? 0 : (dash_pos + 1);
  size_t posStart = no_index;
  size_t alleleStart = no_index;
  FFOR_START (size_t, i, refStart, s. size ())
    if (isDigit (s [i]))
    {
      if (posStart == no_index)
        posStart = i;
    }
    else
    {
      if (posStart != no_index)
        if (alleleStart == no_index)
        {
          alleleStart = i;
          break;
        }
    }
  QC_ASSERT (refStart < posStart);
  QC_ASSERT (posStart < alleleStart);
  QC_ASSERT (alleleStart < s. size ());
  ref    =                s. substr (refStart, posStart - refStart);
  pos    = (size_t) stoi (s. substr (posStart, alleleStart - posStart));
  allele =                s. substr (alleleStart);
  QC_ASSERT (pos);
  pos--;

  // frameshift, ambig 
  if (prot)
  {
    if (ref == "ins")
      ref. clear ();
    if (allele == "del")
      allele. clear ();
    if (allele == "fs")
    {
      frameshift = true;
      allele. clear ();
    }
  #if 0
    for (const char c : allele)
      if (isAmbigAa (c))
        ambig = true;
  #endif
  }
  else
  {
    if (ref == "INS")
      ref. clear ();
    if (allele == "DEL")
      allele. clear ();
  #if 0
    for (const char c : allele)
      if (isAmbigNucl (c))
        ambig = true;
  #endif
  }
  
  setAmbig ();
}



Mutation::Mutation (string geneName_arg,
                    size_t pos_arg,
                    string ref_arg,
                    string allele_arg,
                    bool frameshift_arg)
: prot (true)
, geneName (geneName_arg)
, pos (pos_arg)
, ref (ref_arg)
, allele (allele_arg)
, frameshift (frameshift_arg)
{
  ASSERT (! geneName. empty ());
  setAmbig ();
}



void Mutation::setAmbig ()
{
  for (const char c : allele)
    if (isAmbig (c, prot))
    {
      ambig = true;
      return;
    }
}



void Mutation::qc () const
{
  if (! qc_on)
    return;
    
  QC_IMPLY (! geneName. empty (), isIdentifier (geneName, false));
  QC_ASSERT (pos != no_index);
  
  if (prot)
  {
    for (const char c : ref)
      if (c != *terminator && ! strchr (peptideAlphabet, c))
        throw runtime_error ("Protein mutation cannot have ambiguities in the reference sequence");
    for (const char c : allele)
      QC_ASSERT (strchr (extTermPeptideAlphabet, c));
  }
  else
  {
    for (const char c : ref)
      if (! strchr (dnaAlphabet, c))
        throw runtime_error ("DNA mutation cannot have ambiguities in the reference sequence");
    for (const char c : allele)
      QC_ASSERT (strchr (extDnaAlphabet, c));
  }

  QC_IMPLY (! frameshift, ref != allele);
  QC_IMPLY (frameshift, prot && allele. empty ());
}



bool Mutation::operator< (const Mutation& other) const
{ 
  LESS_PART (*this, other, prot);
  LESS_PART (*this, other, geneName);
  LESS_PART (*this, other, pos);
  LESS_PART (*this, other, ref);
  LESS_PART (*this, other, allele);
  LESS_PART (*this, other, frameshift);
  return false;
}



void Mutation::replace (Dna &refDna) const
{
  ASSERT (! prot);
  ASSERT (stop () <= refDna. seq. size ());
  ASSERT (refDna. seq. substr (pos, ref. size ()) == ref);
  refDna. seq. replace (pos, ref. size (), allele);
}




// KmerIndex

KmerIndex::IdRecord::IdRecord (fstream &f_arg,
                               Addr addr_arg,
                               bool isNew)
: f (f_arg)
, addr (addr_arg)
{ 
  ASSERT (addr != nil);
  
  if (isNew)
  {
    FOR (size_t, i, textSize)
      text [i] = no_char;
    ASSERT (! full ());
  }
  else
  {
    f. seekg ((streamoff) addr);
    readBin (f, prev);
    FOR (size_t, i, textSize)
      f. get (text [i]);
    start = getStart ();
    QC_ASSERT (start);
  }
}



bool KmerIndex::IdRecord::full () const
{
  ASSERT (start <= textSize);
  return start == textSize;
}



void KmerIndex::IdRecord::save () const
{         
  f. seekp ((streamoff) addr);
  writeBin (f, prev);
  FOR (size_t, i, textSize)
    f. put (text [i]);
}



void KmerIndex::IdRecord::renew (Addr addr_arg)
{
  ASSERT (addr_arg != nil);
  
  prev = addr;
  addr = addr_arg;
  ASSERT (prev < addr);
  
  FOR (size_t, i, textSize)
    text [i] = no_char;
    
  start = 0;
  
  ASSERT (! full ());
}



size_t KmerIndex::IdRecord::getStart () const
{ 
  size_t i = 0;
  for (; i < textSize; i++)
    if (text [i] == no_char)
      break;
  return i;
}



void KmerIndex::IdRecord::put (string &s)
{ 
  ASSERT (! s. empty ());
  ASSERT (! contains (s, no_char));
  
  size_t i = s. size ();
  while (i-- > 0)
  {
    if (full ())
    {
      s. erase (i + 1);
      return;
    }
    ASSERT (text [start] == no_char);
    text [start] = s [i];
    start++;
  }
  s. clear ();
}



string KmerIndex::IdRecord::get () const
{
  string s; s. reserve (textSize);
  FOR_REV (size_t, i, start)
    s += text [i];
  return s;
}



KmerIndex::KmerIndex (const string &name_arg,
                      size_t kmer_size_arg)
: Named (name_arg)
, f (name, ios_base::out | ios_base::binary)
, canRead (false)
, kmer_size (kmer_size_arg)
, code_max (powInt (2, 2 * kmer_size))
, addr_new (code2addr (code_max))
{
  QC_ASSERT (kmer_size);
  QC_ASSERT (kmer_size <= kmer_size_max);  
  ASSERT (code_max);
  
  writeBin (f, kmer_size);
  writeBin (f, items);
  ASSERT (f. tellp () == (streamoff) code2addr (0));
  
  {
    Progress prog (code_max, progressSize);  
    const Addr addr = nil;
    FOR (size_t, i, code_max)
    {
      prog ();
      writeBin (f, addr);
    }
  }
  
  QC_ASSERT (f. tellp () == (streamoff) addr_new);
  QC_ASSERT (f. good ());
}



KmerIndex::KmerIndex (const string &name_arg)
: Named (name_arg)
, f (name, ios_base::in | ios_base::out | ios_base::binary)
, canRead (true)
, kmer_size (readKmerSize (f))
, code_max (powInt (2, 2 * kmer_size))
{
  checkFile (name);

  readBin (f, items);

  f. seekg (0, ios_base::end);
  addr_new = (Addr) f. tellg ();

  QC_ASSERT (f. good ());
}



size_t KmerIndex::readKmerSize (fstream &fIn)
{ 
  size_t k;
  readBin (fIn, k);
  return k;
}



void KmerIndex::qc () const
{
  ASSERT (canRead);
  if (! qc_on)
    return;
    
  Named::qc ();
    
  QC_ASSERT (kmer_size);
  QC_ASSERT (kmer_size <= kmer_size_max);  
  QC_ASSERT (f. good ());
  QC_ASSERT (getFileSize (name) == (streamsize) addr_new);
  QC_ASSERT (addr_new < nil);
  // ??
}



size_t KmerIndex::dna2code (const Dna &dna)
{
  size_t code = 0;
  for (const char c : dna. seq)
  {
    code *= 4;
    size_t delta = 0;
    switch (c)
    {
      case 'a': delta = 0; break;
      case 'c': delta = 1; break;
      case 'g': delta = 2; break;
      case 't': delta = 3; break;
      default: throw runtime_error ("Unknown nucleotide: " + to_string (c));
    }
    code += delta;
  }
  return code;
}



size_t KmerIndex::getKmers ()
{
  size_t kmers = 0;
  f. seekg ((streamoff) code2addr (0));
  Progress prog (code_max, progressSize);  
  Addr addr = nil;
  FOR (size_t, i, code_max)
  {
    prog ();
    readBin (f, addr);
    if (addr != nil)
      kmers++;
  }
  return kmers;
}



void KmerIndex::add (const Dna &dna,
                     size_t &kmers,
                     size_t &kmersRejected)
{
  ASSERT (canRead);
  
  const size_t seqSize = dna. seq. size ();
  const string id (dna. getId ());
  kmers = 0;
  size_t kmersUsed = 0;
  FOR (size_t, i, seqSize)
    if (i + kmer_size <= seqSize)
    {
      kmers++;
      const Dna kmer ("x", dna. seq. substr (i, kmer_size), false);
      kmer. qc ();
      ASSERT (kmer. seq. size () == kmer_size);
      if (kmer. getXs ())
        continue;
      kmersUsed++;
      addId (dna2code (kmer), id);
    }
  ASSERT (kmers >= kmersUsed);
  kmersRejected = kmers - kmersUsed;
  
  items++;
  f. seekp (sizeof (kmer_size));
  writeBin (f, items);
}



void KmerIndex::addId (size_t code,
                       const string &id)
{
  ASSERT (canRead);
  ASSERT (code < code_max);
  QC_ASSERT (! id. empty ());
  QC_ASSERT (! contains (id, IdRecord::no_char));

  const streamoff pos = (streamoff) code2addr (code);
  
  Addr addr = 0;
  f. seekg (pos);
  readBin (f, addr);
  QC_ASSERT (addr == nil || (addr >= code2addr (code_max) && addr < addr_new));
  if (addr == nil)
    addr = addr_new;

  // s      
  static_assert (gap != IdRecord::no_char, "gap = no_char");
  string s (id);
  s += gap;
  
  IdRecord rec (f, addr, addr == addr_new);
  while (! s. empty ())
  {
    if (rec. full ())
    {
      ASSERT (rec. addr != addr_new);
      rec. renew (addr_new);
    }
    ASSERT (! rec. full ());
  #ifndef NDEBUG
    const size_t size_old = s. size ();
  #endif
    rec. put (s);
    ASSERT (s. size () < size_old);
    IMPLY (! s. empty (), rec. full ());
    rec. save ();
    if (rec. addr == addr_new)
    {
      addr_new += IdRecord::size;
      QC_ASSERT (addr_new < nil);
    }
  }
  
  f. seekp (pos);
  writeBin (f, rec. addr);

  QC_ASSERT (f. good ());
}



KmerIndex::NumId::NumId (size_t n_arg,
                         const string &id_arg)
: n (n_arg)
, id (id_arg)
{
  ASSERT (n);
  ASSERT (! id. empty ());
}



bool KmerIndex::NumId::operator< (const NumId &other) const
{ 
  LESS_PART (other, *this, n);
  LESS_PART (*this, other, id);  // Tie resolution
  return false;
}



Vector<KmerIndex::NumId> KmerIndex::find (const Dna &dna) 
{
  ASSERT (canRead);
  
  unordered_map<string,size_t> id2num;  id2num. rehash (100000);  // PAR
  const size_t seqSize = dna. seq. size ();
  FOR (size_t, i, seqSize)
    if (i + kmer_size <= seqSize)
    {
      const Dna kmer ("x", dna. seq. substr (i, kmer_size), false);
      kmer. qc ();
      ASSERT (kmer. seq. size () == kmer_size);
      if (kmer. getXs ())
        continue;
      const StringVector ids (code2ids (dna2code (kmer)));
      for (const string& id : ids)
        id2num [id] ++;
    }
    
  Vector<NumId> num2id;  num2id. reserve (id2num. size ());
  for (const auto& it : id2num)
    num2id << std::move (NumId (it. second, it. first));
  num2id. sort ();
  
  return num2id; 
}



StringVector KmerIndex::code2ids (size_t code)
{
  ASSERT (canRead);
  ASSERT (code < code_max);
  
  Addr addr = 0;
  {
    const streamoff pos = (streamoff) code2addr (code);  
    f. seekg (pos);
    readBin (f, addr);
  }
  QC_ASSERT (addr == nil || (addr >= code2addr (code_max) && addr < addr_new));
  
  string s;
  size_t n = 0;
  while (addr != nil)
  {
    const IdRecord rec (f, addr, false);
    s += rec. get ();
    addr = rec. prev;
    n++;
  }
  
  if (s. empty ())
    return StringVector ();
    
  ASSERT (s. back () == gap);
  s. erase (s. size () - 1, 1);

  return StringVector (s, gap, false);
}



}
