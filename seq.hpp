// seq.hpp

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


#ifndef SEQ_HPP
#define SEQ_HPP


#include "common.hpp"
using namespace Common_sp;




namespace Seq_sp
{


// Frame

typedef char Frame;
  // -3, -2, -1, 1, 2, 3
  // If > 0 then strand = 1 else strand = -1

inline bool isFrame (char frame)
  { return
      betweenEqual (frame, char (-3), char (-1)) ||
      betweenEqual (frame, char ( 1), char ( 3));
  }



// Strand

typedef char Strand;
  //  1 - top
  // -1 - bottom
  //  0 - unknown

inline bool isStrand (char strand)
  { return    strand == -1 
           || strand ==  1;
  }



// Alphabets

// DNA

// Lower case: gap is impossible
// Upper case: gap is possible

extern const char* dnaAlphabet;
  // Nucleotides:
  //   'a' - adenine
  //   'c' - cytosine
  //   'g' - guanine
  //   't' - thymine (RNA: 'u' - uracil)

extern const char* dnaWildcards;
/*
  IUPAC:
  m = ac
  r = ag  (purins)
  w = at
  s = cg
  y = ct  (pyramidins)
  k = gt  (ketons)
  v = acg
  h = act
  d = agt
  b = cgt
  n = acgt
*/
inline bool isAmbigNucl (char c)
  { return strchr (dnaWildcards, c); }


extern const char* extDnaAlphabet;
  // = DnaAlphabet + DnaWildcards

extern const char* extSparseDnaAlphabet;
  // extDnaAlphabet + gap

// '-': gap
// ' ': nothing


// Peptides
// Lower case

extern const char* peptideAlphabet;
  // Amino acides

extern const char* peptideWildcards;
inline bool isAmbigAa (char c)
  { return strchr (peptideWildcards, c); }


extern const char* extPeptideAlphabet;
  // = PeptideAlphabet + PeptideWildcards
/*
  X = any
  B = DN
  Z = EQ
  J = IL
  U = selenocysteine ("tga")
  O = Pyrrolysine ("tag") 
*/

extern const char* terminator;
  // strlen() = 1

extern const char* extTermPeptideAlphabet;
  // = ExtPeptideAlphabet + terminator

extern const char* extSparseTermPeptideAlphabet;
  // = extPeptideAlphabet + gap


size_t alphabet2Pos (const char* alphabet,
                     char        c);
  // Return: position of c in alphabet[]
  // Requires: c in alphabet[]
  
  
inline bool isAmbig (char c,
                     bool aa)
  { return aa ? isAmbigAa (c) : isAmbigNucl (c); }


typedef  unsigned char  Gencode;
  // NCBI genetic code
  
  


/////////////////////////// Seq //////////////////////////////////

// "Sparse sequence" is a sequence with '-'

inline size_t countInsertions (const string &seq)
  { return strCountSet (seq, "-"); }
  
inline size_t sparseSeqLen (const string &seq)
  { return seq. size () - countInsertions (seq); }

#if 0
int SparseSeq2SeqPos (const char* SparseSeq,
             		       int         SparseSeqPos);

int Seq2SparseSeqPos (const char* SparseSeq,
		                    int         SeqPos);
  // Return: SparseSeq [Result] != '-'
#endif



extern const size_t fastaLineLen;


typedef char (*GetIintersectChar) (const char* CharSet);
  // Input: CharSet[]
  // Return: Intersection character, ' ' if no intersection



struct Dna;
struct Peptide;



struct Seq : Named
// name: name of sequence
// Positions: 0-based
//            stop: position of the last character + 1
//              start != pos
//              size = stop - start
//            Reverse strand <=> start > stop
{
  string seq;
    // in getSeqAlphabet()
    // May be empty()
  bool sparse {false};
    // true <=> '-' is allowed


protected:
  Seq () = default;
  Seq (const string &name_arg,
       size_t seqLen,
       bool sparse_arg)
    : Named (name_arg)
    , seq (seqLen, '\0')
    , sparse (sparse_arg)
    { qcName (); }
  Seq (const string &name_arg,
       const string &seq_arg,
       bool sparse_arg)
		: Named (name_arg)
		, seq (seq_arg)
    , sparse (sparse_arg)
		{ qcName ();
			if (! sparse)
			  unSparse ();
	  }
  Seq (LineInput &fasta,
       size_t reserveLen,
       bool sparse_arg,
     	 bool makeUpper);
    // Reads 1 sequence and skips blank lines
    // Requires: fasta.line is the first line of the sequence
    // If !makeUpper then make lowercase
	void qcName () const;
	void qcAlphabet () const;
	  // Invokes: getSeqAlphabet()
  static LineInput* getLineInput (const string &fName)
    { LineInput* in = new LineInput (fName);
      in->nextLine ();
      return in;
    }
public:
  virtual Seq* copy () const override = 0;
  void qc () const override
    { if (qc_on)
        return;
      qcName ();
    	qcAlphabet ();
    }


  virtual const Dna* asDna () const
    { return nullptr; }
  virtual const Peptide* asPeptide () const
    { return nullptr; }

  virtual const char* getSeqAlphabet () const = 0;
  virtual bool isAmbiguous (char c) const = 0;
  void printCase (ostream &os,
           	      bool makeUpper) const;    
    // FASTA format
    // Length of lines = fastaLineLen

  size_t getIdSize () const
    { const size_t idSize = name. find (' ');
	    return idSize == string::npos ? name. size () : idSize; 
	  }
	static size_t getTaxonStart (const string &s);
	  // Return: != string::npos => s[Return] == '['
	size_t getTaxonStart () const
	  { return getTaxonStart (name); }
  string getId () const
    { return name. substr (0, getIdSize ()); }
  long /*CSeq_id::TGi*/ getGi () const;
    // If there is no gi then throw
  string getDescription (bool trimTaxon) const;
    // Return: between Id and taxon
  size_t getXs () const;
    // Invokes: isAmbiguous()
  size_t getContiguousXs () const;
    // Return: max. length of contiguous ambiguous characters
  virtual double getComplexityInt (size_t start,
                                   size_t end) const = 0;
    // Input: Subsequence form start to end
    // Return: Entropy of diletters divided by log_2 e; >= 0
    //         If > stdMinComplexity then no repeats
  double getComplexity () const 
    { return getComplexityInt (0, seq. size ()); }
#if 0
  size_t* GetAlphabetCount () const;
    // Return: length = strlen (SeqAlphabet)
    // Invokes: NewUintArray ()
#endif
  void unSparse ()
    { strDeleteSet (seq, "-"); 
      sparse = false;
    }
  map<char,size_t> getCharCount () const;
};



struct Multifasta : Root
/* Usage:
     { Multifasta fa (...);
       while (fa. next ())
         {Peptide|Dna} ... (fa); 
     }  // To close LineInput and Progress
*/
{
  LineInput in;
  const bool aa;
  Progress prog;

	Multifasta (const string &fName,
	            bool aa_arg,
	            size_t displayPeriod = 1000)  
		: in (fName) 
		, aa (aa_arg)
		, prog (0, displayPeriod)  
		{ in. nextLine (); 
			qcNewSeq ();
		}
private:
	void qcNewSeq () const;
public:
	  
	bool next () const
		{ qcNewSeq ();
			return ! in. eof;
		}
};




////////////////////////////// Dna /////////////////////////////////

// acgtb: 4 nucleotides and blank
// Blank is encoded by '-'
// One of them must be true


size_t nuc2num (char wildNucleotide);
  // Return: if wildNucleotide is in dnaAlphabet then position in dnaAlphabet
  //         else if wildNucleotide is in dnaWildcards then 4
  //         else error
  // Case-insensitive

#if 0
typedef PROBABILITY NT_PROBABILITY [5];


void PrintNTProb (const NT_PROBABILITY Prob);


void CompressDna (const char* seq,
                  char*       CompressedSeq,
                  size_t        &CompressedSeqLen);
  // Input: seq in extDnaAlphabet + '-' + uppercase
  // Output: CompressedSeq []
  //         CompressedSeqLen
  // Requires: CompressedSeq [] must have enough space (<= strlen (seq) + 1)

void UncompressDna (size_t        SeqLen,
                    const char* CompressedSeq,
                    char*       &seq,
                    size_t        &CompressedSeqLen);
  // Output: *seq, malloc (), length = SeqLen + 1
  //         CompressedSeqLen
#endif

uchar wild2nucleotides (char wildNucleotide,
 		                    bool acgtb [5]);
  // Output: acgtb []
  // Return: Number of 1's in acgtb []; 0..5
  //         |WildNucleotide|, "Degree of ambiguity"
  // Requires: WildNucleotide in extDnaAlphabet + '-' + ' '

char nucleotides2wild (const bool acgtb [5]);
  // Return: in extDnaAlphabet + '-'; if empty set then ' '

#if 0
void Wild2NucleotideFreq (char           WildNucleotide,
                          NT_PROBABILITY acgtb);
  // Output: acgtb []
  //         Sum_i acgtb [i] = 1
#endif

char complementaryNucleotide (char wildNucleotide);
  // Return: in extDnaAlphabet
  // Possible blank is preserved
  // Requires: wildNucleotide in extDnaAlphabet

char getUnionNucleotide (const string& charSet);
  // Return: in extDnaAlphabet + '-'; if empty set then ' '
  // Requires: CharSet is in extDnaAlphabet + '-'

char getIntersectNucleotide (const char* charSet);
  // GetIintersectChar
  // Return: in extDnaAlphabet + '-'; if empty set then ' '
  // Requires: CharSet is in extDnaAlphabet + '-'

bool nucleotideMatch (char wildNucleotide1,
                      char wildNucleotide2);
  // Return: true if the nucleotides of wildNucleotide1 and wildNucleotide2 intersect

#if 0
bool NucleotideSeqMatch (const char* Seq1,
                         const char* Seq2);
  // Requires: Seq1|2 [] - in extDnaAlphabet

bool MoreGeneralNucleotide (char WildNucleotide1,
                            char WildNucleotide2);
  // Return: true if the nucleotides of WildNucleotide1 are a superset
  //           of those of WildNucleotide2 (i.e., WildNucleotide1 is
  //           more general than or equal to WildNucleotide2)

void SelectSpecificNucleotides (char* CharSet);
  // Delete more general nucleotides
  // Update: CharSet []
  // Requires: CharSet is in extDnaAlphabet + '-'

size_t CountAmbiguousNucleotides (const char* seq);
  // Input: seq - sparse sequence

char* getReverseDna (const char* Source);
  // Return: malloc ()
  // Requires: source is in extDnaAlphabet + '-'
#endif

string& reverseDna (string &seq);
  // Requires: seq is in extDnaAlphabet + '-'

char codon2aa (const char codon [3],
               Gencode gencode,
               bool lowercasePossibleStartCodon);
  // Return: in extTermPeptideAlphabet; lowercasePossibleStartcodon => lower-case
  // Requires: codon[] is in extDnaAlphabet

inline size_t dna2codons_len (size_t dna_len)
  { return dna_len ? (dna_len - 1) / 3 + 1 : 0; }

#if 0
inline bool CodonMatch (const char Codon1 [3],
           		           const char Codon2 [3])
  { return NucleotideSeqMatch (Codon1, Codon2); }
  // Requires: Codon1|2 [] - in extDnaAlphabet

inline bool TerminatorCodon (const char Codon [3])
  { return Codon2AA (Codon) == '*'; }

bool MayBeTerminatorCodon (const char Codon [3]);


size_t Dna2PeptidePos (size_t  DnaPos,
              		     Frame frame,
              		     size_t  DnaLen);

size_t Peptide2DnaPos (size_t  PeptidePos,
              		     Frame frame,
              		     size_t  DnaLen);
#endif



struct Dna : Seq
{
  static constexpr double stdMinComplexity {2.0};  // PAR


  Dna () = default;
  Dna (const string &name_arg,
       size_t seqLen,
       bool sparse_arg)
    : Seq (name_arg, seqLen, sparse_arg)
    {}
  Dna (const string &name_arg,
       const string &seq_arg,
       bool sparse_arg)
    : Seq (name_arg, seq_arg, sparse_arg)
    {}
  Dna (LineInput &fasta,
       size_t reserveLen,
       bool sparse_arg)
    : Seq (fasta, reserveLen, sparse_arg, false)  
    {}
  Dna (Multifasta &fasta,
       size_t reserveLen,
       bool sparse_arg);
  Dna* copy () const final
    { return new Dna (*this); }
	void saveText (ostream& os) const override
	  { printCase (os, false); }


  const Dna* asDna () const final
    { return this; }

  const char* getSeqAlphabet () const final
    { return sparse ? extSparseDnaAlphabet : extDnaAlphabet; }
  bool isAmbiguous (char c) const final
    { return ! strchr (dnaAlphabet, c); }
  double getComplexityInt (size_t start, 
                           size_t end) const final;
#if 0
  void PrintHTML (bool        UpperCase,
                  PHRED_SCORE MinGoodQual) const;
  SEQ_ANNOT_LIST* GetAnnotList (int GoodStart = 0) const;
    // List: "seq", "Score 10", "Score  1"

  // Qual
  void CreateQual ();
    // Requires: Qual == nullptr
    // Postcondition: Qual != nullptr
  bool GoodQual () const;
    // Return: false if there are wrong values in Qual []
  void CopyQual (const PHRED_SCORE* SourceQual);
    // Invokes: CreateQual (), GoodQual ()
  void ReadQual (const char* FName);
    // Input: FName: FASTA-file with phred quality scores
    // Output: Qual
    // Invokes: CreateQual (), GoodQual ()
  // Requires: Qual != nullptr
  void SetQual (PHRED_SCORE DefaultScore);
  PHRED_SCORE GetMinQual () const;
  PHRED_SCORE GetMaxQual () const;
  void Qual2MeanVar (size_t  Start,
                     size_t  End,
                     double &Mean,
                     double &Var) const;
    // Output: Mean, Var (unbiased)
  void Qual2N (PHRED_SCORE MinQual,
               bool        GapCoded);
    // If Qual [i] < MinQual then seq [i] := 'N'|'n'
    // Delete Qual
  void MinScore2Pos (PHRED_SCORE MinScore,
                     size_t        &Start,
                     size_t        &End) const;
    // Output: Start, End; valid if Start < End
  void QualSaveFile (FILE*       F,
                  		 PHRED_SCORE DefaultScore) const;
  char* Qual2String () const;
    // Return: new [], May be nullptr

  void DeleteStart (size_t NewStart);
#endif

  Dna* makeComplementary () const;
    // Return: Reverse and complementary Dna; !nullptr
  void reverse ();
    // Invokes: makeComplementary ()

  Peptide makePeptide (Frame frame,
                       Gencode gencode,
                       bool lowercaseStartcodon,
                       bool firstStartCodon2M,
                       size_t &translationStart) const;
    // Output: translationStart
    // Requires: firstStartCodon2M => lowercaseStartcodon
  Peptide cds2prot (Gencode gencode,
                    bool trunc5,
                    bool trunc3,
	                  bool hasStopCodon,
	                  bool allowExtraStopCodon) const;
	  // Return: if !trunc3 && hasStopCodon then no '*'
	  // Invokes: makePeptide(1,gencode,true,true)
  Vector<Peptide> getOrfs (Frame frame,
                           Gencode gencode,
                           size_t len_min) const;
    // Input: len_min: min. Peptide length without 'X'
    // Return: Peptide: from '*' to '*'
#if 0
  bool ExistsPeptide () const;
    // Return: true if there is peptide w/o an internal stop codon in some frame
  bool LongestCompleteCDS (size_t MinProteinLen,
                           size_t &Start,
                           size_t &End) const;
    // Output: Start, End - position after last nt;
    //         (End - Start) % 3 = 0
    //         End <= strlen (seq)
    //         Valid if true
    // Return: true if found

  void RefineSparse (char* Target,
                     byte* TargetQual) const;
    // Target[]='n' or seq[] is ambiguous
    // Update: Target []: sparse sequence
    //         Qual []; if Target [i] = '-' then TargetQual [i] = 0

  bool ContainsAmbiguity () const;
  size_t GetAmbiguousPrefixEnd () const;
    // = length of ambiguous nt prefix
  size_t GetAmbiguousSuffixStart () const;
  bool DeleteAmbiguousPrefix ();
    // Return: True if non-empty prefix is deleted
  bool DeleteAmbiguousSuffix ();
    // Return: True if non-empty suffix is trimmed

  void TrimN (bool GapCoded);  
  size_t TrimBadStart (size_t        WindowLen,
		                   PROBABILITY MaxAmbigFraction) const;
    // Return: Good sequence start
  size_t TrimBadEnd (size_t        WindowLen,
	                  PROBABILITY MaxAmbigFraction) const;
    // Return: Good sequence end
  size_t TrimLowComplexityStart () const;
    // Return: Good sequence start
  size_t TrimLowComplexityEnd () const;
    // Return: Good sequence end
#endif

  // Poly-nucleotide segment at the end of seq[]
  bool polyNucWindow (size_t start,
                      size_t windowLen,
              		    char nucleotide) const;
    // Return: true if most of the nucleotides in the window are nucleotide
  size_t findPolyNucEnd (size_t windowLen,
                         char nucleotide) const;
    // Find poly-nucleotide segment at the end of seq
    // Return: start of the poly-nucleotide segment; seq.size() if not found
    // Invokes: polyNucWindow ()
  // PolyA
  bool polyAWindow (size_t start,
           		      size_t windowLen) const
    { return polyNucWindow (start, windowLen, 'a'); }
  size_t findPolyA (size_t windowLen) const
    { return findPolyNucEnd (windowLen, 'a'); }
  size_t removePolyA (size_t windowLen)
    { const size_t polyA_start = findPolyA (windowLen);
      const size_t polyA_len = seq. size () - polyA_start;
      seq. erase (polyA_start);
      return polyA_len;
    }
    // Non-idempotent
  
#if 0
  // PolyA inside seq[]
  void FindPolyAWindow (int   &BestStart,
                        int   &BestEnd,
                        double &MaxWeight) const;
    // Output: BestStart = -1 iff MaxWeight < 0

  bool DiNucLowComplexity (int   &BestStart,
                           int   &BestEnd,
                           double &MaxWeight,
                           int   &MaxMonoNucLen) const;
    // Find a low-complexity segement by di-nucleotide statistics
    // Output: Valid if true
    //         MaxWeight: >= 20: significant
    //         MaxMonoNucLen: inside [BestStart, BestEnd]
    // Return: true if found, MaxWeight >= 0

  size_t NucleotideFreq (NT_PROBABILITY acgtb,
                       size_t           Start,
                       size_t           End,
                       uint           Step) const;
    // Return: number of nucleotides couned
    // Output: acgtb []; sum = 1
  size_t NucleotidePairFreq (PROBABILITY acgtb [5] [5]) const;
    // Input: Start, End -??
    // Return: strlen (seq)
    // Output: acgtb [] []; sum = 1

  double GetMeltTemp (double DnaConc  = 300.0,
                      double SaltConc = 50.0,
                      double MgConc   = 1.5) const;
    // Thermodynamic (base-stacking calculated) T_m (melting temperature)
    // Input: DnaConc is the primer concentration, nM
    //        SaltConc is monovalent cation (K+ and Na+) concentration, mM
    //        MgConc, mM
    // Return: The melting temperature, in C^0
    // Requires: strlen (seq) is between 7 and 35
    // From: http://www-genome.wi.mit.edu/ftp/distribution/software/primer.0.5/primer.c, 
    //       now Primer3 (http://frodo.wi.mit.edu/primer3/primer3_code.html)
    // VectorNTI makes the same results
  PRIMER* MakePrimer (size_t        PrimerLen,
                      size_t        SearchStart,
                      size_t        SearchEnd,
                      PHRED_SCORE MinQuality,
                      uint        MaxMaxRepeat,
                      PROBABILITY MaxGCShare,
                      double       DnaConcentration,
                      double       SaltConcentration,
                      double       MgConcentration,
                      double       MinMeltTemp,
                      double       MaxMeltTemp,
                      double       TargetMeltTemp) const;
    // From 5' to 3'
    // Return: Best primer, May be nullptr
    //         Primer inside the segment [SearchStart, SearchEnd]
    // Requires: TargetMeltTemp between MinMeltTemp and MaxMeltTemp
    // Invokes: GetMeltTemp ()

  PROBABILITY NucleotideFreqFrame (size_t        Start,
                     			         Frame       frame,
                     			         size_t        Count,
                     			         const char* Alphabet) const;
  void GetOpenFrameProb (size_t        Start,
		                     size_t        Count,
                     		 PROBABILITY FrameProb [3]) const;
    // Output: FrameProb []

  // Author: Max Troukhan
  bool GetOrf (int                  &CDS_Start,
               bool                 &FL5,
               int                  &CDS_End,
               bool                 &FL3,
               const CodStatistics* StatTable,
               double                Stress) const;
    // Output: CDS_Start|End - may be outside seq
    //         CDS_End - next position after the last known codon of the CDS
    //         FL5|3
    //         Valid if result is true
    //         CDS_Start < CDS_End
    // Return: true if Orf is found
  double CompareOrfs (const CodStatistics* StatTable1,
                      const CodStatistics* StatTable2) const;
    // Return: Positive iff seq is more probable to be of organism of StatTable1
    //         CompareOrfs (a,b) = - CompareOrfs (b,a)
  double Orf2Strand (const CodStatistics* StatTable) const;
    // Return: Positive iff the positive strand of seq is more probable than the negative strand
    //         Orf2Strand (a) = - Orf2Strand (Reverse(a))


  // Input: Neightborhood: to search start/stop codon, in nts
  // Output: if true then CDS_Start < CDS_End and (CDS_End - CDS_Start) % 3 = 0
  // Return: true iff FL
  bool CDS_Start_nr2Orf (int CDS_Start_nr,
                         int Neighborhood,
                         int &CDS_Start,
                         int &CDS_End) const;
  bool CDS_End_nr2Orf (int CDS_End_nr,
                       int Neighborhood,
                       int &CDS_Start,
                       int &CDS_End) const;
    // Output: CDS_Start: minimum possible
#endif
};



struct FastaDna : Dna
{
  FastaDna (const string &fName,
            size_t reserveLen,
            bool sparse_arg)
    : Dna (* unique_ptr<LineInput> (getLineInput (fName)), reserveLen, sparse_arg)
    {}
};



#if 0
class DNA_COLLECTION: public _SEQ_COLLECTION
// of Dna*
{
typedef _SEQ_COLLECTION inherited;

public:
  DNA_COLLECTION (): inherited () {}
  DNA_COLLECTION (const char* FName,
                  bool        GapCoded);


  Dna* GetDna (size_t i) const
    { return (Dna*) GetSeq (i); }

  void SaveQual (const char* FName,
                 uint        DefaultQualScore) const;
};
#endif




/////////////////////////// Peptide ///////////////////////////////

size_t aa2num (char wildAminoacid);
  // Return: if wildAminoacid is in peptideAlphabet then position in peptideAlphabet
  //         else if wildAminoacid is terminator then 20
  //         else if wildAminoacid is in peptideWildcards then 21
  //         else error
  // Case-insensitive

bool moreGeneralAminoacid (char wildAminoacid1,
                           char wildAminoacid2);
  // Return: true if the aminoacids of wildAminoacid1 is a superset of those of wildAminoacid2 
  //           (i.e., wildAminoacid1 is more general than or equal to wildAminoacid2)
  // Requires: wildAminoacid1|2 in extTermPeptideAlphabet

inline bool aaMatch (char aa1,
                     char aa2)
  { return    moreGeneralAminoacid (aa1, aa2)
           || moreGeneralAminoacid (aa2, aa1);
  }

typedef  double  AlignScore;



struct SubstMat : Root
// Substitution matrix for protein alignment
// NCBI: e.g., /am/ftp-blast/matrices/BLOSUM62
/*           ------ default -----
  matrix     gap_open  gap_extent   
  --------   --------  ----------   
	BLOSUM90         10           1
  BLOSUM80         10           1
  BLOSUM62         11           1 
  BLOSUM50         13           2
	BLOSUM45         14           2
  PAM250           14           2
  PAM70            10           1
  PAM30             9           1
  IDENTITY         15           2
*/
{
  static constexpr const char* chars {"ARNDCQEGHILKMFPSTWYVBZX*"};
  static constexpr size_t sim_size = 128;
  AlignScore sim [sim_size] [sim_size];  
  
  
  explicit SubstMat (const string &fName);
  void qc () const override;
    // Print anomalies
  void saveText (ostream& os) const override;
    
    
  bool goodChar (size_t i) const
    { return sim [i] [i] == sim [i] [i]; }  // Not a NaN
  AlignScore getSubstitutionDist (size_t row,
                                  size_t col) const
    { return sim [row] [row] + sim [col] [col] - 2 * sim [row] [col]; }
  AlignScore getDeletionDist (size_t row,
                              AlignScore gap_sim) const
    { return sim [row] [row] - 2 * gap_sim; }
};

	
	
struct PeptideOrf : Root
// DNA -> translation of a segment -> ORF
{
  size_t translationStart {no_index};
    // Start of the peptide translation
    // !strand => !translationStart
  Strand strand {0};
  // In the translation which starts from translationStart
  size_t start {no_index};
    // May be: no_index
  bool startM {false};
    // Peptide::isStartAa(peptide->seq[start])
  size_t stop {no_index};
    // >= start
  bool stopTerminator {false};
    // peptide->seq[stop] == *terminator


  PeptideOrf (size_t translationStart_arg,
              Strand strand_arg,
      		     const Peptide* peptide,
             size_t start_arg);
    // Input: peptide: DNA segment translation
  PeptideOrf (size_t translationStart_arg,
              Strand strand_arg,
      		    size_t start_arg,
      		    bool startM_arg,
      		    size_t stop_arg,
      		    bool stopTerminator_arg)
		: translationStart (translationStart_arg)
		, strand (strand_arg)
		, start (start_arg)
		, startM (startM_arg)
		, stop (stop_arg)
		, stopTerminator (stopTerminator_arg)
		{}
		// In peptide translated from the DNA
  explicit PeptideOrf (istream &is)
    { int iStrand;
    	is >> translationStart >> iStrand >> start >> stop >> startM >> stopTerminator; 
    	strand = (char) iStrand;
    }
  PeptideOrf () = default;
  void qc () const override;
	void saveText (ostream &os) const override
	  { os         << translationStart 
	  	   << '\t' << (int) strand 
	  	   << '\t' << start 
	  	   << '\t' << stop 
	  	   << '\t' << startM
	  	   << '\t' << stopTerminator;
	  }
	  // Requires: !empty()
    
    
  bool empty () const override
    { return start == no_index; }
  bool good (size_t size_min) const
    { return    ! empty ()
    	       && startM        
		         && stopTerminator
		         && size () >= size_min;
    }
  // Requires; !empty()
  size_t size () const
    { return stop - start; }
    // Without the terminator
  size_t dnaPos (size_t pos) const
    { const size_t dnaLen = 3 * pos;
    	return strand == 1 ? (translationStart + dnaLen) : (translationStart - dnaLen); 
    }
    // Requires: (bool)frame 
  size_t cdsStart () const
    { return dnaPos (start); }
  size_t cdsStop () const
    { return dnaPos (stop + stopTerminator); }
  Peptide* toPeptide (const Peptide* peptide) const;
    // Input: peptide: used in Peptide::getOrfs()
};



struct Peptide : Seq
{
  static constexpr const size_t stdAveLen {400};  // PAR
  static constexpr double stdMinComplexity {2.5};  // PAR
  bool pseudo {false};


  Peptide () = default;
  Peptide (const string &name_arg,
           size_t seqLen,
           bool sparse_arg)
    : Seq (name_arg, seqLen, sparse_arg) 
    {}
  Peptide (const string &name_arg,
       	   const string &seq_arg,
       	   bool sparse_arg)
    : Seq (name_arg, seq_arg, sparse_arg) 
    {}
  Peptide (LineInput &fasta,
           size_t reserveLen,
           bool sparse_arg)
    : Seq (fasta, reserveLen, sparse_arg, true) 
    {}
  Peptide (Multifasta &fasta,
           size_t reserveLen,
           bool sparse_arg);
  Peptide* copy () const final
    { return new Peptide (*this); }
  void qc () const override;
	void saveText (ostream& os) const override
	  { printCase (os, true); }


  const Peptide* asPeptide () const final
    { return this; }

  const char* getSeqAlphabet () const final
    { return sparse ? extSparseTermPeptideAlphabet : extTermPeptideAlphabet; }
  bool isAmbiguous (char c) const final
    { return c == 'X'; }
  bool hasInsideStop () const
    { const size_t pos = seq. find ('*');
	    return pos != string::npos && pos != seq. size () - 1;
    }
  bool trimStop ()
    { return trimSuffix (seq, "*"); }
  bool isDescriptionPartial () const;
  static bool isStartAa (char aa) 
    { return aa != '*' && isLower (aa); }
    // Requires: after codon2aa(,,true)
  Vector<PeptideOrf> getPeptideOrfs (size_t translationStart,
                                     Strand strand,
                                     bool includeInitial,
                                     bool longestOnly,
                                     size_t len_min) const;
  double getComplexityInt (size_t start,
                           size_t end) const final;
  // Return: >= 0; -1 <=> seq has an amino acid not in *peptideAlphabet
  double getSelfSimilarity (const SubstMat &mat,
                            size_t start = 0,
                            size_t stop = 0) const;
    // Input: start, stop = 0 => whole sequence
  double getSimilarity (const Peptide &other,
                        const SubstMat &mat,
                        double gapOpenCost,
                        double gapCost) const;
    // Requires: *this and other are aligned
#if 0
  size_t GetLeftMPos (size_t Start) const;
    // Return: may be no_index
  size_t GetClosestMPos (size_t Start) const;
    // Return: may be no_index
#endif
  size_t ambig2X ();
    // Return: # ambiguous characters converted to 'X'
  void toGBMR4 ();
    // Reduced alphabet
};



struct FastaPeptide : Peptide
{
  FastaPeptide (const string &fName,
                size_t reserveLen,
                bool sparse_arg)
    : Peptide (* unique_ptr<LineInput> (getLineInput (fName)), reserveLen, sparse_arg)
    {}
};



inline Seq* makeSeq (const string &name,
				       	     const string &seq,
				       	     bool aa,
				       	     bool sparse)
  { Seq* s = nullptr;
  	if (aa)
  		s = new Peptide (name, seq, sparse);
  	else
  		s = new Dna (name, seq, sparse); 
  	s->qc ();
  	return s;
  }

inline Seq* makeSeq (Multifasta &fa,
                     bool sparse)
  { Seq* s = nullptr;
  	if (fa. aa)
  		s = new Peptide (fa, Peptide::stdAveLen, sparse);
  	else
  		s = new Dna (fa, 500/*PAR*/, sparse); 
  	s->qc ();
  	return s;
  }




//

struct Cds : Root
// All lengths are in na
{
  // PAR
  static constexpr size_t peptideSize_min = 20;  // aa
  static constexpr size_t promoter_min = 60;
  //
  static constexpr size_t size_min = 3 * (peptideSize_min + 1/*stop codon*/);

	// Key
	size_t start {0};  
	size_t stop {0};  
	string refProt;
	double positives {0.0};
	  // 0..1

  // Output	  
	Vector<const Cds*> prevs;
	  // != nullptr
	const Cds* bestPrev {nullptr};
	long sumLen {0};
	  // To be maximized
	
	
  Cds (size_t start_arg,
		   size_t stop_arg)
    : start (start_arg)
    , stop (stop_arg)
    , positives (1.0)
    {}
	Cds (size_t start_arg,
			 size_t stop_arg,
			 const string &refProt_arg,
			 double positives_arg)
    : start (start_arg)
    , stop (stop_arg)
    , refProt (refProt_arg)
    , positives (positives_arg)
    {}
	Cds () = default;
  void qc () const override;
  void saveText (ostream &os) const override
    { os << left () << '\t' << right () << '\t' << (int) strand () << '\t' << ((size () - 1) / 3); }


  bool strand () const
    { return start < stop; }
  size_t left () const
    { return min (start, stop); }
  size_t right () const
    { return max (start, stop); }
  size_t start_human () const
    { return (strand () ? start : stop) + 1; }
  size_t stop_human () const
    { return strand () ? stop : start; }
  size_t size () const
    { return right () - left (); }
  int sizeEffective () const
    { return (int) size () - ((int) size_min - 3/*PAR*/); }
    // Return: > 0
  size_t getOverlap (const Cds &other) const
    // Symmetric
    { if (other. left () >= right ())
    		return 0;
    	if (left () >= other. right ())
    		return 0;
    	if (left () <= other. left ())
    		return getOverlap_ (other);
    	return other. getOverlap_ (*this);
    }
private:
  size_t getOverlap_ (const Cds &other) const;
public:
  bool coexists (const Cds &next) const;
    // Symmetric
  int getLenIncrease (const Cds &prev) const;
  bool operator< (const Cds &other) const;
  bool worse (const Cds &than) const
    { return    than. left () == left ()
    	       && than. size () == size ()
    	       && (   than. positives > positives
    	           || (   than. positives == positives
    	               && than. refProt < refProt  // tie
    	              )
    	          );
    }
};



struct DnaAnnot : Root
{
  Vector<Cds> cdss;  
  
  
  DnaAnnot ()
    { cdss. reserve (10000); } // PAR
  
  
  const Cds* run ();    
    // Return: last best Cds; may be nullptr
    //         to get all best Cds's: loop by cds->bestPrev
    // Input: cdss
    
};




// Align

struct Align : Root
// Needleman-Wunsch algorithm 
{
	// Output:
	AlignScore score {0.0};
	  // To be maximized
	string transformations;
	  // Of sequence 1
	  // '-': insertion
	  // '_': deletion
	  // '|': match
	  // ' ': substitution
	size_t matches {0};
	  // = identities
	  // Not ambiguities
	size_t substitutions {0};
	size_t insertions {0};
		// In seq1
	size_t deletions {0};
		// = insertions in seq2
	
	size_t start1 {0};
	size_t start2 {0};
  size_t stop1 {0};
	size_t stop2 {0};
	
	AlignScore self_score1 {0.0};
	AlignScore self_score2 {0.0};
	
	
  // Default NCBI alignment parameters
  // Input: match_len_min: valid if semiglobal, otherwise 0
  // !semiglobal = global
	Align (const Peptide &pep1,
         const Peptide &pep2,
         const SubstMat &substMat,
         AlignScore gapOpenCost,
         AlignScore gapCost,
         size_t semiglobalMatchLen_min);
    // Global alignment <=> !semiglobalMatchLen_min
#if 0
  ??
	Align (const Dna &dna1,
	       const Dna &dna2,
	       bool semiglobal,
	       size_t match_len_min);
private:
	void finish (const Seq &seq1,
	             const Seq &seq2,
	             bool semiglobal,
	             size_t match_len_min);
public:
#endif
  void qc () const override;
	
	
	// Sizes of the input sequences
  size_t size1 () const
    { return transformations. size () - insertions; }
  size_t size2 () const
    { return transformations. size () - deletions; }
  //
	AlignScore getMinEditDistance () const;
    // Return: >= 0
  void printAlignment (const string &seq1,
	                     const string &seq2,
	                     size_t line_len) const;
};



struct Mutation : Root
{
  bool prot {false};
  string geneName;
    // May be empty()
  // !prot => positive strand
  size_t pos {no_index};
    // In reference
    // != no_index
  string ref;
  string allele;
  bool frameshift {false};
    // => prot
  
  bool ambig {false};
    // Function of allele
    

  Mutation (bool prot_arg,
            const string& line);
  Mutation (string geneName_arg,
            size_t pos_arg,
            string ref_arg,
            string allele_arg)
    : prot (! geneName_arg. empty ())
    , geneName (geneName_arg)
    , pos (pos_arg)
    , ref (ref_arg)
    , allele (allele_arg)
    { setAmbig (); }
  Mutation (string geneName_arg,
            size_t pos_arg,
            string ref_arg,
            string allele_arg,
            bool frameshift_arg);
private:
  void setAmbig ();
public:
  Mutation () = default;
  Mutation (const Mutation&) = default;
  Mutation& operator= (const Mutation&) = default;
  Mutation (Mutation&&) = default;
  Mutation& operator= (Mutation&&) = default;
  void qc () const override;
  void saveText (ostream &os) const override
    { if (prot)
        os << geneName << '-' << nvl (ref, "ins") << pos + 1 << (frameshift ? "fs" : nvl (allele, "del")); 
      else
        os <<                    nvl (ref, "INS") << pos + 1 <<                      nvl (allele, "DEL"); 
    }
    

  bool operator== (const Mutation& other) const
    { return    prot       == other. prot
             && geneName   == other. geneName
             && pos        == other. pos
             && ref        == other. ref
             && allele     == other. allele
             && frameshift == other. frameshift;
    }
  bool operator< (const Mutation& other) const;
  struct Hash
  { size_t operator() (const Mutation &mut) const
      { static hash<string> strHash;
        return   (size_t) mut. prot
               ^ strHash (mut. geneName) 
               ^          mut. pos 
               ^ strHash (mut. ref) 
               ^ strHash (mut. allele)
               ^ (size_t) mut. frameshift;
      }
  };  
  size_t stop () const
    { return pos + ref. size (); }

  // Requires: !prot
  bool isFrameshift () const
    { return difference (ref. size (), allele. size ()) % 3; }
  void replace (Dna &refDna) const;
    // Update: refDna
};

    

struct KmerIndex : Named, Singleton<KmerIndex>
// DNA k-mer index
// Assumptions for time: DNA sequence length is O(1)
//                       DNA identifier length is O(1)
//                       kmer_size = O(log(items))
{
  typedef  size_t  Addr;
  static_assert (! numeric_limits<Addr>::is_signed, "addr is signed");
  static_assert (sizeof (Addr) <= sizeof (streamsize), "Addr is larger than streamsize");
  static constexpr Addr nil {(Addr) -1};


  struct IdRecord
  {
    static constexpr size_t size {32};
  private:
    fstream& f;
  public:
    Addr addr {nil};
      // != nil
    Addr prev {nil};
    static constexpr size_t textSize {size - sizeof (prev)};
    static constexpr char no_char {'\0'};
    char text [textSize];
      // Reversed text
  private:
    size_t start {0};
      // <= textSize
  public:
      
    IdRecord (fstream &f_arg,
              Addr addr_arg,
              bool isNew);
      
    bool full () const;
    void save () const;
    void renew (Addr addr_arg);
  private:
    size_t getStart () const;
  public:
    void put (string &s);
      // Reversed s --> text[]
    string get () const;
      // Iteration over text[] from start to 0
  };
  
  
  static constexpr size_t progressSize {1000000};  // PAR
private:
  fstream f;
    // binary
public:
  const bool canRead;
  static constexpr size_t header_size {2};
  static constexpr size_t kmer_size_max {16};  // PAR
  const size_t kmer_size;
  const size_t code_max;
    // = 2 ^ (2 * kmer_size)
    // Use hashes of k-mers to make code_max O(1) ??!
  size_t items {0};
  Addr addr_new {0};
    // = file size of f
private:
  static constexpr char gap {' '};
public:
  
  
  KmerIndex (const string &name_arg,
             size_t kmer_size_arg);
    // Time: O(1) ??!
  explicit KmerIndex (const string &name_arg);
private:
  static size_t readKmerSize (fstream &fIn);
public:
  void qc () const final;
    // Requires: canRead
    
    
private:
  static size_t dna2code (const Dna &dna);
    // Time: O(dna.seq.size())
  static Addr code2addr (size_t code)
    { return header_size * sizeof (size_t) + code * sizeof (Addr); }
public:
  size_t getIdRecords () const
    { return (addr_new - code2addr (code_max)) / IdRecord::size; }
  size_t getKmers ();
    // Time: O(code_max)
  void add (const Dna &dna,
            size_t &kmers,
            size_t &kmersRejected);
    // Update: items++
    // Time: O(kmer_size)
    // For long sequences use only k-mers preceded by a specific prefix ??!
private:
  void addId (size_t code,
              const string &id);
public:
              
  struct NumId
  {
    size_t n {0};
      // > 0
    string id;  // DNA::getId()
      // !empty()
    NumId (size_t n_arg,
           const string &id_arg);
    bool operator< (const NumId &other) const;
  };
  Vector<NumId> find (const Dna &dna);
    // Return: sorted by NumId::n descending
    // Time: O(kmer_size) 
private:
  StringVector code2ids (size_t code);
    // Time: O(1)
};



}



#endif
