// amr_report.cpp

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
*   Identification of AMR genes using the AMR database and BLAST and HMM search results
*
*/
   
   
#undef NDEBUG 
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;
#include "gff.hpp"
using namespace GFF_sp;



static constexpr bool useCrossOrigin = false;  // GPipe: true



namespace 
{


constexpr static double frac_delta = 1e-5;  // PAR      




struct Batch;  // forward



struct Fam 
// Table PROTEUS.VF..FAM
{
	const Fam* parent {nullptr};
	  // Tree
  string id; 
  string genesymbol;
  string familyName; 
  bool reportable {false}; 

  // HMM
  string hmm; 
    // May be empty()
  double tc1 {NAN}; 
  double tc2 {NAN}; 


  Fam (Fam* parent_arg,
       const string &id_arg,
       const string &genesymbol_arg,
       const string &hmm_arg,
       double tc1_arg,
       double tc2_arg,
       const string &familyName_arg,
       bool reportable_arg)
    : parent (parent_arg)
    , id (id_arg)
    , genesymbol (genesymbol_arg)
    , familyName (familyName_arg)
    , reportable (reportable_arg)
    , hmm (hmm_arg)
    , tc1 (tc1_arg)
    , tc2 (tc2_arg)
    { if (genesymbol == "-")
        genesymbol. clear ();
      if (hmm == "-")
        hmm. clear ();
      if (hmm. empty () != ! tc1)
      { cout << id << ' ' << hmm << ' ' << tc1 << endl;
        ERROR;
      }
      ASSERT (hmm. empty () == ! tc2); 
    //IMPLY (! hmm. empty (), tc2 > 0);
      if (familyName == "NULL")
        familyName. clear ();
      ASSERT (tc2 >= 0);        
      ASSERT (tc2 <= tc1);
    }
  Fam ()
    {}
  void saveText (ostream &os) const 
    { os << hmm << " " << tc1 << " " << tc2 << " " << familyName << " " << reportable; }


	bool descendantOf (const Fam* ancestor) const
	  { if (! ancestor)
	  	  return true;
	  	if (this == ancestor)
	  		return true;
	  	if (parent)
	  		return parent->descendantOf (ancestor);
	  	return false;
	  }
  const Fam* getHmmFam () const
    // Return: most specific HMM
    { const Fam* f = this;
      while (f && f->hmm. empty ())
        f = static_cast <const Fam*> (f->parent);
      return f;
    }
};


map<string/*famId*/,const Fam*> famId2fam;



struct HmmAlignment 
// Query: AMR HMM
{
  string sseqid; 
  double score1 {NAN}; 
  double score2 {NAN}; 
    // May be different from max(Domain::score)
  const Fam* fam {nullptr};
    // Query
//ali_from, ali_to ??
  
  
  HmmAlignment (const string &line,
                const Batch &batch);
    // Update: batch.domains
  HmmAlignment ()
    {}
  void saveText (ostream &os) const 
    { os << sseqid << ' ' << score1 << ' ' << score2 << ' ' << (fam ? fam->hmm : string ()); }
      
      
  bool good () const
    { return    fam 
             && ! fam->hmm. empty ()
             && score1 >= fam->tc1
             && score2 >= fam->tc2
           //&& fam->reportable
             ; 
    }
private:
  bool betterEq (const HmmAlignment &other,
                 unsigned char criterion) const
    // Reflexive
    // For one sseqid: one HmmAlignment is better than all others
    { ASSERT (good ());
      ASSERT (other. fam);
      if (sseqid != other. sseqid)
        return false;
      switch (criterion)
      {
        case 0: return fam->descendantOf (other. fam); 
        case 1: 
          {
            LESS_PART (other, *this, score1);
          //return score2 >= other. score2;  // GP-16770
            LESS_PART (other, *this, fam->tc1);
            LESS_PART (*this, other, fam->id);  // Tie resolution
            return true;
          }
        default: ERROR;
      }
      throw logic_error ("Never call");
    }
public:
  bool better (const HmmAlignment &other,
               unsigned char criterion) const
    { return betterEq (other, criterion) && ! other. betterEq (*this, criterion); }


  typedef  pair<string/*sseqid*/,string/*FAM.id*/>  Pair;
  
  
  struct Domain  
  {
    double score {0};
    size_t hmmLen {0};
    size_t hmmStart {0}; 
    size_t hmmStop {0};
    size_t seqLen {0}; 
    size_t seqStart {0};
    size_t seqStop {0}; 
    

    Domain (const string &line,
            Batch &batch);
      // Input: line: hmmsearch -domtable line
    Domain ()
      {} 
  };
};



double ident_min = NAN;
double complete_cover_min = NAN;
double partial_cover_min = NAN;
bool cdsExist = false;
bool print_fam = false;

bool reportPseudo = false; 
const string stopCodonS ("[stop]");
const string frameShiftS ("[frameshift]");



struct BlastAlignment 
{
  // BLAST alignment
  size_t length {0}, nident {0}  // aa
       ,    refStart {0},    refEnd {0},    refLen {0}
       , targetStart {0}, targetEnd {0}, targetLen {0};
    // Positions are 0-based
    // targetStart < targetEnd

  // target    
  string targetName; 
  bool targetProt {true};
    // false <=> DNA
  bool targetStrand {true}; 
    // false <=> negative
  size_t targetAlign {0};
  size_t targetAlign_aa {0};
  bool partialDna {false};
  bool stopCodon {false}; 
  bool frameShift {false};
  
  // Reference (AMR) protein
  bool refStrand {true};
  long gi {0};  
    // 0 <=> HMM method
  string accessionProt; 
  string accessionDna;  
  size_t part {1};    
    // >= 1
    // <= parts
  size_t parts {1};  
    // >= 1
  // Table FAM
  string famId;  
  string gene;   
    // FAM.class  
  string mechanism; 
    // FAM.resistance
  const bool mutant {false}; 
  string product;  

  vector<Cds> cdss;
  
  static constexpr size_t mismatchTail_aa = 10;  // PAR
  size_t mismatchTailTarget {0};


  BlastAlignment (const string &line,
                  bool targetProt_arg)
    : targetProt (targetProt_arg)
    {
	    istringstream iss (line);
      string refName, targetSeq;
	    iss >> targetName >> refName >> length >> nident >> targetStart >> targetEnd >> targetLen >> refStart >> refEnd >> refLen >> targetSeq;
	  // format:  qseqid      sseqid    length    nident         qstart         qend         qlen      sstart      send      slen         sseq
    // blastp:  ...         ...          663       169              2          600          639           9       665       693          ...
    // blastx:  ...         ...          381       381          13407        14549        57298           1       381       381          ...
    // blastn:  ...         ...          733       733          62285        63017        88215         105       837       837          ...
	    ASSERT (! targetSeq. empty ());	    

      // refName	    
	    product       =                     rfindSplit (refName, '|'); 
	  //mutant        = //       str2<int> (rfindSplit (refName, '|')); // ??
	    mechanism     =                     rfindSplit (refName, '|');
	    gene          =                     rfindSplit (refName, '|');  // Reportable_vw.class
	    famId         =                     rfindSplit (refName, '|');  // Reportable_vw.fam
	    parts         = (size_t) str2<int> (rfindSplit (refName, '|'));
	    part          = (size_t) str2<int> (rfindSplit (refName, '|'));
	    accessionDna  =                     rfindSplit (refName, '|');
	    accessionProt =                     rfindSplit (refName, '|');
	    gi            = str2<long> (refName);
	    ASSERT (gi > 0);
	    	    
	    replace (product, '_', ' ');
	    
	    ASSERT (refStart != refEnd);
	    refStrand = refStart < refEnd;  
	    ASSERT (refStrand);  // ??
	    if (! refStrand)
	      swap (refStart, refEnd);

	    ASSERT (targetStart != targetEnd);
	    targetStrand = targetStart < targetEnd;  
	    IMPLY (targetProt, targetStrand);
	    if (! targetStrand)
	      swap (targetStart, targetEnd);
	      
	    ASSERT (refStart >= 1);
	    ASSERT (targetStart >= 1);
	    ASSERT (refStart < refEnd);
	    ASSERT (targetStart < targetEnd);
	    refStart--;
	    targetStart--;

	    partialDna = false;
	    constexpr size_t mismatchTailDna = 10;  // PAR
	    if (! targetProt && targetEnd - targetStart >= 30)  // PAR, PD-671
	    {
	           if (refStart > 0      && targetTail (true)  <= mismatchTailDna)  partialDna = true;
	      else if (refEnd   < refLen && targetTail (false) <= mismatchTailDna)  partialDna = true;
	    }

      setTargetAlign ();	    
	    
	    if (contains (targetSeq, "*"))  
	      stopCodon  =  true;
	  //frameShift = contains (targetSeq, "/");  // Needs "blastall -p blastx ... "
	    if (! targetProt && (targetEnd - targetStart) % 3 != 0)
  	    frameShift = true;	  
	  	  
	    // For BLASTX
	  	// PD-1280
	    if (   ! targetProt 
	        && refStart == 0 
	        && charInSet (targetSeq [0], "LIV") 
	        && nident < targetAlign_aa
	       )
	      nident++;
	    	    
	    mismatchTailTarget = mismatchTail_aa;
	    if (! targetProt)
	      mismatchTailTarget *= 3;
    }
  explicit BlastAlignment (const HmmAlignment& other)
    : targetName (other. sseqid)     
    , famId      (other. fam->id)   
    , gene       (other. fam->id)   
    , product    (other. fam->familyName)   
    { if (allele ())
        ERROR_MSG (famId + " " + gene);
    }
  void qc () const
    {
      if (! qc_on)
        return;
	    ASSERT (! famId. empty ());
	    ASSERT (! gene. empty ());
	    ASSERT (part >= 1);
	    ASSERT (part <= parts);
	    ASSERT (! product. empty ());
	    ASSERT (refStrand);
	    IMPLY (targetProt, targetStrand);  
	    ASSERT (targetStart < targetEnd);
	    ASSERT (targetEnd <= targetLen);
	    ASSERT ((bool) gi == (bool) length);
	    ASSERT ((bool) gi == (bool) refLen);
	    ASSERT ((bool) gi == (bool) nident);
	    ASSERT ((bool) gi == ! accessionProt. empty ());
	    IMPLY (! gi, getFam () -> getHmmFam ());
	    IMPLY (accessionProt. empty (), accessionDna. empty ());
	    IMPLY (targetProt, ! partialDna);
	  //IMPLY (! targetProt, (targetEnd - targetStart) % 3 == 0);
	    ASSERT (targetAlign);
	    IMPLY (targetProt, targetAlign == targetAlign_aa);
	    IMPLY (! targetProt, targetAlign == 3 * targetAlign_aa);
	    ASSERT (nident <= targetAlign_aa);
	    IMPLY (! targetProt, cdss. empty ());
	    if (gi)
	    {
	      ASSERT (refStart < refEnd);
  	    ASSERT (nident <= refEnd - refStart);
  	    ASSERT (refEnd <= refLen);
  	    ASSERT (refEnd - refStart <= length);	    
  	    ASSERT (targetAlign_aa <= length);
	    }
    }
  void saveText (ostream& os) const 
    { // PD-736, PD-774, PD-780, PD-799
      const string method (getMethod ());
      const string na ("NA");
      const string proteinName (refExactlyMatched () || ! gi ? product : nvl (getFam () -> familyName, na));
      ASSERT (! contains (proteinName, '\t'));
      if (! verbose ())  // Otherwise data may be not completely prepared
	      if (targetProt && ! (! cdsExist == cdss. empty ()))
	      {
	      	cout << endl << targetName << '!' << endl;
	      	ERROR;
	      }
      vector<Cds> cdss_ (cdss);
      if (cdss_. empty ())
        cdss_. push_back (Cds ());
      for (const Cds& cds : cdss_)
      {
        TabDel td (2, false);
      //if (targetProt)
          td << targetName;
        if (cdsExist)
          td << (cds. contig. empty () ? targetName  : cds. contig)
             << (cds. contig. empty () ? targetStart : cds. start) + 1
             << (cds. contig. empty () ? targetEnd   : cds. stop)
             << (cds. contig. empty () ? (/*refStrand*/ targetStrand ? '+' : '-') : (cds. strand ? '+' : '-'));
        td << (print_fam 
                 ? famId
                 : (method == "ALLELE" ? famId : nvl (getFam () -> genesymbol, na))
              )
           << proteinName + ifS (reportPseudo, ifS (stopCodon, " " + stopCodonS) + ifS (frameShift, " " + frameShiftS))
           << method
           << (targetProt ? targetLen : targetAlign_aa);  
        if (gi)
          td << refLen
             << refCoverage () * 100  
             << (/*targetProt ? pRefEffectiveLen () :*/ pIdentity ()) * 100  // refIdentity
             << length
             << accessionProt
             << product
           //<< accessionDna
             ;
        else
          td << na 
             << na
             << na
             << na
             << na
             << na;
        // PD-775
  	    if (const Fam* f = getFam () -> getHmmFam ())
          td << f->hmm
             << f->familyName;
        else
        {
          td << na
             << na;
          ASSERT (method != "HMM");
        }
        if (cdsExist)
        {
        	IMPLY (cds. crossOriginSeqLen, useCrossOrigin);
        	if (useCrossOrigin)
        	{
	          if (cds. crossOriginSeqLen) 
	            td << cds. crossOriginSeqLen;
	          else 
	            td << na;        	
	        }
        }
        os << td. str () << endl;
      }
    }
    

  bool allele () const
    { return famId != gene && parts == 1; }
  size_t targetTail (bool upstream) const
    { return targetStrand == upstream ? targetStart : (targetLen - targetEnd); }
  size_t refEffectiveLen () const
    { return partialDna ? refEnd - refStart : refLen; }
#if 0
  double pRefEffectiveLen () const
    { ASSERT (nident);
      return (double) nident / (double) refEffectiveLen ();
    }
#endif
  double pIdentity () const
    { return (double) nident / (double) length; }
  double refCoverage () const
    { return (double) (refEnd - refStart) / (double) refLen; }
  bool refExactlyMatched () const
    { return    refLen   
             && nident == refLen 
             && refLen == length;
	  }
  bool partial () const
    // Requires: good()
    { return refCoverage () < complete_cover_min - frac_delta; }
	string getMethod () const
	  { return refExactlyMatched () 
	             ? allele () && (! targetProt || refLen == targetLen)
	               ? "ALLELE"
	               : "EXACT"  // PD-776
	             : gi
	                ? partial ()
	                  ? "PARTIAL"
	                  : "BLAST"
	                : "HMM"; 
	  }
	  // PD-736
  bool good () const
    { if (! reportPseudo)
      {
        if (stopCodon)
          return false; 
        if (frameShift)
          return false; 
      }
    #if 0
	    if (targetProt)
	    { if (pRefEffectiveLen () < ident_min - frac_delta)  
  	      return false;
  	  }
  	  else
  	#endif
  	  { // PD-1032
  	    if (   pIdentity ()   < ident_min - frac_delta
  	        || refCoverage () < partial_cover_min - frac_delta
  	       )
  	      return false;
  	    if (parts > 1 && refCoverage () < complete_cover_min - frac_delta)
  	    	return false;
  	  }
	    return true;
    }
private:
  bool insideEq (const BlastAlignment &other) const
    { return    targetStrand                     == other. targetStrand
             && targetStart + mismatchTailTarget >= other. targetStart 
             && targetEnd                        <= other. targetEnd + mismatchTailTarget;
    }
    // Requires: same targetName
  bool descendantOf (const BlastAlignment &other) const
    { return    ! other. allele ()
             && getFam () -> descendantOf (other. getFam ());
    }
  bool betterEq (const BlastAlignment &other) const
    // Reflexive
    { if (targetName != other. targetName)
        return false;
      if (targetProt != other. targetProt)  // ??
        return false;
      // PD-727
    /*if (getFam () != other. getFam ())
      { if (descendantOf (other))
          return true;
        if (other. descendantOf (*this))
          return false;
      }*/
      // PD-807
    //if (! targetProt && ! other. insideEq (*this))
      if (   ! other. insideEq (*this)
      	  && !        insideEq (other)
      	 )
        return false;
      LESS_PART (other, *this, refExactlyMatched ());  // PD-1261, PD-1678
      LESS_PART (other, *this, nident);
      LESS_PART (*this, other, refEffectiveLen ());
      return true;
    }
public:
  const Fam* getFam () const
    { const Fam* fam = famId2fam [famId];
      if (! fam)
        fam = famId2fam [gene];
      ASSERT (fam);
      return fam;
    }
  bool better (const BlastAlignment &other) const
    { return    betterEq (other) 
    	       && (   ! other. betterEq (*this) 
    	           || accessionProt < other. accessionProt  // Tie resolution; PD-1245
    	          );
    }
  bool better (const HmmAlignment& other) const
    { ASSERT (other. good ());
      if (targetName != other. sseqid)
        return false;
      return    refExactlyMatched () 
           //|| getFam () -> getHmmFam () == other. fam
             || getFam () -> descendantOf (other. fam)
             ;
    }
  bool operator< (const BlastAlignment &other) const
    { LESS_PART (*this, other, targetName);
      LESS_PART (*this, other, targetStart);
      LESS_PART (*this, other, famId);
      LESS_PART (*this, other, accessionProt);
      return false;
    }
//size_t lengthScore () const
  //{ return refLen - (refEnd - refStart); } 
  void setTargetAlign ()
    { targetAlign = targetEnd - targetStart;
      targetAlign_aa = targetAlign;
      if (! targetProt)
      {
        ASSERT (targetAlign % 3 == 0);
        targetAlign_aa = targetAlign / 3;
      }
    }
};




struct Batch
{
  // Reference input
  map<string/*hmm*/,const Fam*> hmm2fam;
  Fam* root {nullptr};
    // Not delete'd in ~Batch()

  // Target input
  typedef  List<BlastAlignment>  BlastAls; 
  typedef  List<HmmAlignment>  HmmAls;   

  BlastAls blastAls;   
  map<HmmAlignment::Pair, HmmAlignment::Domain> domains;  // Best domain  
  HmmAls hmmAls;  
  
  // Output
  BlastAls goodBlastAls; 
  
  
  explicit Batch (const string &famFName)
	  : root (new Fam ())
	  {
	    if (famFName. empty ())
	    	throw runtime_error ("fam (protein family hierarchy) file is missing");

	  	// Tree of Fam
	  	// Pass 1  
	    {
	      LineInput f (famFName);  
	  	  while (f. nextLine ())
	  	  {
	  	    trim (f. line);
	      //cout << f. line << endl; 
	  	    const string famId               (findSplit (f. line, '\t'));
	  	    /*const string parentFamName =*/  findSplit (f. line, '\t');
	  	    const string genesymbol          (findSplit (f. line, '\t'));
	  	    const string hmm                 (findSplit (f. line, '\t'));
	  	    const double tc1 = str2<double>  (findSplit (f. line, '\t'));
	  	    const double tc2 = str2<double>  (findSplit (f. line, '\t'));
	  	    const int reportable = str2<int> (findSplit (f. line, '\t'));
	  	    ASSERT (   reportable == 0 
	  	            || reportable == 1
	  	           );
	  	    const auto fam = new Fam (root, famId, genesymbol, hmm, tc1, tc2, f. line, reportable);
	  	    famId2fam [famId] = fam;
	  	    if (! fam->hmm. empty ())
	  	      hmm2fam [fam->hmm] = fam;
	  	  }
	  	}
	  	// Pass 2
	    {
	      LineInput f (famFName);  
	  	  while (f. nextLine ())
	  	  {
	  	    trim (f. line);
	  	  //cout << f. line << endl;  
	  	    Fam* child = const_cast <Fam*> (famId2fam [findSplit (f. line, '\t')]);
	  	    ASSERT (child);
	  	    const string parentFamId (findSplit (f. line, '\t'));
	  	    Fam* parent = nullptr;
	  	    if (! parentFamId. empty ())
	  	      { EXEC_ASSERT (parent = const_cast <Fam*> (famId2fam [parentFamId])); }
	  	    child->parent = parent;
	  	  }
	  	}
	  }
	  	  

	void process (bool skip_hmm_check) 
  // Input: root, blastAls, domains, hmmAls
	// Output: goodBlastAls
	{
		ASSERT (root);
		
    if (! (ident_min >= 0 && ident_min <= 1))
    	throw runtime_error ("ident_min must be between 0 and 1");
    if (! (complete_cover_min >= 0 && complete_cover_min <= 1))
    	throw runtime_error ("complete_cover_min must be between 0 and 1");
    if (! (partial_cover_min >= 0 && partial_cover_min <= 1))
    	throw runtime_error ("partial_cover_min must be between 0 and 1");
    if (partial_cover_min > complete_cover_min)
    	throw runtime_error ("partial_cover_min msut be less than or equal to complete_cover_min");

    // Filtering by ::good() has been done above
        
    // Group by targetName and process each targetName separately for speed ??    

    // Pareto-better()  
	  for (const auto& blastAl : blastAls)
    {
    	ASSERT (blastAl. good ());
  	  bool found = false;
  	  for (const auto& goodBlastAl : goodBlastAls)
  	    if (goodBlastAl. better (blastAl))
	      {
	        found = true;
	        break;
	      }
	    if (found)
	      continue;	      
      for (Iter<BlastAls> goodIter (goodBlastAls); goodIter. next ();)
        if (blastAl. better (*goodIter))
          goodIter. erase ();          
      goodBlastAls << blastAl;	    
    }
  	if (verbose ())
  	{
  	  cout << "# Best Blasts: " << goodBlastAls. size () << endl;
  	  if (! goodBlastAls. empty ())
  	    goodBlastAls. front (). saveText (cout);
  	}
  #if 0
	  for (const auto& blastAl1 : goodBlastAls)
		  for (const auto& blastAl2 : goodBlastAls)
		  {
		  	cout << endl;
		  	blastAl1. saveText (cout);
		  	blastAl2. saveText (cout);
		  	cout        << blastAl1. better (blastAl2) 
		  	     << ' ' << blastAl2. better (blastAl1) 
		  	     << endl;
		  }
	#endif

    // Pareto-better()  
    HmmAls goodHmmAls; 
    FOR (unsigned char, criterion, 2)
    {
      // hmmAls --> goodHmmAls
      goodHmmAls. clear ();
  	  for (const auto& hmmAl : hmmAls)
      {
      	ASSERT (hmmAl. good ());
    	  bool found = false;
    	  for (const auto& goodHmmAl : goodHmmAls)
    	    if (goodHmmAl. better (hmmAl, criterion))
  	      {
  	        found = true;
  	        break;
  	      }
  	    if (found)
  	      continue;  
        for (Iter<HmmAls> goodIter (goodHmmAls); goodIter. next ();)
          if (hmmAl. better (*goodIter, criterion))
            goodIter. erase ();
        goodHmmAls << hmmAl;
      }
      //
      hmmAls = goodHmmAls;
      if (verbose ())
      {
        cout << "Pareto-better HMMs: (Criterion " << (int) criterion << "): " << hmmAls. size () << endl;
        for (const HmmAlignment& al : hmmAls)
        {
          al. saveText (cout);
          cout << endl;
        }
      }
    }

    // PD-741
  	if (! skip_hmm_check && ! goodHmmAls. empty ())
      for (Iter<BlastAls> iter (goodBlastAls); iter. next ();)
        if (   /*! iter->refExactlyMatched () */
        	     iter->pIdentity () < 0.98 - frac_delta  // PAR  // PD-1673
            && ! iter->partial ()
           )
	        if (const Fam* fam = iter->getFam () -> getHmmFam ())  
	        {
	          bool found = false;
	      	  for (const auto& hmmAl : goodHmmAls)
	            if (   iter->targetName == hmmAl. sseqid
	                && fam == hmmAl. fam
	               )
	            {
	              found = true;
	              break;
	            }
	          if (! found)   // BLAST is wrong
	            iter. erase ();
	        }
  	if (verbose ())
  	  cout << "# Best Blasts left: " << goodBlastAls. size () << endl;

    for (Iter<HmmAls> hmmIt (goodHmmAls); hmmIt. next ();)
  	  for (const auto& blastAl : goodBlastAls)
  	    if (blastAl. better (*hmmIt))
	      {
          hmmIt. erase ();
	        break;
	      }

    // Output 
    
    // goodHmmAls --> goodBlastAls
  	for (const auto& hmmAl : goodHmmAls)
  	{
  	  BlastAlignment al (hmmAl);
  	  if (verbose ())
  	    cout << al. targetName << " " << al. gene << endl;  
  	  const HmmAlignment::Domain domain = domains [HmmAlignment::Pair (al. targetName, al. gene)];
  	  if (! domain. hmmLen)  
  	    continue;  // domain does not exist
  	/*al. refLen      = domain. hmmLen;
  	  al. refStart    = domain. hmmStart;
  	  al. refEnd      = domain. hmmStop; */
  	  al. targetLen   = domain. seqLen;
  	  al. targetStart = domain. seqStart;
  	  al. targetEnd   = domain. seqStop;
  	  al. setTargetAlign ();
  	  ASSERT (! al. refExactlyMatched ());
  	  al. qc ();
  	  goodBlastAls << al;
  	}
  
    goodBlastAls. sort ();
	}
	
	
	
	void report (ostream &os) const
	// Input: goodBlastAls
	{
    // PD-283, PD-780
    // Cf. BlastAlignment::saveText()

    {
    	// Cf. BlastAlignment::saveText()
	    TabDel td;
	    td << "Target identifier";  // targetName
	    if (cdsExist)  
	      // Contig
	      td << "Contig id"
	         << "Start"  // targetStart
	         << "Stop"  // targetEnd
	         << "Strand";  // targetStrand
	    td << (print_fam ? "FAM.id" : "Gene symbol")
	       << "Protein name"
	       << "Method"
	       << "Target length" 
	       //
	       << "Reference protein length"         // refLen
	       << "% Coverage of reference protein"  // queryCoverage
	       << "% Identity to reference protein"  
	       << "Alignment length"                 // length
	       << "Accession of closest protein"     // accessionProt
	       << "Name of closest protein"
	       //
	       << "HMM id"
	       << "HMM description"
	       ;
      if (cdsExist)
	    	if (useCrossOrigin)
	      	 td << "Cross-origin length";
	    os << td. str () << endl;
	  }

  	for (const auto& blastAl : goodBlastAls)
  	  if (blastAl. getFam () -> reportable)
    	{
    	  blastAl. qc ();
    	  blastAl. saveText (os);
    	}
	}



	void printTargetIds (ostream &os) const
	{
		ASSERT (os. good ());
  	for (const auto& blastAl : goodBlastAls)
  	  if (blastAl. getFam () -> reportable)
        os << blastAl. targetName << endl;
	}
};




// HmmAlignment

HmmAlignment::HmmAlignment (const string &line,
                            const Batch &batch)
{
  istringstream iss (line);
  string hmm, dummy;
  //                                                        --- full sequence ---  --- best 1 domain --
  //     target name  accession  query name     accession   E-value  score     bias     E-value  score  bias   exp reg clu  ov env dom rep inc description of target
  iss >> sseqid >>    dummy      >> dummy      >> hmm >>    dummy >> score1 >> dummy >> dummy >> score2;
  ASSERT (score1 > 0);
  ASSERT (score2 > 0)
  find (batch. hmm2fam, hmm, fam);
}




// Domain

HmmAlignment::Domain::Domain (const string &line,
                              Batch &batch)
{
  istringstream iss (line);
  string target_name, accession, query_name, query_accession;
  size_t n, of, env_from, env_to;
  double eValue, full_score, full_bias, cValue, i_eValue, domain_bias, accuracy;
  iss >> target_name >> accession >> seqLen 
      >> query_name >> query_accession >> hmmLen 
      >> eValue >> full_score >> full_bias 
      >> n >> of >> cValue >> i_eValue >> score >> domain_bias 
      >> hmmStart >> hmmStop 
      >> seqStart >> seqStop 
      >> env_from >> env_to
      >> accuracy;
  ASSERT (accession == "-");
  ASSERT (hmmStart);
  ASSERT (seqStart);
  ASSERT (env_from);
  hmmStart--;
  seqStart--;
  env_from--;
  ASSERT (hmmStart < hmmStop);
  ASSERT (seqStart < seqStop);
  ASSERT (hmmStop <= hmmLen);
  ASSERT (seqStop <= seqLen);
  ASSERT (full_score > 0);
  ASSERT (n >= 1);
  ASSERT (n <= of);
  ASSERT (score > 0);

  const Fam* fam = batch. hmm2fam [query_accession/*query_name*/];
  if (! fam)
    return;
  const HmmAlignment::Pair p (target_name, fam->id);
  const HmmAlignment::Domain domain_old (batch. domains [p]);
  if (domain_old. score > score)
    return;
  batch. domains [p] = *this;
}




// ThisApplication

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Report AMR proteins")
    {
      // Input
      addKey ("fam", "Table FAM");
      const string blastFormat ("qseqid sseqid length nident qstart qend qlen sstart send slen sseq. qseqid format: gi|Protein accession|DNA accession|fusion part|# fusions|FAM.id|FAM.class|resistance mechanism|Product name");
      addKey ("blastp", "blastp output in the format: " + blastFormat);  
      addKey ("blastx", "blastx output in the format: " + blastFormat);  
      addKey ("gff", ".gff assembly file");
      addKey ("gff_match", ".gff-FASTA matching file for \"locus_tag\": \"<FASTA id> <locus_tag>\"");
      addKey ("hmmdom", "HMM domain alignments");
      addKey ("hmmsearch", "Output of hmmsearch");
      addKey ("ident_min", "Min. identity to the reference protein (0..1)", "0.9");
      addKey ("complete_cover_min", "Min. coverage of the reference protein (0..1) for a complete hit", "0.9");
      addKey ("partial_cover_min", "Min. coverage of the reference protein (0..1) for a partial hit", "0.5");
      addFlag ("skip_hmm_check", "Skip checking HMM for a BLAST hit");
      // Output
      addKey ("out", "Identifiers of the reported input proteins");
      addFlag ("print_fam", "Print the FAM.id instead of gene symbol"); 
      addFlag ("pseudo", "Indicate pseudo-genes in the protein name as \"" + stopCodonS + "\" or \"" + frameShiftS + "\""); 
      // Testing
      addFlag ("nosame", "Exclude the same reference ptotein from the BLAST output (for testing)"); 
      addFlag ("noblast", "Exclude the BLAST output (for testing)"); 
    #ifdef SVN_REV
      addFlag ("version", "Print the version and exit");
    #endif
    }



  void body () const final
  {
    const string famFName      = getArg ("fam");
    const string blastpFName   = getArg ("blastp");
    const string blastxFName   = getArg ("blastx");
    const string gffFName      = getArg ("gff");
    const string gffMatchFName = getArg ("gff_match");
    const string hmmDom        = getArg ("hmmdom");
    const string hmmsearch     = getArg ("hmmsearch");  
                 ident_min          = str2<double> (getArg ("ident_min"));  
                 complete_cover_min = str2<double> (getArg ("complete_cover_min"));  
                 partial_cover_min  = str2<double> (getArg ("partial_cover_min")); 
    const bool skip_hmm_check  = getFlag ("skip_hmm_check"); 
    const string outFName      = getArg ("out");
                 print_fam     = getFlag ("print_fam");
                 reportPseudo  = getFlag ("pseudo");
    const bool nosame          = getFlag ("nosame");
    const bool noblast         = getFlag ("noblast");
  #ifdef SVN_REV
    const bool version         = getFlag ("version");
  #endif
    
    ASSERT (hmmsearch. empty () == hmmDom. empty ());
    
    
  #ifdef SVN_REV
    if (version) 
    {
      cout << "version: " << SVN_REV << endl;
      exit(1);
    }    
  #endif
    
    cdsExist =    ! blastxFName. empty ()
               || ! gffFName. empty ();


    Batch batch (famFName);  
  
  
    // Fusion proteins, see PD-283 ??
    
    
    // Input 

    // batch.blastAls
  	// BlastAlignment::good()
    if (! noblast)
    {
      if (! blastpFName. empty ())
      {
        LineInput f (blastpFName);
    	  while (f. nextLine ())
    	  {
    	    { 
    	      Unverbose unv;
    	      if (verbose ())
    	        cout << f. line << endl;  
    	    }
    	    const BlastAlignment al (f. line, true);
    	    al. qc ();  
    	    if (nosame && toString (al. gi) == al. targetName)
    	      continue;
    	    if (al. good ())
    	      batch. blastAls << al;
    	    else
      	    { ASSERT (! al. refExactlyMatched ()); }
    	  }
    	}

      if (! blastxFName. empty ())
      {
        LineInput f (blastxFName);
    	  while (f. nextLine ())
    	  {
    	    { 
    	      Unverbose unv;
    	      if (verbose ())
    	        cout << f. line << endl;  
    	    }
    	    const BlastAlignment al (f. line, false);
    	    al. qc ();  
    	    if (nosame && toString (al. gi) == al. targetName)
    	      continue;
    	    if (al. good ())
    	      batch. blastAls << al;
    	    else
      	    { ASSERT (! al. refExactlyMatched ()); }
    	  }
    	}
    }
  	if (verbose ())
  	  cout << "# Good Blasts: " << batch. blastAls. size () << endl;
  	
  
    // batch.domains
    if (! hmmDom. empty ())
    {
      LineInput f (hmmDom);
  	  while (f. nextLine ())
  	  {
  	    trim (f. line);
  	    if (   f. line. empty () 
  	        || f. line [0] == '#'
  	       )
  	      continue;
  	    const HmmAlignment::Domain domain (f. line, batch);
  	  }
    }


    // batch.hmmAls 
  	// HmmAlignment::good()
  	if (! hmmsearch. empty ())  // redundant file ??
  	{
      LineInput f (hmmsearch);
  	  while (f. nextLine ())
  	  {
  	    if (verbose ())
  	      cout << f. line << endl;  
  	    if (f. line. empty () || f. line [0] == '#')
  	      continue;
  	    const HmmAlignment hmmAl (f. line, batch);
  	    if (hmmAl. good ())
  	      batch. hmmAls << hmmAl;
  	  }
  	}
  	if (verbose ())
  	  cout << "# Good HMMs: " << batch. hmmAls. size () << endl;
  
  
    // Output
    batch. process (skip_hmm_check);    


    // For Batch::report()
    if (! gffFName. empty ())
    {
	    const Gff gff (gffFName, ! gffMatchFName. empty ());
	    map<string/*seqid*/,string/*locusTag*/> seqId2locusTag;
	    if (! gffMatchFName. empty ())
	    {
	      LineInput f (gffMatchFName);
	  	  while (f. nextLine ())
	  	  {
	  	  	istringstream iss (f. line);
	  	  	string seqId, locusTag;
	  	  	iss >> seqId >> locusTag;
	  	  	ASSERT (! locusTag. empty ());
	  	  	seqId2locusTag [seqId] = locusTag;
	  	  }
	    }
    	for (auto& al : batch. goodBlastAls)
    	  if (al. targetProt)
      	{
      	  ASSERT (al. cdss. empty ());
      	#if 0
      	  string s (al. targetName);
      	  trimSuffix (s, "|");
      	  string locusTag (rfindSplit (s, '|'));
      	#else
      	  string locusTag = al. targetName;
      	#endif
      	  if (! gffMatchFName. empty ())
      	  	locusTag = seqId2locusTag [locusTag];
      	  if (const Set<Cds>* cdss = findPtr (gff. seqid2cdss, locusTag))
	      	  insertAll (al. cdss, *cdss);
      	  else
      	    throw runtime_error ("Locus tag \"" + locusTag + "\" is misssing in .gff file. Protein name: " + al. targetName);
      	  al. qc ();
      	}
    }
    
    
    batch. report (cout);


    if (! outFName. empty ())
    {
	    OFStream ofs (outFName);
      batch. printTargetIds (ofs);
    }
  }
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



