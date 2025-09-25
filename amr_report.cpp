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
*   Identification of AMR genes using BLAST and HMM search results vs. the AMRFinderPlus database
*
*/
   
   

#undef NDEBUG 

#include "common.hpp"
#include "tsv.hpp"
using namespace Common_sp;
#include "seq.hpp"
using namespace Seq_sp;
#include "gff.hpp"
using namespace GFF_sp;
#include "alignment.hpp"
using namespace Alignment_sp;
#include "columns.hpp"

#include "common.inc"



namespace 
{


// Global

// PAR
constexpr bool useCrossOrigin = false;  // GPipe: true
constexpr double frac_delta = 1e-5;  
  
bool ident_min_user = false;
bool equidistant = false;
bool cdsExist = false;
bool print_node = false;
bool print_node_raw = false;
//bool reportPseudo = false; 
bool targetProt = true;
string input_name;

//const string stopCodonS ("[stop]");
//const string frameShiftS ("[frameshift]");

map <string/*accession*/, Vector<AmrMutation>>  accession2mutations;



struct BlastRule final : Root
// PD-2310
{
  // 0 <=> undefined
  // 0 .. 1
  double ident {0.0};
    // Of alignment
  double target_coverage {0.0};  // Not used
  double ref_coverage {0.0};


  BlastRule (double ident_arg,
             double ref_coverage_arg)
    : ident        (ident_arg)
    , ref_coverage (ref_coverage_arg)
    {
      QC_ASSERT (ident > 0.0);
      QC_ASSERT (ident <= 1.0);
      QC_ASSERT (target_coverage >= 0.0);
      QC_ASSERT (target_coverage <= 1.0);
      QC_ASSERT (ref_coverage > 0.0);
      QC_ASSERT (ref_coverage <= 1.0);
    }  
  BlastRule () = default;
  void saveText (ostream &os) const final
    { os << ident << ' ' << ref_coverage; }
  bool empty () const final
    { return ! ident; }
};

BlastRule defaultCompleteBR;
BlastRule defaultPartialBR;



struct Fam 
// Table PROTEUS.VF..FAM
{
	const Fam* parent {nullptr};
	  // Tree
  string id; 
  string genesymbol;
  string familyName; 
  uchar reportable {0}; 
    // 0,1,2

  // HMM
  string hmm; 
    // May be empty()
  double tc1 {NaN}; 
  double tc2 {NaN}; 

  // BlastRule's
  BlastRule completeBR;
  BlastRule partialBR;
  
  string type;
  string subtype;
  string classS;
  string subclass;


  Fam (const string &id_arg,
       const string &genesymbol_arg,
       const string &hmm_arg,
       double tc1_arg,
       double tc2_arg,
       const BlastRule &completeBR_arg,
       const BlastRule &partialBR_arg,
       const string &type_arg,
       const string &subtype_arg,
       const string &class_arg,
       const string &subclass_arg,
       const string &familyName_arg,
       uchar reportable_arg)
    : id (id_arg)
    , genesymbol (genesymbol_arg)
    , familyName (familyName_arg)
    , reportable (reportable_arg)
    , hmm (hmm_arg)
    , tc1 (tc1_arg)
    , tc2 (tc2_arg)
    , completeBR (completeBR_arg)
    , partialBR (partialBR_arg)
    , type (type_arg)
    , subtype (subtype_arg)
    , classS (class_arg)
    , subclass (subclass_arg)
    { if (genesymbol == "-")
        genesymbol. clear ();
      if (hmm == "-")
        hmm. clear ();
      QC_ASSERT (hmm. empty () == ! tc1);
      QC_ASSERT (hmm. empty () == ! tc2); 
    //IMPLY (! hmm. empty (), tc2 > 0);
      if (familyName == "NULL")
        familyName. clear ();
      QC_ASSERT (tc2 >= 0.0);        
      QC_ASSERT (tc2 <= tc1);
      QC_ASSERT (! completeBR. empty ())
      QC_ASSERT (! partialBR. empty ());
    }
  Fam () = default;
  void saveText (ostream &os) const 
    { os << hmm << " " << tc1 << " " << tc2 << " " << familyName << " " << (int) reportable; }


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
        f = f->parent;
      return f;
    }
};



map<string/*famId*/,const Fam*> famId2fam;
  // Value: !nullptr
  //        Not delete'd



struct Batch; 
struct BlastAlignment;



struct HmmAlignment 
// Query: AMR HMM
{
  string sseqid; 
  double score1 {NaN}; 
  double score2 {NaN}; 
    // May be different from max(Domain::score)
  const Fam* fam {nullptr};
    // Query
    // !nullptr
//ali_from, ali_to ??
  unique_ptr<const BlastAlignment> blastAl;
    // (bool)get()
  
  
  HmmAlignment (const string &line,
                const Batch &batch);
    // Update: batch.domains
  void saveText (ostream &os) const 
    { os << sseqid << ' ' << score1 << ' ' << score2 << ' ' << (fam ? fam->hmm : noString); }
      
      
  bool good () const
    { QC_ASSERT (fam);
      QC_ASSERT (! fam->hmm. empty ());
    	return    score1 >= fam->tc1
             && score2 >= fam->tc2; 
    }
private:
  bool betterEq (const HmmAlignment &other,
                 unsigned char criterion) const
    // Reflexive
    // Must: transitive
    // For one sseqid: one HmmAlignment is better than all others
    { ASSERT (good ());
      ASSERT (other. good ());
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
            if (! equidistant)
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
    { return             betterEq (other, criterion) 
             && ! other. betterEq (*this, criterion); }
  bool better (const BlastAlignment &other) const;


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
    Domain () = default;
  };
};



struct Susceptible final : Root
{
	// !empty()
	string genesymbol;
	double cutoff;
	  // 0 <=> broken gene confers resistance
	string classS;
	string subclass;
	string name;
	
	
	Susceptible (const string &genesymbol_arg,
	             double cutoff_arg,
  			 			 const string &class_arg,
  				 		 const string &subclass_arg,
  						 const string &name_arg)
    : genesymbol (genesymbol_arg)
    , cutoff (cutoff_arg / 100.0)
    , classS (class_arg)
    , subclass (subclass_arg)
    , name (name_arg)
    { 
      QC_ASSERT (cutoff >= 0.0);
      QC_ASSERT (cutoff < 1.0);
    	QC_ASSERT (! name. empty ());
      QC_ASSERT (! contains (name, '\t'));
      replace (name, '_', ' ');
      QC_ASSERT (! contains (name, "  "));
    }
  Susceptible () = default;
  Susceptible& operator= (Susceptible&& other) = default;
  void saveText (ostream &os) const final
    { os         << genesymbol 
         << '\t' << cutoff 
         << '\t' << classS
         << '\t' << subclass
         << '\t' << name
         << endl;
    }
};


map <string/*accession*/, Susceptible>  accession2susceptible;



struct BlastAlignment final : Alignment
// BLASTP or BLASTX
{
  // Target  
  bool partialDna {false};
  
  // Reference protein
  const bool fromHmm;
  string refAccession; 
    // empty() <=> HMM method
  size_t part {1};    
    // >= 1
    // <= parts
  size_t parts {1};  
    // >= 1
    // > 1 <=> fusion protein
  VectorPtr<BlastAlignment> fusions;
  bool fusionRedundant {false};
  // Table FAM
  string famId;  
  string gene;   
    // FAM.class  
  string resistance;
  uchar reportable {0};
  string classS;
  string subclass;

  const Fam* brFam {nullptr};
  // Valid if !fromHmm and inFam()
  BlastRule completeBR;  
  BlastRule partialBR;   
  
  string product;  
  Vector<Locus> cdss;
  static constexpr size_t mismatchTail_aa = 10;  // PAR
  
  const Susceptible* susceptible {nullptr};
    // In accession2susceptible
    // Only for the right organism
  
  const HmmAlignment* hmmAl {nullptr};


  BlastAlignment (const string &line,
                  bool sProt_arg)
    : Alignment (line, true, sProt_arg)
    , fromHmm (false)
    {
    	try 
    	{
        try
        {	
		      // qseqid	    
			    product                     =                     rfindSplit (qseqid, '|'); 
			    classS                      =                     rfindSplit (qseqid, '|'); 
			    subclass                    =                     rfindSplit (qseqid, '|'); 
			    reportable                  = (uchar) str2<int>  (rfindSplit (qseqid, '|')); 
			    resistance                  =                     rfindSplit (qseqid, '|'); 
			    gene                        =                     rfindSplit (qseqid, '|');  // Reportable_vw.class
			    famId                       =                     rfindSplit (qseqid, '|');  // Reportable_vw.fam
			    parts                       = (size_t) str2<int> (rfindSplit (qseqid, '|'));
			    part                        = (size_t) str2<int> (rfindSplit (qseqid, '|'));
			    refAccession                =                                 qseqid;  // rfindSplit (qseqid, '|');
			  //str2<long> (qseqid);  // gi: dummy
  		    if (contains (refAccession, ':'))  
  		    {
  		      QC_ASSERT (isMutationProt ());
  		      const string geneMutation =               rfindSplit (refAccession, ':');
  		      const size_t pos          = str2<size_t> (rfindSplit (refAccession, ':'));
  		      ASSERT (refMutation. empty ());
  		      refMutation = std::move (AmrMutation (pos, geneMutation));
  		      QC_ASSERT (! refMutation. empty ());
  		      refMutation. qc ();
  		    }
			  }
			  catch (const exception &e)
			  {
			  	throw runtime_error (string ("Bad AMRFinder database\n") + e. what () + "\n" + line);
			  }
		    QC_ASSERT (! refAccession. empty ());
		  //qseqid = refAccession;		    	    
		    replace (product,  '_', ' ');
		    replace (classS,   '_', ' ');
		    replace (subclass, '_', ' ');
        if (isSusceptibleProt ())
		      susceptible = findPtr (accession2susceptible, refAccession);
  	    if (isMutationProt ())
  		    if (const Vector<AmrMutation>* refMutations = findPtr (accession2mutations, refAccession))
  		    {
  		    	if (verbose ())
  		        cout << "AmrMutation protein found: " << refAccession << endl << line << endl;
    	      setSeqChanges (*refMutations, 0);
    	      if (verbose ())
    	        cout << endl;
  		    }
  		  if (sProt)
  		    finish ();
		  }
		  catch (const exception &e)
		  {
		  	throw runtime_error (line + "\n" + e. what ());;
		  }
    }
  BlastAlignment (const HmmAlignment& hmmAl_arg,
                  const BlastAlignment* best)
    : fromHmm    (true)
    , famId      (hmmAl_arg. fam->id)   
    , gene       (hmmAl_arg. fam->id)   
    , product    (hmmAl_arg. fam->familyName) 
    , hmmAl      (& hmmAl_arg)  
    { ASSERT (hmmAl_arg. good ());
      sProt = true;
      if (best)
      {
        Hsp::operator= (*best);
        refAccession = best->refAccession;
      }
      sseqid = hmmAl_arg. sseqid;
      finishHsp (false, false);
      ASSERT (sProt);
      ASSERT (! allele ());
    }
  void finish ()
    {
      ASSERT (! fromHmm);      
      // BlastRule
	    // PD-2310
      if (inFam ())
      {
		    // PD-4856
		    ASSERT (! brFam);
		    // brFam
        EXEC_ASSERT (brFam = getFam ()); 
        while (brFam)
        {
          if (partial ())
          {
            if (passBlastRule (brFam->partialBR))
          	  break;
          }
          else
            if (passBlastRule (brFam->completeBR))
          	  break;
          brFam = brFam->parent;
        }
		  }
		  else
		  {
	      completeBR = defaultCompleteBR;
	      partialBR  = defaultPartialBR;
	    }
	    //
	    partialDna = false;
	    constexpr size_t mismatchTailDna = 10;  // PAR
	    if (! sProt && sAbsCoverage () >= 30)  // PAR, PD-671
	    {
	           if (qInt. start > 0   && sTail (true)  <= mismatchTailDna)  partialDna = true;
	      else if (qInt. stop < qlen && sTail (false) <= mismatchTailDna)  partialDna = true;
	    }
      //
	    if (! sProt)
	      cdss. emplace_back (0, sseqid, sInt. start, sInt. stop, sInt. strand == 1, partialDna, 0, noString, noString);
    }
  void qc () const override
    {
      if (! qc_on)
        return;
      Alignment::qc ();
      if (empty ())
        return;
      QC_ASSERT (qProt);
	    QC_ASSERT (! famId. empty ());
	    QC_ASSERT (! gene. empty ());
	    QC_ASSERT (part >= 1);
	    QC_ASSERT (part <= parts);
	  //QC_IMPLY (part > 1, ! fusions. empty () || fusionRedundant);
	    QC_IMPLY (! fusions. empty (), parts >= 2);
	    QC_IMPLY (isMutationProt (), ! isSusceptibleProt ());
	    QC_IMPLY (susceptible, isSusceptibleProt ());
	    QC_IMPLY (isStrongSusceptibleProt (), isSusceptibleProt ());
	    QC_IMPLY (isSusceptibleProt (), fusions. empty () && ! fusionRedundant);
	    QC_IMPLY (isMutationProt (), fusions. empty () && ! fusionRedundant);
	    QC_IMPLY (! isMutationProt () && ! seqChanges. empty () && ! seqChanges. front (). empty (), isStrongSusceptibleProt ());
	    QC_IMPLY (hasDeclarativeFrameshift (), isMutationProt ());
	    QC_IMPLY (! fusions. empty (), fusions. size () >= 2);
	    QC_IMPLY (! fusions. empty (), fusions. front () == this);
	    QC_IMPLY (fusionRedundant, fusions. empty ());
	    QC_ASSERT (! product. empty ());
	    QC_IMPLY (! fromHmm, ! refAccession. empty ());
	    QC_ASSERT (refAccession. empty () == sseq. empty ());
	    QC_ASSERT (! refAccession. empty () == (bool) qlen);
	    QC_ASSERT (! refAccession. empty () == (bool) nident);
	    QC_IMPLY (refAccession. empty () && inFam (), ! getFam () -> hmm. empty ());
	    QC_IMPLY (sProt, ! partialDna);
	    QC_IMPLY (! refMutation. empty (), isMutationProt ());
	    if (! sProt)
    	  for (const Locus& cds : cdss)
    	    QC_ASSERT (cds. contig == sseqid);
	    QC_IMPLY (! seqChanges. empty (), isMutationProt () || isStrongSusceptibleProt ());
    }
  void report (TsvOut &td,
               const string &targetName,
               bool mutationAll) const 
    { // PD-736, PD-774, PD-780, PD-799
      const string proteinName (susceptible
                                  ? susceptible->name
                                  : isMutationProt () 
                                    ? product 
                                    : refProtExactlyMatched (true) || parts >= 2 || refAccession. empty ()   // PD-3187, PD-3192
                                      ? product 
                                      : nvl (checkPtr (getMatchFam ()) -> familyName, na)
                               );      
      ASSERT (! contains (proteinName, '\t'));
      Vector<Locus> cdss_ (cdss);
      if (cdss_. empty ())
        cdss_ << Locus ();
      Vector<SeqChange> seqChanges_ (seqChanges);
      if (seqChanges_. empty ())
        seqChanges_ << SeqChange ();
      for (const Locus& cds : cdss_)
      {
        if (sProt)
        {
          if (targetProt)
            { ASSERT (sseqid == targetName); }
          else
            if (cds. contig != targetName)
              continue;
        }
        else
          { ASSERT (sseqid == targetName); }
	      for (const SeqChange& seqChange : seqChanges_)
	      {
          const string method (empty () ? na : getMethod (cds));
          VectorPtr<AmrMutation> mutations (seqChange. mutations);
          if (mutations. empty ())
            mutations << nullptr;
          for (const AmrMutation* mut : mutations)
          {
            QC_IMPLY (! verbose () && isMutationProt (), ! seqChange. empty () || mut);
  	        const bool isMutation = ! seqChange. empty () && mut && ! seqChange. replacement;
  	        if (mutationAll)
  	        {
  	          if (! isMutationProt ())
  	            continue;
  	        }
  	        else if (isMutationProt () && ! isMutation)
  	          continue;
      	    if (! input_name. empty ())
      	      td << input_name;
            td << (sProt ? nvl (sseqid, na) : na);  // PD-2534
  	        if (cdsExist)
  	          td << (empty () ? na : (cds. contig. empty () ? nvl (sseqid, na)  : cds. contig))
  	             << (empty () ? 0  : (cds. contig. empty () ? sInt. start : cds. start) + 1)
  	             << (empty () ? 0  : (cds. contig. empty () ? sInt. stop  : cds. stop))
  	             << (empty () ? na : string (1, cds. contig. empty () ? strand2char (sInt. strand) : (cds. strand ? '+' : '-')));
  	        td << (isMutationProt ()
  			             ? mut 
  			               ? seqChange. empty ()
                         ? mut->wildtype ()
  			                 : mut->geneMutation
  			               : gene + "_" + seqChange. getMutationStr ()
 			               : seqChange. empty ()
                       ? fusion2geneSymbols ()
                       : gene + disruption_delim + seqChange. getMutationStr ()  // isStrongSusceptibleProt()
                  )
  	           << (isMutationProt ()
  	                 ? mut 
                       ? seqChange. empty ()
                         ? proteinName + " [WILDTYPE]"
                         : mut->name
                       : proteinName + " [UNKNOWN]"
  	                 : proteinName 
  	              )
  	           << (isMutationProt () || fusion2core () ? "core" : "plus");  // PD-2825
            // PD-1856
  	        if (isMutationProt ())
  	          td << "AMR"
  	             << "POINT"
  	             << (mut ? nvl (mut->classS,   na) : na)
  	             << (mut ? nvl (mut->subclass, na) : na);
  	        else
  	          td << nvl (fusion2type (), na)  
    	           << (isStrongSusceptibleProt () ? "POINT_DISRUPT" : nvl (fusion2subtype (), na))
    	           << nvl (fusion2class (), na)
    	           << nvl (fusion2subclass (), na);
  	        td << method
  	           << (sProt ? slen : (sAbsCoverage () / 3));  // Approximate for disrupted genes
  	        if (refAccession. empty ())
  	          td << na 
  	             << na
  	             << na
  	             << na
  	             << na
  	             << na;
  	        else
  	          td << qlen
  	             << qRelCoverage () * 100.0  
  	             << relIdentity () * 100.0  // refIdentity
  	             << sseq. size ()
  	             << refAccession
  	             << product 
  	             ;
  	        // PD-775
  	        if (isMutationProt ())
  	          td << na
  	             << na;
  	        else 
  	        {
  	          if (const HmmAlignment* hmmAl_ = fusion2hmmAl ())
    	          td << hmmAl_->fam->hmm
    	             << hmmAl_->fam->familyName;
    	        else
    	        {
    	          td << na
    	             << na;
    	          ASSERT (method != "HMM");
    	        }
    	      }
  	        if (cdsExist)
  	        {
  	        	IMPLY (cds. crossOrigin, useCrossOrigin);
  	        	if (useCrossOrigin)
  	        	{
  		          if (cds. crossOrigin) 
  		            td << cds. contigLen;
  		          else 
  		            td << na;        	
  		        }
  	        }
  	        if (print_node)
  	        {
  	          if (! print_node_raw && allele () && ! refProtExactlyMatched (true))
  	            td << gene;
  	          else
  	            td << fusion2famIds ();
  	        }
  	        IMPLY (isMutationProt () && isMutation /*! seqChange. empty () && mut && ! seqChange. replacement*/, hasMutation ());
  	        td. newLn ();
  	      }
	      }
	    }
    }
    

  // Sequence mutations confer resistsance
  bool isMutationProt () const
    { return resistance == "mutation"; }
  bool isSusceptibleProt () const
    { return resistance == "susceptible"; }  
  bool isStrongSusceptibleProt () const
    { return    susceptible 
             && ! susceptible->cutoff; 
    }  
  //
  bool inFam () const
    { return    ! isMutationProt () 
             && ! isSusceptibleProt (); 
    } 
#if 0
  bool notInFam () const
    { return ! inFam (); }
#endif
  Set<string> getMutationSymbols () const
    { Set<string> mutationSymbols;
      for (const SeqChange& seqChange : seqChanges)
        if (! seqChange. empty ())
          for (const AmrMutation* mutation : seqChange. mutations)
            mutationSymbols << mutation->geneMutation;
      return mutationSymbols;
    }        
  bool allele () const
    { return famId != gene && parts == 1; }
  size_t refEffectiveLen () const
    { return partialDna ? qInt. len () : qlen; }
  bool partial () const
    // Requires: good()
    { return qRelCoverage () < defaultCompleteBR. ref_coverage - frac_delta; }  
  bool getTargetStrand (const Locus &cds) const
    { return sProt
               ? cds. empty ()
                 ? true
                 : cds. strand
               : (sInt. strand == 1);
    }
  size_t missedDnaStart (const Locus &cds) const
    { return 3 * (getTargetStrand (cds)
                    ? qInt. start
                    : qlen - qInt. stop
                 );
    }
  size_t missedDnaStop (const Locus &cds) const
    { return 3 * (getTargetStrand (cds)
                    ? qlen - qInt. stop
                    : qInt. start
                 );
    }
  bool truncated (const Locus &cds) const
    { return    (missedDnaStart (cds) > 0 && (sProt ? (cds. empty () ? false : cds. atContigStart ()) : sInt. start       <= Locus::end_delta))
             || (missedDnaStop  (cds) > 0 && (sProt ? (cds. empty () ? false : cds. atContigStop  ()) : slen - sInt. stop <= Locus::end_delta));
    }
  bool truncatedCds () const
    { for (const Locus& cds : cdss)
        if (truncated (cds))
          return true;
      return false;
    }
  bool refProtExactlyMatched (bool targetComplete) const
    { return    qProt
             && qlen   
             && nident == qlen 
             && nident == sseq. size ()
             && c_complete != efalse
             && (! sProt || ! targetComplete || qlen + (c_complete == etrue ? 1 : 0) == slen)
             && disrs. empty ();
	  }
  bool alleleMatch () const
    { return refProtExactlyMatched (true) && allele (); }
  bool alleleReportable () const  // PD-3583
    { return alleleMatch () && reportable >= 2; }
  uchar getReportable () const
    { if (const Fam* f = getMatchFam ())
        return f->reportable;
      return 0;
    }
  uchar fusion2reportable () const
    { ASSERT (! isMutationProt ());
      if (susceptible)
        return 2;
      if (fusions. empty ())
        return getReportable ();
      uchar reportable_max = 0;
      for (const BlastAlignment* fusion : fusions)
        maximize (reportable_max, fusion->getReportable ());
      return reportable_max;
    }
  bool fusionOverrides (const BlastAlignment &other) const
    { ASSERT (sProt == other. sProt);
      return    parts > 1
  	         && other. parts == 1
  	         && other. insideEq (*this)
  	         && (   sInt. start + Cds::peptideSize_min * a2s <= other. sInt. start 
                 || sInt. stop                               >= other. sInt. stop + Cds::peptideSize_min * a2s
                );
    }
private:
  string getGeneSymbol () const
    { return susceptible
               ? susceptible->genesymbol
               : alleleMatch () 
                 ? famId 
                 : nvl (checkPtr (getMatchFam ()) -> genesymbol, na);
    }
public:
  string fusion2geneSymbols () const
    { ASSERT (! isMutationProt ());
      if (fusions. empty ())
        return getGeneSymbol ();
      string s;
      for (const BlastAlignment* fusion : fusions)
        add (s, "/" /*fusion_infix*/, fusion->getGeneSymbol ());  // PD-5155 ??
      return s;
    }
  string fusion2famIds () const
    { if (isMutationProt () || fusions. empty ())
        return famId;
      string s;
      for (const BlastAlignment* fusion : fusions)
        add (s, fusion_infix, fusion->famId);
      return s;
    }
private:
  bool isCore () const 
    { return fusion2reportable () >= 2 || alleleReportable (); }
public:
  bool fusion2core () const
    { ASSERT (! isMutationProt ());
      if (fusions. empty ())
        return isCore ();
      for (const BlastAlignment* fusion : fusions)
        if (fusion->isCore ())
          return true;
      return false;
    }
private:
  string getType () const
    { return susceptible ? "AMR" : checkPtr (getMatchFam ()) -> type; }
public:
  string fusion2type () const
    { ASSERT (! isMutationProt ());
      if (fusions. empty ())
        return getType ();
      StringVector vec;
      for (const BlastAlignment* fusion : fusions)
        vec << fusion->getType ();
      vec. sort ();
      vec. uniq ();
      return vec. toString ("/");
    }    
private:
  string getSubtype () const
    { return susceptible ? "AMR" : checkPtr (getMatchFam ()) -> subtype; }
public:
  string fusion2subtype () const
    { ASSERT (! isMutationProt ());
      if (fusions. empty ())
        return getSubtype ();
      StringVector vec;
      for (const BlastAlignment* fusion : fusions)
        vec << fusion->getSubtype ();
      vec. sort ();
      vec. uniq ();
      return vec. toString ("/");
    }    
private:
	string getClass () const
	  { if (alleleMatch () && ! subclass. empty ())  // class may be empty()
	      return classS;
	    return susceptible ? susceptible->classS : checkPtr (getMatchFam ()) -> classS;
	  }
public:
  string fusion2class () const
    { ASSERT (! isMutationProt ());
      if (fusions. empty ())
        return getClass ();
      StringVector vec;
      for (const BlastAlignment* fusion : fusions)
        vec << std::move (StringVector (fusion->getClass (), '/', true));
      vec. sort ();
      vec. uniq ();
      return vec. toString ("/");
    }   
private: 
	string getSubclass () const
	  { if (alleleMatch () && ! subclass. empty ())
	      return subclass;
	    return susceptible ? susceptible->subclass : checkPtr (getMatchFam ()) -> subclass;
	  }
public:
  string fusion2subclass () const
    { ASSERT (! isMutationProt ());
      if (fusions. empty ())
        return getSubclass ();
      StringVector vec;
      for (const BlastAlignment* fusion : fusions)
        vec << std::move (StringVector (fusion->getSubclass (), '/', true));
      vec. sort ();
      vec. uniq ();
      return vec. toString ("/");
    }    
  const HmmAlignment* fusion2hmmAl () const
    { if (fusions. empty ())
        return hmmAl;
      for (const BlastAlignment* fusion : fusions)
        if (fusion->hmmAl)
          return fusion->hmmAl;
      return nullptr;
    }    
	string getMethod (const Locus &cds) const
	  { //IMPLY (refExactlyMatched () && ! mutation_all. get (), ! isMutationProt ())
	    string method (fromHmm
	                     ? "HMM"
	                     : isMutationProt () || isStrongSusceptibleProt ()
    	                   ? "POINT"
    	                   : refProtExactlyMatched (true) 
            	             ? alleleMatch () 
            	               ? "ALLELE"
            	               : "EXACT"  // PD-776
        	                 : partial ()
          	                 ? truncated (cds)
          	                   ? partialContigEnd_Name  // PD-2267
          	                   : partial_Name
        	                   : "BLAST"
        	           );
      // PD-2088, PD-2320
      bool suffix = true;
	    if (   (   method == "BLAST" 
    	        || method == partial_Name
    	        || method == partialContigEnd_Name
    	       ) 
	      //&& ! sProt  // Redundant
	       )
	    {
	      if (sInternalStop)
	      {
	        method = internalStop_Name;	
  	      suffix = false;
	      }
	    #if 0
	      else if (findDisruption (Disruption::eFrameshift))
	      {
	        method = frameshift_Name;
	        suffix = false;
	      }
	      else if (hasLongDisruption (10))  // PAR 
	      {
	        method = "DISRUPTED";
	        suffix = false;
	      }
	    #endif
	    }
	    if (suffix && method != "HMM")
	      method += (sProt ? "P" : "X");	  
	    return method;
	  }
	  // PD-736
	bool passBlastRule (const BlastRule &br) const
	  { ASSERT (! br. empty ());
	    return    relIdentity ()  >= br. ident        - frac_delta
  	         && qRelCoverage () >= br. ref_coverage - frac_delta;
	  }	
	bool partialPseudo () const
	  { return    partial () 
             && ! cdss. empty ()
             && ! truncatedCds ();
	  }
#if 0
	bool pseudo () const
	  { return    sInternalStop
	           || partialPseudo ();
	  }
#endif
	StringVector getContigs () const
	  { ASSERT (! cdss. empty ());
	    StringVector contigs;  contigs. reserve (cdss. size ());
	    for (const Locus& locus : cdss)
	      contigs << locus. contig;
	    contigs. sort ();
	    contigs. uniq ();
	    return contigs;
	  }
  bool good () const
    { ASSERT (! fromHmm);
      if (refAccession. empty ())
        return true;
    //if (! refMutation. empty ())   // PD-4981: drop
      //return true;      
    #if 0
      if (refAccession == "WP_104009840.1" && sseqid == "blaPDC-114_added")
      {
      //PRINT (reportPseudo);
      //PRINT (pseudo ());
        PRINT (isMutationProt ());
        PRINT (seqChanges. empty ());
        PRINT (isSusceptibleProt ());
        PRINT (susceptible);
        PRINT (disrs. empty ());
        PRINT (partial ());
        PRINT (parts);
        PRINT (truncatedCds ());
        PRINT (brFam);
        PRINT (inFam ());
      //PRINT (passBlastRule (partialBR));
      //PRINT (passBlastRule (completeBR));
        cout << endl;
        saveText (cout);        
      }
    #endif
    #if 0
      if (! reportPseudo && pseudo ())
        return false;
    #endif
      if (isMutationProt () && seqChanges. empty ())  // PD-4981
        return false;
      if (isSusceptibleProt () && ! susceptible)  // Alien organism
        return false;
		  if (susceptible)
		  {
        if (susceptible->cutoff)
        {
          if (relIdentity () > susceptible->cutoff + frac_delta) 
            return false;
        }
        else if (disrs. empty ())
          return false;
	    }
  	  // PD-1032
	    if (   partial ()
	        && (   parts > 1 
    	        || (! truncatedCds () && qInt. len () <= 35)  // PAR  // PD-3287
    	       )
    	   )
  	    return false;
  	  if (brFam)
  	    return true;  	  
  	  if (inFam ())  // !brFam
  	    return false;
	    if (partial ())
        return passBlastRule (partialBR);
      return passBlastRule (completeBR);
    }
private:
  bool insideEq (const BlastAlignment &other) const
    { ASSERT (sProt  == other. sProt);
      ASSERT (sseqid == other. sseqid);
    	if (sInt. strand != other. sInt. strand)
    	  return false;
    	if (   sInt. start + mismatchTail_aa * a2s >= other. sInt. start 
          && sInt. stop                          <= other. sInt. stop + mismatchTail_aa * a2s
         )
        return true;
      if (   ! partialPseudo ()
          || isMutationProt ()
          || other. isMutationProt ()
          || fusion2geneSymbols () != other. fusion2geneSymbols ()
         )
        return false;
      // PD-4698
      // Most probably a repeat protein
      const size_t pseudoOverlap = mismatchTail_aa * a2s * 2;  // PAR
      return    (sInt. start < other. sInt. start && sInt. stop                  >= other. sInt. start + pseudoOverlap)
             || (sInt. stop  > other. sInt. stop  && sInt. start + pseudoOverlap <= other. sInt. stop);
    }
  bool matchesCds (const BlastAlignment &other) const
    { ASSERT (sProt);
    	ASSERT (! other. sProt);
    	ASSERT (! cdss. empty ());
    	for (const Locus& cds : cdss)
    		if (   cds. contig == other. sseqid
    			  && cds. strand == (other. sInt. strand == 1)
    			  && ! cds. crossOrigin
    			 )
    		{ size_t protStart = cds. start;
    			size_t protStop  = cds. stop;
    			if (cds. strand)
    			{
    				ASSERT (protStop > 3);
    				protStop -= 3;
    			}
    			else
    			  protStart += 3;
    			ASSERT (protStart < protStop);
    			const size_t dnaStart  = other. sInt. start;
    			const size_t dnaStop   = other. sInt. stop;
    			ASSERT (dnaStart < dnaStop);
    			const size_t intersectionStart = max (protStart, dnaStart);
    			const size_t intersectionStop  = min (protStop,  dnaStop);
    			const size_t unionStart        = min (protStart, dnaStart);
    			const size_t unionStop         = max (protStop,  dnaStop);
    			if (   (   intersectionStart < intersectionStop
    				      && double (intersectionStop - intersectionStart) / double (protStop - protStart) > 0.75  // PAR, PD-2320
    				     )
    				  || (intersectionStart - unionStart) + (unionStop - intersectionStop) <= 60 * 3  // PAR, PD-4169
    				  || (   protStart <= dnaStart + mismatchTail_aa * 3 
    				      && dnaStop   <= protStop + mismatchTail_aa * 3 
    				      && other. partial ()
    				     )
    				  || (   dnaStart <= protStart + mismatchTail_aa * 3 
    				      && protStop <= dnaStop   + mismatchTail_aa * 3 
    				      && partial ()
    				     )
    				 )
    				return true;
    		}
    	return false;
    }
  bool betterEq (const BlastAlignment &other) const
    // Reflexive
    // Must: transitive
    { if (this == & other)
        return true;
      // PD-4981
      if (isMutationProt () != other. isMutationProt ())  
        return false;
      if (refMutation. empty () != other. refMutation. empty ())  
        return false;
      if (isSusceptibleProt () != other. isSusceptibleProt ())          
        return false;
      //
      if (sProt == other. sProt)  
      {
	    	if (sseqid != other. sseqid)
	        return false;
	    #if 0
	      if (sseqid == "WP_063839878.1.last91" && qseqid == "WP_063839878.1")
	      {
	        PRINT (other. insideEq (*this));
        	PRINT (insideEq (other));
        	PRINT (sProt);
        	PRINT (inFam ());
        	PRINT (other. inFam ());
        	PRINT (fusion2geneSymbols ());
        	PRINT (other. fusion2geneSymbols ());
        	PRINT (refAccession);
        	PRINT (other. refAccession);
        	//
        	PRINT (nident);
        	PRINT (sInt. start);
        	PRINT (sInt. stop);
        	PRINT (qInt. start);
        	PRINT (qInt. stop);
        	PRINT (hasMutation ());
        	PRINT (refProtExactlyMatched (false)); 
  	      PRINT (refEffectiveLen ());
	        //
        	PRINT (other. nident);
        	PRINT (other. sInt. start);
        	PRINT (other. sInt. stop);
        	PRINT (other. qInt. start);
        	PRINT (other. qInt. stop);
        	PRINT (other. hasMutation ());
        	PRINT (other. refProtExactlyMatched (false)); 
  	      PRINT (other. refEffectiveLen ());
  	      cout << endl;
	      }
	    #endif
        // PD-807, PD-4277
        if (   ! other. insideEq (*this)
        	  && !        insideEq (other)
        	 )
          if (   ! sProt 
              || (   inFam () /*! isMutationProt ()*/   // PD-4722
                  && other. inFam () /*! other. isMutationProt ()*/   // PD-4755
                  && fusion2geneSymbols () != other. fusion2geneSymbols ()
                 )
             )  // PD-4687
            return false;
        if (   inFam () /*! isMutationProt ()*/
            && ! refAccession. empty () 
            && refAccession == other. refAccession  // PD-4013
           )  
        {
  	      LESS_PART (other, *this, nident);
  	      LESS_PART (*this, other, sInt. start);
  	      LESS_PART (other, *this, sInt. stop);
  	      LESS_PART (*this, other, qInt. start);
  	      LESS_PART (other, *this, qInt. stop);
        }
        else
        {
  	      LESS_PART (other, *this, /*notInFam ()*/ hasMutation ());
  	      // PD-5427
  	      if (fusionOverrides (other))
  	        return true;
  	      if (other. fusionOverrides (*this))
  	        return false;
  	      //
  	      LESS_PART (other, *this, refProtExactlyMatched (false));  // PD-1261, PD-1678
  	      LESS_PART (other, *this, nident);
  	      LESS_PART (*this, other, refEffectiveLen ());  
  	    }
	    }
	    else
	    { 
        // PD-1902, PD-2139, PD-2313, PD-2320
	    	if (sProt && ! matchesCds (other))
	    	  return false;
	    	if (! sProt && ! other. matchesCds (*this))
	    	  return false;
	      LESS_PART (other, *this, refProtExactlyMatched (false));  
	    //LESS_PART (other, *this, allele ());  // PD-2352
	      LESS_PART (other, *this, alleleMatch ());  
	    //LESS_PART (*this, other, partial ());  // sProt => supposed to be a correct annotation
	    }
      if (isMutationProt () && other. isMutationProt ())  
      {
			  const Set<string> mutationSymbols      (       getMutationSymbols ());
			  const Set<string> otherMutationSymbols (other. getMutationSymbols ());
			  if (mutationSymbols == otherMutationSymbols && sProt != other. sProt)
			    return sProt;
			  if (   mutationSymbols. containsAll (otherMutationSymbols)
			      && ! otherMutationSymbols. containsAll (mutationSymbols)
			     )
			    return true;
			  if (   otherMutationSymbols. containsAll (mutationSymbols)
			      && ! mutationSymbols. containsAll (otherMutationSymbols)
			     )
			    return false;
			}
	    LESS_PART (other, *this, sProt);
      return true;
    }
  const Fam* getFam () const
    { ASSERT (inFam ());
      const Fam* fam = famId2fam [famId];
      if (! fam)
        fam = famId2fam [gene];
      if (! fam)
      	throw runtime_error ("Cannot find hierarchy for: " + famId + " (genesymbol: " + gene + ")");
      return fam;
    }
public:
  const Fam* getMatchFam () const
    { ASSERT (inFam ());
      if (fromHmm)
        return getFam ();
      return brFam;
    }
  bool better (const BlastAlignment &other) const
    { return    betterEq (other) 
    	       && (   ! other. betterEq (*this)
    	           || (! isMutationProt () /*inFam () --PD-5014*/ && ! equidistant && refAccession < other. refAccession)  // Tie resolution: PD-1245
    	          );
    }
  bool better (const HmmAlignment& other) const
    { ASSERT (other. good ());
    	ASSERT (other. blastAl. get ());
    	ASSERT (inFam ());
    	if (sProt)
    	{ if (sseqid != other. sseqid)
	        return false;
	    }
	    else
	    	if (! other. blastAl->matchesCds (*this))  // ??
	    		return false;
      return    refProtExactlyMatched (true) 
             || (checkPtr (getMatchFam ()) -> descendantOf (other. fam))  
             ;
    }
  size_t getCdsStart () const
    { if (cdss. empty ())
        return 0;
      size_t start = numeric_limits<size_t>::max ();
      for (const Locus& cds : cdss)
        minimize (start, cds. start);
      ASSERT (start != numeric_limits<size_t>::max ());
      return start;
    }
  size_t getCdsStop () const
    { if (cdss. empty ())
        return 0;
      size_t stop = 0;
      for (const Locus& cds : cdss)
        maximize (stop, cds. stop);
      ASSERT (stop);
      return stop;
    }
  void setCdss (const Annot &annot)
    { ASSERT (sProt);
    	ASSERT (cdss. empty ());
    	const Set<Locus> &cdss_ = annot. findLoci (sseqid);
    	ASSERT (! cdss_. empty ());
  	  insertAll (cdss, cdss_);
  	  qc ();
    }
  static bool less (const BlastAlignment* a,
                    const BlastAlignment* b) 
    { ASSERT (a);
      ASSERT (b);
      LESS_PART (*a, *b, sseqid);
      LESS_PART (*a, *b, sInt. start);
      LESS_PART (*a, *b, sInt. stop);
      LESS_PART (*a, *b, getCdsStart ());
      LESS_PART (*a, *b, getCdsStop ());
      LESS_PART (*a, *b, refAccession);
      LESS_PART (*a, *b, part);
    //LESS_PART (*a, *b, famId);
    //LESS_PART (*b, *a, relIdentity ());
      return false;
    }
  bool sameMatch (const BlastAlignment* other) const
    // Cf. less()
    { ASSERT (other);
      return    sseqid         == other->sseqid
             && sInt           == other->sInt
             && getCdsStart () == other->getCdsStart ()
             && getCdsStop ()  == other->getCdsStop ()
             && refAccession   == other->refAccession;
    }
};




bool HmmAlignment::better (const BlastAlignment& other) const
{ 
  ASSERT (good ());
  ASSERT (other. good ());
	ASSERT (other. inFam ());
	if (! other. sProt)
	  return false;
	if (sseqid != other. sseqid)
    return false;
  return    fam != other. getMatchFam ()  
         && fam->descendantOf (other. getMatchFam ());  
}




// Batch

struct Batch
{
  // Reference input
  map<string/*hmm*/,const Fam*> hmm2fam;
  uchar reportable_min {0};
  StringVector suppress_prots;  // of accessions
  StringVector alien_prots;  // of accessions

  // Target input
  VectorOwn<BlastAlignment> blastAls;
  VectorOwn<HmmAlignment> hmmAls;
  bool hmmExist {false};
  map<HmmAlignment::Pair, HmmAlignment::Domain> domains;  // Best domain  
  
  // Output
  //  targetProt => accession is protein 
  // !targetProt => accession is DNA 
  map<string/*accession*/,VectorPtr<BlastAlignment>> target2blastAls;
  map<string/*accession*/,VectorPtr<BlastAlignment>> target2goodBlastAls;
  //
  map<string/*protein accession*/,VectorPtr<HmmAlignment>> target2hmmAls;
  map<string/*protein accession*/,VectorPtr<HmmAlignment>> target2goodHmmAls;
    
  
  Batch (const string &famFName,
         const string &organism, 
         const string &mutation_tab,
         const string &susceptible_tab,
         const string &suppress_prot_FName,
         bool non_reportable,
         bool report_core_only)
    : reportable_min (non_reportable 
                        ? 0 
                        : report_core_only
                          ? 2
                          : 1
                     ) 
	  {
	  	const Chronometer_OnePass cop ("Batch", cerr, false, Chronometer::enabled);  
  
	    if (famFName. empty ())
	    	throw runtime_error ("fam (protein family hierarchy) file is missing");
	    	
	  	// Tree of Fam
	  	// Pass 1  
	    {
	    	if (verbose ())
	    		section ("Reading " + famFName + " Pass 1", true);
	      LineInput f (famFName);  
	  	  while (f. nextLine ())
	  	    try
  	  	  {
  	  	    trim (f. line);
  	  	    if (   f. line. empty () 
  	  	        || f. line [0] == '#'
  	  	       )
  	  	      continue;
  	      //cout << f. line << endl; 
  	  	    const string famId               (findSplit (f. line, '\t'));
  	  	    const string parentFamId         (findSplit (f. line, '\t'));
  	  	    const string genesymbol          (findSplit (f. line, '\t'));
  	  	    const string hmm                 (findSplit (f. line, '\t'));
  	  	    const double tc1 = str2<double>  (findSplit (f. line, '\t'));
  	  	    const double tc2 = str2<double>  (findSplit (f. line, '\t'));
  	  	    BlastRule completeBR;
  	  	    BlastRule partialBR;
  	  	    {
              completeBR. ident           = str2<double> (findSplit (f. line, '\t')); 
              completeBR. target_coverage = str2<double> (findSplit (f. line, '\t')); 
              completeBR. ref_coverage    = str2<double> (findSplit (f. line, '\t')); 
              partialBR.  ident           = str2<double> (findSplit (f. line, '\t')); 
              partialBR.  target_coverage = str2<double> (findSplit (f. line, '\t')); 
              partialBR.  ref_coverage    = str2<double> (findSplit (f. line, '\t')); 
      		    toProb (completeBR. ident);
      		    toProb (completeBR. target_coverage);
      		    toProb (completeBR. ref_coverage);
      		    toProb (partialBR.  ident);
      		    toProb (partialBR.  target_coverage);
      		    toProb (partialBR.  ref_coverage);
      		    QC_ASSERT (completeBR. empty () == partialBR. empty ());
      		    if (completeBR. empty ())
      		    {
        		    completeBR = defaultCompleteBR;
        		    partialBR  = defaultPartialBR;
      		    }
      		    else
      		    {
        		    completeBR. ref_coverage = defaultCompleteBR. ref_coverage;
        		    partialBR.  ref_coverage = defaultPartialBR.  ref_coverage;
        		    if (ident_min_user)
        		    {
        		      completeBR. ident = defaultCompleteBR. ident;
        		      partialBR.  ident = defaultPartialBR.  ident;
        		    }
        		  }
      		  }
  	  	    const uchar reportable = (uchar) str2<int> (findSplit (f. line, '\t'));
  	  	    const string type     (findSplit (f. line, '\t'));
  	  	    const string subtype  (findSplit (f. line, '\t'));
  	  	    const string classS   (findSplit (f. line, '\t'));
  	  	    const string subclass (findSplit (f. line, '\t'));
  	  	    const auto fam = new Fam (famId, genesymbol, hmm, tc1, tc2, completeBR, partialBR, type, subtype, classS, subclass, f. line, reportable);
  	  	    if (famId2fam [famId])
  	  	      throw runtime_error ("Family " + famId + " is duplicated");
  	  	    famId2fam [famId] = fam;
  	  	    if (! fam->hmm. empty ())
  	  	      hmm2fam [fam->hmm] = fam;
  	  	  }
  	  	  catch (const exception &e)
  	  	  {
  	  	    throw runtime_error ("Cannot read " + famFName +", " + f. lineStr () + "\n" + e. what ());
  	  	  }
	  	}
	    {
	    	if (verbose ())
	    		section ("Reading " + famFName + " Pass 2", true);
	      LineInput f (famFName);  
	  	  while (f. nextLine ())
	  	  {
	  	    trim (f. line);
	  	  //cout << f. line << endl;  
	  	    if (   f. line. empty () 
	  	        || f. line [0] == '#'
	  	       )
	  	      continue;
	  	    const Fam* child = famId2fam [findSplit (f. line, '\t')];
	  	    const string parentFamId (findSplit (f. line, '\t'));
	  	    QC_ASSERT (child);
	  	    const Fam* parent = nullptr;
	  	    if (! parentFamId. empty ())
	  	    { 
	  	    	parent = famId2fam [parentFamId]; 
	  	    	if (! parent)
	  	    	  throw runtime_error ("parentFamId " + strQuote (parentFamId) + " is not found in famId2fam for child " + strQuote (child->id));
	  	    }
	  	    var_cast (child) -> parent = parent;
	  	  }
	  	}
	  	
	  	if (qc_on)
	  	{
	  	  size_t roots = 0;
	  	  for (const auto& it : famId2fam)
	  	    if (! it. second->parent)
	  	      roots++;
	  	  QC_ASSERT (roots == 1);
	  	}
	  	
	    if (! organism. empty ())
	    {
	      {
	        if (mutation_tab. empty ())
	          throw runtime_error ("mutation_tab is empty");
  	    	if (verbose ())
  	    		cout << "Reading " << mutation_tab << endl;
  	      LineInput f (mutation_tab);
  	      Istringstream iss;
  	  	  while (f. nextLine ())
  	  	  {
  	  	    if (isLeft (f. line, "#"))
  	  	      continue;
  	  	    try
  	  	    {
    	  	  	string organism_, accession, geneMutation_std, geneMutation_report, classS, subclass, name;
    					int pos;
        	  	iss. reset (f. line);
    	  	  	iss >> organism_ >> accession >> pos >> geneMutation_std >> geneMutation_report >> classS >> subclass >> name;
    	  	  	QC_ASSERT (pos > 0);
    	  	  	QC_ASSERT (! name. empty ());
    	  	  	replace (organism_, '_', ' ');
    	  	  	if (organism_ == organism)
    	  	  		accession2mutations [accession]. emplace_back ((size_t) pos, geneMutation_std, geneMutation_report, classS, subclass, name);
    	  	  	else
    	  	  	  alien_prots << accession;
    	  	  }
    	  	  catch (const exception &e)
    	  	  {
    	  	    throw runtime_error ("Reading " + strQuote (mutation_tab) + " line:\n" + f. line + "\n" + e. what ());
    	  	  }
  	  	  }
  	  	  for (auto& it : accession2mutations)
  	  	  {
     	  	  it. second. sort ();
  	  	    if (! it. second. isUniq ())
  	  	  	  throw runtime_error ("Duplicate mutations for " + it. first);
  	  	  }
  	  	}
  	  	if (! susceptible_tab. empty ())
	      {
  	    	if (verbose ())
  	    		cout << "Reading " << susceptible_tab << endl;
  	      LineInput f (susceptible_tab);
  	      Istringstream iss;
  	  	  while (f. nextLine ())
  	  	  {
  	  	    if (isLeft (f. line, "#"))
  	  	      continue;
  	  	    try
  	  	    {
    	  	  	string organism_, genesymbol, accession, classS, subclass, name;
    					double cutoff = 0.0;
        	  	iss. reset (f. line);
    	  	  	iss >> organism_ >> genesymbol >> accession >> cutoff >> classS >> subclass >> name;
    	  	  	QC_ASSERT (! name. empty ());
    	  	  	replace (organism_, '_', ' ');
    	  	  	if (organism_ == organism)
    	  	  	{
    	  	  	  if (contains (accession2susceptible, accession))
    	  	  	    throw runtime_error ("Duplicate protein accession " + accession + " in " + susceptible_tab);
    	  	      accession2susceptible [accession] = std::move (Susceptible (genesymbol, cutoff, classS, subclass, name));
    	  	    }
    	  	  	else
    	  	  	  alien_prots << accession;
    	  	  }
    	  	  catch (const exception &e)
    	  	  {
    	  	    throw runtime_error ("Reading " + strQuote (susceptible_tab) + " line:\n" + f. line + "\n" + e. what ());
    	  	  }
  	  	  }
  	  	  if (verbose ())
  	  	    PRINT (accession2susceptible. size ());
  	  	}
	    }
	    alien_prots. sort ();
	  	  
  	  if (! suppress_prot_FName. empty ())
  	  {
  	    IFStream f (suppress_prot_FName);  	  	    
	  	  while (! f. eof ())
	  	  {
	  	    string accver;
	  	    f >> accver;
	  	    if (accver. empty ())
	  	      break;
	  	    suppress_prots << accver;
  	    }
  	  }  	  
  	  suppress_prots. sort ();
	  }


private:
  static void toProb (double &x)
  { 
    QC_ASSERT (x >= 0.0);
    QC_ASSERT (x <= 100.0);
    x /= 100.0;
  }
    
    
  void blastParetoBetter ()
  // Input: target2blastAls
  // Output: target2goodBlastAls
  {
    ASSERT (target2goodBlastAls. empty ());
	  for (const auto& it : target2blastAls)
	  {
  	  VectorPtr<BlastAlignment>& vec = target2goodBlastAls [it. first];
  	  for (const BlastAlignment* blastAl : it. second)
      {
        ASSERT (blastAl);
      	ASSERT (blastAl->good ());
    	  bool found = false;
    	  for (const BlastAlignment* goodBlastAl : vec)
    	    if (goodBlastAl->better (*blastAl))
  	      {
  	        found = true;
  	        break;
  	      }
  	    if (found)
  	      continue;	      
        for (Iter<VectorPtr<BlastAlignment>> goodIter (vec); goodIter. next ();)
          if (blastAl->better (**goodIter))
            goodIter. erase ();       
        ASSERT (! blastAl->sseqid. empty ());
        vec << blastAl;
      }
    }
    reportDebug ("Pareto-better");
  }


  void setStopCodon (BlastAlignment &blastAlP)  
  { 
    ASSERT (blastAlP. sProt);
    for (const Locus& cds : blastAlP. cdss)
      if (const VectorPtr<BlastAlignment>* blastAlsX = findPtr (target2blastAls, cds. contig))
        for (const BlastAlignment* blastAlX : *blastAlsX)
          if (   blastAlX->blastx ()
              && blastAlX->sInternalStop
              && blastAlP. better (*blastAlX)
             )
            blastAlP. sInternalStop = true;
  }
    
    
#if 0
  struct Frameshift 
  {
    string refName;
    size_t targetPos {no_index};
    size_t refPos {no_index};
    long insertion {0};
      // %3
  };


  static bool frameshiftSort (const BlastAlignment* al1,
                              const BlastAlignment* al2)
  {
    ASSERT (al1);
    ASSERT (al2);
    LESS_PART (*al1, *al2, sseqid);
    LESS_PART (*al1, *al2, qseqid);
    LESS_PART (*al1, *al2, sInt. strand);
    LESS_PART (*al1, *al2, sInt. start);
    return false;
  }
#endif
                         

  VectorPtr<BlastAlignment> processDisruptions (const VectorPtr<BlastAlignment> &origAls)
  // Return: new 
  {
    ASSERT (! origAls. empty ());
      
    Unverbose unv;

    VectorPtr<BlastAlignment> als;  

    VectorPtr<Hsp> origHsps;  origHsps. reserve (origAls. size ());
    for (const BlastAlignment* al : origAls)
    {
      ASSERT (al);
      if (   al->blastx ()
          && ! al->isMutationProt ()
          && al->isStrongSusceptibleProt ()  
         )
        origHsps << al;
      else
        als << al;
    }

    Hsp::Merge merge (origHsps, nullptr/*sm*/, 20, true/*bacteria*/);  // PAR
    for (;;)
    { 
      const Hsp* origHsp = nullptr;
      AlignScore score = - score_inf;
      Hsp hsp (merge. get (origHsp, score)); 
      if (hsp. empty ()) 
        break; 
      ASSERT (origHsp);
      auto al = new BlastAlignment (* static_cast <const BlastAlignment*> (origHsp));
      ASSERT (al->refMutation. empty ());
      ASSERT (al->seqChanges. empty ());
      blastAls << al;
      * static_cast <Hsp*> (al) = std::move (hsp);
      al->qc ();
      als << al;
    }

    return als;
  }
  
    
public:
	void process (bool retainBlasts,
	              bool skip_hmm_check) 
  // Input: target2blastAls, domains, target2hmmAls
	// Output: target2goodBlastAls
	{
    ASSERT (target2goodBlastAls. empty ());
    ASSERT (target2goodHmmAls. empty ());

		const Chronometer_OnePass cop ("process", cerr, false, Chronometer::enabled);  

		  		  
		// Disruption's
    for (auto& it : target2blastAls)
    {
      VectorPtr<BlastAlignment>& als = it. second;
      ASSERT (! als. empty ());
      if (qc_on)
      {
        als. sort ();  
        ASSERT (als. isUniq ());  
      }
      als. sort (Hsp::less); 
      VectorPtr<BlastAlignment> mergedAls;  
      {  
        VectorPtr<BlastAlignment> hsps;  // Subset of als
        const BlastAlignment* prev = nullptr;
    	  for (const BlastAlignment* al : als)
    	  {
    	    ASSERT (al);
    	  //ASSERT (al->sseqid == it. first);  // Locus::contig != al->sseqid
          if (   prev
              && ! (   al->sseqid       == prev->sseqid
                    && al->sInt. strand == prev->sInt. strand
                    && al->qseqid       == prev->qseqid
                   )
             )
    	    {
    	      mergedAls << processDisruptions (hsps);
    	      hsps. clear ();
    	    }
 	        hsps << al;
          prev = al;     
    	  }
    	  ASSERT (prev);
        mergedAls << processDisruptions (hsps);
      }
      ASSERT (! mergedAls. empty ());
      als = mergedAls;      
      ASSERT (! als. empty ());
      for (const BlastAlignment* al : als)
        if (al->isStrongSusceptibleProt ())  
          for (const Disruption& disr : al->disrs)
          {
            var_cast (al) -> seqChanges << std::move (SeqChange (al, disr));
            break;  // PD-5394
          }
    }
 	  reportDebug ("Unframeshifted Blasts");
 	  
 	  
    for (auto& it : target2blastAls)
    {
      VectorPtr<BlastAlignment>& als = it. second;
      for (const BlastAlignment* al : als)
      {
        ASSERT (al);        
        if (   ! al->fromHmm 
            && ! al->sProt
           )
          var_cast (al) -> finish ();
  	    QC_ASSERT ((al->fromHmm || al->inFam ()) == al->completeBR. empty ());
  	    QC_ASSERT ((al->fromHmm || al->inFam ()) == al->partialBR.  empty ());
      }
    }
 	  
 	  
 	#if 0
    // Declarative frameshifts for mutation proteins  // --> processed as susceptible 
    for (auto& it : target2blastAls)
    {
      VectorPtr<BlastAlignment>& als = it. second;
      ASSERT (! als. empty ());
      als. sort (frameshiftSort);      
      
      Vector<Frameshift> frameshifts;
      {
        // Cf. tblastn2frameshift.cpp
        constexpr size_t diff_max_aa = 30;  // PAR
        const BlastAlignment* al_prev = nullptr;
    	  for (const BlastAlignment* al : als)
    	  {
    	    ASSERT (al);
    	  //ASSERT (al->sseqid == it. first);  // Locus::contig != al->sseqid
    	    if (al->hasDeclarativeFrameshift ())
    	      continue;
    	    if (   al_prev
    	        && al_prev->blastx ()
    	        && al->blastx ()
    	        && al_prev->sseqid       == al->sseqid
    	        && al_prev->qseqid       == al->qseqid
    	        && al_prev->sInt. strand == al->sInt. strand
    	        && difference (al_prev->sInt. stop, al->sInt. start) <= 3 * diff_max_aa
    	        && (   (al->sInt. strand ==  1 && difference (al_prev->qInt. stop,  al->qInt. start) <= diff_max_aa)
      	          || (al->sInt. strand == -1 && difference (al_prev->qInt. start, al->qInt. stop)  <= diff_max_aa)
      	         )
      	     )
          {
            const long diff = al->sStartGlobal () - al_prev->sStartGlobal ();
            if (diff % 3)  
              frameshifts << Frameshift { al->qseqid
                                        , al->sInt. strand == 1 ? al_prev->sInt. stop : (al->sInt. start - al->a2s)
                                        , al->sInt. strand == 1 ? al_prev->qInt. stop : al->qInt. stop
                                        , diff
                                        };
          }
    	    al_prev = al;
    	  }
    	}

  	  for (const BlastAlignment* al : als)
  	    if (al->hasDeclarativeFrameshift ())
    	  {
    	    ASSERT (al->seqChanges. size () == 1);
    	    const SeqChange& seqChange = al->seqChanges [0];
    	    ASSERT (seqChange. hasFrameshift ());
          const AmrMutation* mut = seqChange. mutations [0];
          ASSERT (mut);
          bool found = false;
          for (const Frameshift& fs : frameshifts)
            if (   al->qseqid                       == fs. refName
                && seqChange. start_target          == fs. targetPos
                && (size_t) mut->pos_std            == fs. refPos
                && (long) mut->frameshift_insertion == fs. insertion
               )
            {
              found = true;
              break;
            }
          if (! found)
            var_cast (al) -> seqChanges. clear ();
    	  }
    }
  #endif


    for (auto& it : target2blastAls)
      for (Iter<VectorPtr<BlastAlignment>> iter (it. second); iter. next ();)
        if (alien_prots. containsFast ((*iter)->refAccession))
          iter. erase ();
	  reportDebug ("Non-alien Blasts");


    for (auto& it : target2blastAls)
      for (Iter<VectorPtr<BlastAlignment>> iter (it. second); iter. next ();)
        if ((*iter)->good ())
        {
          ASSERT ((bool) (*iter)->susceptible == (*iter)->isSusceptibleProt ());
        }
        else
          iter. erase ();
 	  reportDebug ("Good Blasts");
        

	  // PD-2322
	  // setStopCodon()
    for (const auto& it : target2blastAls)
   	  for (const BlastAlignment* blastAlP : it. second)
        if (blastAlP->sProt)
          setStopCodon (* var_cast (blastAlP));
 	  for (const HmmAlignment* hmmAl : hmmAls)
 	  {
 	    const BlastAlignment* blastAl = hmmAl->blastAl. get ();
 	    ASSERT (blastAl);
 	    setStopCodon (* var_cast (blastAl));  
 	  }
 	  if (verbose (-1))
 	  {
 	    cout << "After setStopCodon():" << endl;
   	  for (const auto& it : target2blastAls)
  		  for (const BlastAlignment* blastAl : it. second)
  		  {
  		    blastAl->saveText (cout);
  		    cout << endl;
  		  }
  	}
 	  

    if (retainBlasts)
      for (const auto& it : target2blastAls)
        target2goodBlastAls [it. first] = it. second;
    else
      blastParetoBetter ();


    // Cf. dna_mutation.cpp
    for (const auto& it : target2goodBlastAls)
      for (const BlastAlignment* blastAl1 : it. second)
        for (const SeqChange& seqChange1 : blastAl1->seqChanges)
        {
          ASSERT (seqChange1. al == blastAl1);
        //ASSERT (seqChange1. mutation);
          for (const BlastAlignment* blastAl2 : it. second)
          {
            if (   blastAl2->sProt        == blastAl1->sProt
                && blastAl2->sseqid       == blastAl1->sseqid
                && blastAl2->sInt. strand == blastAl1->sInt. strand
                && blastAl2               != blastAl1
               )  
            //for (Iter<Vector<SeqChange>> iter (var_cast (blastAl2) -> seqChanges); iter. next (); )
              for (SeqChange& seqChange2 : var_cast (blastAl2) -> seqChanges)
              {
              //SeqChange& seqChange2 = *iter;
              //ASSERT (seqChange2. mutation);
                ASSERT (seqChange2. al == blastAl2);
                if (   seqChange1. start_target     == seqChange2. start_target 
                    && seqChange1. hasFrameshift () == seqChange2. hasFrameshift ()
                    && seqChange1. better (seqChange2)
                   )
                //iter. erase ();
                  seqChange2. replacement = & seqChange1;
              }
          }
        }    


    // HMM: Pareto-better()  
    // target2hmmAls --> target2goodHmmAls
    for (auto& it : target2hmmAls)
    {
      VectorPtr<HmmAlignment> hmmAls_ (it. second);
      VectorPtr<HmmAlignment>& goodHmmAls = target2goodHmmAls [it. first];
      goodHmmAls. reserve (hmmAls_. size ());  
      FOR (unsigned char, criterion, 2)
      {
        // hmmAls_ --> goodHmmAls
        ASSERT (goodHmmAls. empty ());
    	  for (const HmmAlignment* hmmAl : hmmAls_)
        {
          ASSERT (hmmAl);
        	ASSERT (hmmAl->good ());
      	  bool found = false;
      	  for (const HmmAlignment* goodHmmAl : goodHmmAls)
      	    if (goodHmmAl->better (*hmmAl, criterion))
    	      {
    	        found = true;
    	        break;
    	      }
    	    if (found)
    	      continue;  
          for (Iter<VectorPtr<HmmAlignment>> goodIter (goodHmmAls); goodIter. next ();)
            if (hmmAl->better (**goodIter, criterion))
              goodIter. erase ();
          goodHmmAls << hmmAl;
        }
        reportDebug ("Pareto-better HMMs (Criterion " + to_string ((int) criterion) + ")");
        if (criterion < 1)
          hmmAls_ = std::move (goodHmmAls);
      }
    }


    // PD-741
  	if (hmmExist && ! skip_hmm_check)
  	  for (auto& it : target2goodBlastAls)
        for (Iter<VectorPtr<BlastAlignment>> iter (it. second); iter. next ();)
          if (   (*iter) -> inFam ()
          	  && (*iter) -> sProt
              && ! (*iter) -> partial ()
          	  && (*iter) -> relIdentity () < 0.98 - frac_delta  // PAR  // PD-1673  
             )
  	        if (const Fam* fam = checkPtr ((*iter) -> getMatchFam ()) -> getHmmFam ())    
  	        {
  	          bool found = false;
  	      	  for (const HmmAlignment* hmmAl : target2goodHmmAls [it. first])
  	            if (   (*iter) -> sseqid == hmmAl->sseqid
  	                && fam == hmmAl->fam
  	               )
  	            {
  	              found = true;
  	              break;
  	            }
  	          if (! found)   // BLAST is wrong
  	          {
  	            if (verbose ())
  	            {
                  cout << "HMM-wrong:  ";
                  (*iter) -> saveText (cout); 
                  cout << endl;
                }
  	            iter. erase ();
  	          }
  	        }
    reportDebug ("Best Blasts left");


    // PD-2783
    for (auto& it : target2goodHmmAls)
      for (Iter<VectorPtr<HmmAlignment>> hmmIt (it. second); hmmIt. next ();)
    	  for (const BlastAlignment* blastAl : target2goodBlastAls [it. first])
    	    if (blastAl->inFam () && blastAl->better (**hmmIt))
  	      {
  	        var_cast (blastAl) -> hmmAl = *hmmIt;
            hmmIt. erase ();
  	        break;
  	      }
    reportDebug ("HMMs non-suppressed by BLAST");


 	  for (auto& it : target2goodBlastAls)
      for (Iter<VectorPtr<BlastAlignment>> blastIt (it. second); blastIt. next ();)
        if ((*blastIt) -> inFam ())
      	  for (const HmmAlignment* hmmAl : target2goodHmmAls [it. first])
      	    if (hmmAl->better (**blastIt))
    	      {
              blastIt. erase ();
    	        break;
    	      }
    reportDebug ("Best HMMs left");


    // Output 
    
    // goodHmmAls --> target2goodBlastAls
    for (const auto& it : target2goodHmmAls)
    	for (const HmmAlignment* hmmAl : it. second)
    	{
    	  const BlastAlignment* blastAl = hmmAl->blastAl. get ();
    	  ASSERT (blastAl);
    	  target2goodBlastAls [it. first] << blastAl;
    	}
  	
    // PD-2394
    // BlastAlignment::{fusions,fusionRedundant}
 	  for (auto& it : target2goodBlastAls)
 	  {
 	    auto& goodBlastAls_ = it. second;
      goodBlastAls_. sort (BlastAlignment::less);  
      {
        const BlastAlignment* prev = nullptr;
        const BlastAlignment* fusionMain = nullptr;
        for (const BlastAlignment* blastAl : goodBlastAls_)
        {
          if (! blastAl->inFam () /*blastAl->isMutationProt ()*/)
            continue;
          if (blastAl->fromHmm)
            continue;
          if (   prev 
              && prev->sameMatch (blastAl)
             )
          {
            QC_ASSERT (blastAl->parts == prev->parts);
            QC_ASSERT (blastAl->part  >= prev->part);
            if (blastAl->part == prev->part)  // multi-domain tccP problem, PD-4217
              fusionMain = nullptr;
            else
            {
              QC_ASSERT (blastAl->parts >= 2);
              if (! fusionMain)
              {
                fusionMain = prev;
                var_cast (fusionMain) -> fusions << fusionMain;
              }
              ASSERT (fusionMain);
              var_cast (fusionMain) -> fusions << blastAl;
              var_cast (blastAl) -> fusionRedundant = true;
            }
          }
          else
            fusionMain = nullptr;
          prev = blastAl;
        }
      }
    }

    reportDebug ("After process()");
	}
	
	
	void reportDebug (const string &header) const
  {
    if (! verbose ())
      return;
    
    cout << endl;
    
    cout << header << " (target2blastAls):" << endl;
 	  for (const auto& it : target2blastAls)
		  for (const BlastAlignment* blastAl : it. second)
		  {
		    cout << it. first << ": ";
		    blastAl->saveText (cout);
		    cout << endl;
		  }
		cout << endl;
		  
    cout << header << " (target2goodBlastAls):" << endl;
 	  for (const auto& it : target2goodBlastAls)
		  for (const BlastAlignment* blastAl : it. second)
		  {
		    cout << it. first << ": ";
		    blastAl->saveText (cout);
		    cout << endl;
		  }
		cout << endl;
		  
    cout << header << " (target2hmmAls):" << endl;
 	  for (const auto& it : target2hmmAls)
		  for (const HmmAlignment* hmmAl : it. second)
		  {
		    cout << it. first << ": ";
		    hmmAl->saveText (cout);
		  }
		cout << endl;

    cout << header << " (target2goodHmmAls):" << endl;
 	  for (const auto& it : target2goodHmmAls)
		  for (const HmmAlignment* hmmAl : it. second)
		  {
		    cout << it. first << ": ";
		    hmmAl->saveText (cout);
		  }
		cout << endl;

		cout << endl;
	}
		
	
	void report (TsvOut &td,
	             bool mutationAll) const
	// Input: target2goodBlastAls
	{
	  ASSERT (td. empty ());
	  	  
		const Chronometer_OnePass cop ("report", cerr, false, Chronometer::enabled);  

    // PD-283, PD-780
  	// Cf. BlastAlignment::report()
    if (! input_name. empty ())
      td << "Name";
    // Target
    td << prot_colName;                          //  1  // sseqid 
    if (cdsExist)  
      // Contig (target)
      td << contig_colName                       //  2
         << start_colName                        //  3
         << stop_colName                         //  4
         << strand_colName;                      //  5
    // Reference
    td << genesymbol_colName                     //  6 or 2  
       << elemName_colName                       //  7 or 3
       << scope_colName                          //  8 or 4       
       << type_colName                           //  9 or 5
       << subtype_colName                        // 10 or 6
       << class_colName                          // 11 or 7
       << subclass_colName                       // 12 or 8
       << method_colName                         // 13 or 9
       // target
       << targetLen_colName                      // 14 or 10       
       //
       << refLen_colName                         // 15 or 11  // qlen
       << refCov_colName                         // 16 or 12  // queryCoverage
       << refIdent_colName                       // 17 or 13 
       << alignLen_colName                       // 18 or 14  // sseq.size()
       << closestRefAccession_colName            // 19 or 15  // refAccession
       << closestRefName_colName                 // 20 or 16
       //
       << hmmAccession_colName                   // 21 or 17
       << hmmDescr_colName                       // 22 or 18
       ;
    if (cdsExist)
    	if (useCrossOrigin)
      	 td << "Cross-origin length";
    if (print_node)
      td << hierarchyNode_colName; 
    td. newLn ();

 	  for (const auto& it : target2goodBlastAls)
    	for (const BlastAlignment* blastAl : it. second)
    	{
    	  ASSERT (blastAl);
     	  blastAl->qc ();
   	  	if (   ! blastAl->seqChanges. empty ()
   	  	    || (   (   blastAl->fusion2reportable () >= reportable_min
     	              || blastAl->alleleReportable ()
     	             )
     	          && ! suppress_prots. containsFast (blastAl->refAccession)
     	          && ! blastAl->fusionRedundant
     	         )
     	     )
       	  blastAl->report (td, it. first, mutationAll);
      }
	}


	void printTargetIds (ostream &os) const
	{
		QC_ASSERT (os. good ());
 	  for (const auto& it : target2goodBlastAls)
    	for (const BlastAlignment* blastAl : it. second)
    	  if (   blastAl->sProt
    	  	  && ! blastAl->isMutationProt ()
    	  	  && blastAl->fusion2reportable () >= reportable_min
    	  	 )
          os << blastAl->sseqid << endl;
	}
};




// HmmAlignment

HmmAlignment::HmmAlignment (const string &line,
                            const Batch &batch)
{
  static Istringstream iss;
  iss. reset (line);
  string hmm, dummy;
  //                                                        --- full sequence ---  --- best 1 domain --
  //     target name  accession  query name     accession   E-value  score     bias     E-value  score  bias   exp reg clu  ov env dom rep inc description of target
  iss >> sseqid >>    dummy      >> dummy      >> hmm >>    dummy >> score1 >> dummy >> dummy >> score2;
  QC_ASSERT (score1 > 0);
  QC_ASSERT (score2 > 0)
  find (batch. hmm2fam, hmm, fam);
  if (! fam)
    throw runtime_error ("No family for HMM " + hmm);
}




// Domain

HmmAlignment::Domain::Domain (const string &line,
                              Batch &batch)
{
  static Istringstream iss;
  iss. reset (line);
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
  QC_ASSERT (accession == "-");
  QC_ASSERT (hmmStart);
  QC_ASSERT (seqStart);
  QC_ASSERT (env_from);
  hmmStart--;
  seqStart--;
  env_from--;
  QC_ASSERT (hmmStart < hmmStop);
  QC_ASSERT (seqStart < seqStop);
  QC_ASSERT (hmmStop <= hmmLen);
  QC_ASSERT (seqStop <= seqLen);
  QC_ASSERT (full_score > 0);
  QC_ASSERT (n >= 1);
  QC_ASSERT (n <= of);
  QC_ASSERT (score > 0);

  const Fam* fam = batch. hmm2fam [query_accession];
  if (! fam)
    return;
  const HmmAlignment::Pair p (target_name, fam->id);
  const HmmAlignment::Domain domain_old (batch. domains [p]);
  if (domain_old. score > score)
    return;
  batch. domains [p] = *this;
}




// ThisApplication

constexpr double ident_min_def = 0.9;
constexpr double complete_coverage_min_def = 0.9;
constexpr double partial_coverage_min_def  = 0.5;



struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Report AMR proteins")
    {
      // Input
      addKey ("fam", "Table FAM");
      const string blastFormat (string (Hsp::format [true]) + ". qseqid format: gi|Protein accession|fusion part|# fusions|FAM.id|FAM.class|Product name");
      addKey ("blastp", "blastp output in the format: " + blastFormat);  
      addKey ("blastx", "blastx output in the format: " + blastFormat);  

      addKey ("gff", ".gff assembly file");
      addKey ("gfftype", "Type of GFF file: " + Gff::names. toString (", "), "genbank");
      addKey ("gff_prot_match", ".gff-protein FASTA matching file: \"<protein FASTA id> <protein GFF id>\"");
      addKey ("gff_dna_match", ".gff-DNA FASTA matching file: \"<DNA FASTA id> <DNA GFF id>\"");
      addFlag ("lcl", "Nucleotide FASTA created by PGAP has \"lcl|\" prefix in accessions");  
      addFlag ("bed", "Browser Extensible Data format of the <gff> file");

      addKey ("dna_len", "File with lines: <dna id> <dna length>");
      addKey ("hmmdom", "HMM domain alignments");
      addKey ("hmmsearch", "Output of hmmsearch");
      addKey ("organism", "Taxonomy group for mutations");
      addKey ("mutation", "Mutations table");
      addKey ("susceptible", "Table of susceptible proteins with resistance cutoffs");
      addKey ("mutation_all", "File to report all mutations");
      addKey ("suppress_prot", "File with protein accessions to suppress");
      addKey ("ident_min", "Min. identity to the reference protein (0..1). -1 means use a curated threshold if it exists and " + toString (ident_min_def) + " otherwise", "-1");
      addKey ("coverage_min", "Min. coverage of the reference protein (0..1) for partial hits", toString (partial_coverage_min_def));
      addFlag ("skip_hmm_check", "Skip checking HMM for a BLAST hit");
      addFlag ("report_equidistant", "Report all equidistant BLAST and HMM matches");  // PD-3772
      
      // Output
      addKey ("out", "Identifiers of the reported input proteins");
      addFlag ("print_node", "Print FAM.id replaced by FAM.parent for non-exact allele matches"); 
      addFlag ("print_node_raw", "Print FAM.id"); 
    //addFlag ("pseudo", string ("Indicate pseudo-genes as method ") + internalStop_Name + " or " + frameshift_Name);
      addFlag ("force_cds_report", "Report contig/start/stop/strand even if this information does not exist");
      addFlag ("non_reportable", "Report non-reportable families");
      addFlag ("core", "Report only core reportale families");
      addKey ("name", "Text to be added as the first column \"name\" to all rows of the report");
      
      // Testing
      addFlag ("nosame", "Exclude the same reference protein accessions from the BLAST output (for testing)"); 
      addFlag ("noblast", "Exclude the BLAST output (for testing)"); 
      addFlag ("nohmm", "Exclude the HMMer output (for testing)"); 
      addFlag ("retain_blasts", "Retain all blast hits (for testing)");
      
 	    version = SVN_REV;  
    }
    
    

  void body () const final
  {
    static_assert (complete_coverage_min_def > partial_coverage_min_def, "complete_coverage_min_def > partial_coverage_min_def");
    
    const string  famFName             = getArg ("fam");
    const string  blastpFName          = getArg ("blastp");
    const string  blastxFName          = getArg ("blastx");
    const string  gffFName             = getArg ("gff");
    const Gff::Type gffType            = Gff::name2type (getArg ("gfftype"));
    const string  gffProtMatchFName    = getArg ("gff_prot_match");
    const string  gffDnaMatchFName     = getArg ("gff_dna_match");
    const bool    lcl                  = getFlag ("lcl");
    const bool    bedP                 = getFlag ("bed");
    const string  dnaLenFName          = getArg ("dna_len");
    const string  hmmDom               = getArg ("hmmdom");
    const string  hmmsearch            = getArg ("hmmsearch");  
          string  organism             = getArg ("organism");  
    const string  mutation_tab         = getArg ("mutation");  
    const string  susceptible_tab      = getArg ("susceptible");  
    const string  mutation_all_FName   = getArg ("mutation_all");
    const string  suppress_prot_FName  = getArg ("suppress_prot");
          double  ident_min            = str2<double> (getArg ("ident_min"));  
    const double  partial_coverage_min = str2<double> (getArg ("coverage_min"));  
    const bool    skip_hmm_check       = getFlag ("skip_hmm_check"); 
                  equidistant          = getFlag ("report_equidistant");
    const string  outFName             = getArg ("out");
                  print_node           = getFlag ("print_node");
                  print_node_raw       = getFlag ("print_node_raw");
                //reportPseudo         = getFlag ("pseudo");
    const bool    force_cds_report     = getFlag ("force_cds_report");
    const bool    non_reportable       = getFlag ("non_reportable");
    const bool    report_core_only     = getFlag ("core");
                  input_name           = getArg ("name");
    const bool    nosame               = getFlag ("nosame");
    const bool    noblast              = getFlag ("noblast");
    const bool    nohmm                = getFlag ("nohmm");
    const bool    retainBlasts         = getFlag ("retain_blasts");
    
    replace (organism, '_', ' ');
    
    QC_ASSERT (hmmsearch. empty () == hmmDom. empty ());
    QC_IMPLY (! outFName. empty (), ! blastpFName. empty ());
    QC_IMPLY (! gffFName. empty (), ! blastpFName. empty ());
    if (! blastpFName. empty () && ! blastxFName. empty () && gffFName. empty ())
    	throw runtime_error ("If BLASTP and BLASTX files are present then a GFF file must be present");
    	
    if (gffFName. empty ())
    {
      QC_ASSERT (gffProtMatchFName. empty ());
      QC_ASSERT (gffDnaMatchFName. empty ());
      QC_ASSERT (gffType == Gff::genbank);
      QC_ASSERT (! lcl);
      QC_ASSERT (! bedP);
    }
    
    QC_IMPLY (print_node_raw, print_node);
           			  
    
    // defaultCompleteBR, defaultPartialBR, ident_min_user
    {
      if (ident_min == -1.0)
        ident_min = ident_min_def;
      else
        ident_min_user = true;
        
      if (! (ident_min >= 0.0 && ident_min <= 1.0))
      	throw runtime_error ("-ident_min must be -1 or between 0 and 1");
      if (! (partial_coverage_min >= 0.0 && partial_coverage_min <= 1.0))
      	throw runtime_error ("-coverage_min must be -1 or between 0 and 1");    	
      if (partial_coverage_min > complete_coverage_min_def)
        throw runtime_error ("-coverage_min must be less than " + toString (complete_coverage_min_def) + " - threshod for complete matches");
        
      defaultCompleteBR = BlastRule (ident_min, complete_coverage_min_def);  
      defaultPartialBR  = BlastRule (ident_min, partial_coverage_min);
    }
    
    
    cdsExist =    force_cds_report
               || ! blastxFName. empty ()
               || ! gffFName. empty ();
    
    
		const Chronometer_OnePass cop1 ("amr_report", cerr, false, Chronometer::enabled);  
  

    Batch batch (famFName, organism, mutation_tab, susceptible_tab, suppress_prot_FName, non_reportable, report_core_only);  
  
  
    // Input 

    // batch.{blastAls,target2blastAls}
  	ASSERT (batch. blastAls. empty ());
  	ASSERT (batch. target2blastAls. empty ());
    if (   ! noblast
        && ! blastpFName. empty ()
       )
    {
  		const Chronometer_OnePass cop ("blastp", cerr, false, Chronometer::enabled);
      LineInput f (blastpFName);  
  	  while (f. nextLine ())
  	  {
  	    { 
  	      Unverbose unv;
  	      if (verbose ())
  	        cout << f. line << endl;  
  	    }
  	    unique_ptr<BlastAlignment> al (new BlastAlignment (f. line, true));
	      al->qc ();  
  	    if (nosame && al->refAccession == al->sseqid)
  	      continue;
        batch. blastAls << al. get ();
        ASSERT (! al->sseqid. empty ());
        batch. target2blastAls [al->sseqid] << al. release ();
  	  }
  	}
  	if (verbose ())
  	  cout << "# Blasts: " << batch. blastAls. size () << endl;
  	
  
    if (! nohmm)
    {
      // batch.domains
      if (! hmmDom. empty ())
      {
    		const Chronometer_OnePass cop ("hmmDom", cerr, false, Chronometer::enabled);  
      	batch. hmmExist = true;
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

      // batch.{hmmAls,target2hmmAls}
    	// HmmAlignment::good()
    	if (! hmmsearch. empty ())  
    	{
    		const Chronometer_OnePass cop ("hmmsearch", cerr, false, Chronometer::enabled); 
        LineInput f (hmmsearch);
    	  while (f. nextLine ())
    	  {
    	    if (verbose ())
    	      cout << f. line << endl;  
    	    if (f. line. empty () || f. line [0] == '#')
    	      continue;
    	    unique_ptr<HmmAlignment> hmmAl (new HmmAlignment (f. line, batch));
    	    if (! hmmAl->good ())
    	    {
    	      if (verbose ())
    	      {
    	        cout << "  Bad HMM: " << endl;
    	        hmmAl->saveText (cout);
    	      }
    	    	continue;
    	    }
   	      const BlastAlignment* bestBlastAl = nullptr;  // PD-3475
    	    if (const VectorPtr<BlastAlignment>* blastAls_ = findPtr (batch. target2blastAls, hmmAl->sseqid))
    	    {
    	      size_t nident_max = 0;
    	      for (const BlastAlignment* blastAl : *blastAls_)
    	      {
    	        ASSERT (blastAl->sseqid == hmmAl->sseqid);
    	        if (maximize (nident_max, blastAl->nident))
    	          bestBlastAl = blastAl;
    	      }
    	    }
    	    auto al = new BlastAlignment (*hmmAl, bestBlastAl);
  	  	  hmmAl->blastAl. reset (al);
  	  	  if (verbose ())
  	  	    cout << al->sseqid << " " << al->gene << endl;  
  	  	  const HmmAlignment::Domain domain = batch. domains [HmmAlignment::Pair (al->sseqid, al->gene)];
  	  	  if (! domain. hmmLen)  
  	  	    continue;  // domain does not exist
 	  	    if (! bestBlastAl)  // Stand-alone HMM hit
  	  	  {
    	  	/*al->qlen        = domain. hmmLen;
    	  	  al->qstart      = domain. hmmStart;
    	  	  al->qend        = domain. hmmStop; */
    	  	  al->slen        = domain. seqLen;
    	  	  al->sInt. start = domain. seqStart;
    	  	  al->sInt. stop  = domain. seqStop;
    	  	//ASSERT (! al->refExactlyMatched ());
    	  	//ASSERT (! al->partial ());
      	  }
  	  	  al->qc ();
  	  	  batch. target2hmmAls [hmmAl->sseqid] << hmmAl. get ();  
  	      batch. hmmAls                        << hmmAl. release ();
    	  }
    	}
    }
   	if (verbose ())
   	  cout << "# Good initial HMMs: " << batch. hmmAls. size () << endl;
  
  
    if (   ! noblast
        && ! blastxFName. empty ()
       )       
    {
  		const Chronometer_OnePass cop ("blastx", cerr, false, Chronometer::enabled);  
      LineInput f (blastxFName);
  	  while (f. nextLine ())
  	  {
  	    { 
  	      Unverbose unv;
  	      if (verbose ())
  	        cout << f. line << endl;  
  	    }
  	    unique_ptr<BlastAlignment> al (new BlastAlignment (f. line, false));
  	    al->qc ();  
  	    if (nosame && al->refAccession == al->sseqid)
  	      continue;
 	      batch. blastAls << al. release ();
  	  }
  	}
  	if (verbose ())
  	  cout << "# Blasts: " << batch. blastAls. size () << endl;


    // For Batch::report()
    if (! gffFName. empty ())
    {
      const Chronometer_OnePass cop ("gff", cerr, false, Chronometer::enabled);  
    	unique_ptr<const Annot> annot;
    	if (bedP)
    	{
    	  QC_ASSERT (gffProtMatchFName. empty ());
    	  QC_ASSERT (gffDnaMatchFName. empty ());
    	  QC_ASSERT (! lcl);
		    annot. reset (new Annot (gffFName));
		  }
    	else
    	{
		    annot. reset (new Annot (gffFName, gffType, ! gffProtMatchFName. empty (), lcl));
  	    if (! gffProtMatchFName. empty ())
    		  var_cast (annot. get ()) -> load_fasta2gff_prot (gffProtMatchFName);	
  	    if (! gffDnaMatchFName. empty ())
    		  var_cast (annot. get ()) -> load_fasta2gff_dna (gffDnaMatchFName);	
		  }
		  ASSERT (annot);
	    {
	      // Locus::contigLen
  	    map<string,size_t> contig2len;
  	    if (! dnaLenFName. empty ())
  	    {
  	      LineInput f (dnaLenFName);
  	      string contig;
  	      size_t len;
  	      Istringstream iss;
  	      while (f. nextLine ())
  	      {
  	        iss. reset (f. line);
  	        len = 0;
  	        iss >> contig >> len;
  	        QC_ASSERT (len);
  	        contig2len [contig] = len;
  	      }
  	      for (const auto& it : annot->prot2loci)
  	        for (const Locus& locus : it. second)
  	          if (! locus. contigLen)
  	            var_cast (locus). contigLen = contig2len [locus.contig];
  	    }
  	  }
    	for (const BlastAlignment* al : batch. blastAls)
    	  if (al->sProt)
     	  	var_cast (al) -> setCdss (*annot);
	    for (const HmmAlignment* hmmAl : batch. hmmAls)
        var_cast (hmmAl->blastAl. get ()) -> setCdss (*annot);
    }
    
      
    if (! blastxFName. empty ())
    {
  		const Chronometer_OnePass cop ("target2...", cerr, false, Chronometer::enabled);  
      targetProt = false;
      // batch.target2blastAls
      // Proteins in BLASTP which are not in GFF are ignored
      batch. target2blastAls. clear ();  
    	for (const BlastAlignment* al : batch. blastAls)
        if (al->sProt)
    	    for (const string& contig : al->getContigs ())
    	      batch. target2blastAls [contig] << al;
    	  else
          batch. target2blastAls [al->sseqid] << al;
      // batch.target2hmmAls
      batch. target2hmmAls. clear ();
    	for (const HmmAlignment* hmmAl : batch. hmmAls)
  	  {
  	    ASSERT (hmmAl);
  	    ASSERT (hmmAl->blastAl);
   	    for (const string& contig : hmmAl->blastAl->getContigs ())
  	      batch. target2hmmAls [contig] << hmmAl;
  	  }
    }      

  
    batch. process (retainBlasts, skip_hmm_check);    


    // Output
    {
      TsvOut td (cout, 2, false);
      td. usePound = false;
      batch. report (td, false);
    }
    if (! mutation_all_FName. empty ())
    {
      OFStream f (mutation_all_FName);
      TsvOut td (f, 2, false);
      td. usePound = false;
      batch. report (td, true);
    }
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



