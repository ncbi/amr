// dna_point_mut.cpp

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
*   Identification of point mutations at DNA level
*
* Release changes:
*   3.2.4 11/15/2019 PD-3191  neighborhoodMismatch <= 0.04; good(): length >= min (refLen, 2 * flankingLen + 1)
*
*/
   
   
#undef NDEBUG 
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;
#include "point_mut.hpp"
using namespace PointMut_sp;




// PD-3096
// PAR!
#ifdef SVN_REV
  #define SOFTWARE_VER SVN_REV
#else
  #define SOFTWARE_VER "3.4.1"
#endif



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



#if 0
struct PointMut
{
  string gene;
	size_t pos {0};
	  // In reference sequence
	  // >= 0
	char alleleChar {' '};
	  // Upper-case
	// !empty()
	string geneMutation;
	  // Depends on the above, except for pos
	string classS;
	string subclass;
	string name;
	  // Species binomial + resistance
	double neighborhoodMismatch {0.0};

	
	PointMut (const string &gene_arg,
	          size_t pos_arg,
						const string &geneMutation_arg,
						const string &class_arg,
						const string &subclass_arg,
						const string &name_arg)
		: gene (gene_arg)
		, pos (pos_arg)
		, alleleChar (geneMutation_arg. back ())  // ??
		, geneMutation (geneMutation_arg)
		, classS (class_arg)
		, subclass (subclass_arg)
		, name (name_arg)
		{ 
			QC_ASSERT (! gene. empty ());
			QC_ASSERT (pos > 0);
			pos--;
			QC_ASSERT (alleleChar != ' ');
			alleleChar = toUpper (alleleChar);
			QC_ASSERT (! geneMutation. empty ());
			QC_ASSERT (! name. empty ());
      QC_ASSERT (! contains (name, '\t'));
      replace (name, '_', ' ');
      QC_ASSERT (geneMutation. back () == alleleChar);
		}
	PointMut () = default;


  bool empty () const
    { return gene. empty (); }
  bool better (const PointMut &other) const  
    { if (geneMutation != other. geneMutation)
        return false;
    /*if (pos != other. pos)
        return false; */
      LESS_PART (*this, other, neighborhoodMismatch);
      // Tie
      LESS_PART (*this, other, name); 
      LESS_PART (*this, other, geneMutation);  
      return false;
    }
  bool operator< (const PointMut &other) const
    { LESS_PART (*this, other, pos);
      LESS_PART (*this, other, geneMutation);
      return false;
    }
  bool operator== (const PointMut &other) const
    { return geneMutation == other. geneMutation; }
};
#endif



map <string/*accession*/, Vector<PointMut>>  accession2pointMuts;
unique_ptr<OFStream> point_mut_all;  // ??



struct BlastAlignment 
{
	// PD-2001
 	static constexpr const size_t flankingLen = 200;  // PAR

  // BLASTN alignment
  size_t length {0}, nident {0}  // aa
       ,    refStart {0},    refEnd {0},    refLen {0}
       , targetStart {0}, targetEnd {0}, targetLen {0};
    // Positions are 0-based
    // targetStart < targetEnd

  // target    
  string targetName; 
  string refAccessionFrag;
  string product;
  bool targetStrand {true}; 
    // false <=> negative
  
  // Reference
  string refName; 

//map<size_t/*targetPos*/,PointMut> targetPos2pointMut;
  Vector<PointMut> pointMuts;
  

  explicit BlastAlignment (const string &line)
    {
      try
      {
  	    static Istringstream iss;
  	    iss. reset (line);
        string targetSeq, refSeq;
  	    iss >> targetName >> refName >> length >> nident >> targetStart >> targetEnd >> targetLen >> refStart >> refEnd >> refLen >> targetSeq >> refSeq;
  	  // format:  qseqid      sseqid    length    nident         qstart         qend         qlen      sstart      send      slen         sseq
      // blastn:  ...         ...          733       733          62285        63017        88215         105       837       837          ...
  	    QC_ASSERT (! targetSeq. empty ());	
  	    QC_ASSERT (targetSeq. size () == refSeq. size ());    

		    normalizeAlignment (targetSeq, refSeq);

  	    QC_ASSERT (refStart != refEnd);
  	    bool refStrand = refStart < refEnd;  
  	    if (! refStrand)
  	      swap (refStart, refEnd);
  	    
  	    QC_ASSERT (targetStart < targetEnd);
  	      
  	    QC_ASSERT (refStart >= 1);
  	    QC_ASSERT (targetStart >= 1);
  	    QC_ASSERT (refStart < refEnd);
  	    QC_ASSERT (targetStart < targetEnd);
  	    refStart--;
  	    targetStart--;
  	    
  	    if (! refStrand)
  	    {
  	    	reverseDna (targetSeq);
  	    	reverseDna (refSeq);
  	    	toggle (refStrand);
  	    	toggle (targetStrand);
  	    }
  	    ASSERT (refStrand);
  	    
      	QC_ASSERT (targetEnd - targetStart <= targetSeq. size ());
      	QC_ASSERT (refEnd - refStart <= refSeq. size ());

        // refAccessionFrag, product
  	    {
  	      string s (refName);
  	      refAccessionFrag = findSplit (s, '@');
  	      product = findSplit (s, ':');
  	      replace (product, '_', ' ');
  	      QC_ASSERT (! s. empty ());
  	      refAccessionFrag += ":" + s;
  	    }

        // Cf. amr_report.cpp
  	    if (const vector<PointMut>* pointMuts_all = findPtr (accession2pointMuts, refName))
  	    {
  	    	if (verbose ())
  	        cout << "PointMut DNA found: " << refName << endl 
  	             << targetStart << ' ' << targetEnd << ' ' << targetStrand << ' ' << targetSeq << endl 
  	             << refStart    << ' ' << refEnd    << ' ' << refStrand    << ' ' << refSeq << endl;
  	      ASSERT (! pointMuts_all->empty ());
  	      
  	      
  	      // pmGene
  	      string s (pointMuts_all->front (). geneMutation);
  	      rfindSplit (s, pm_delimiter);  // Only gene symbol
  	      const string pmGene (s);

    		  Vector<SeqChange> seqChanges; 
    		  {
  	  		  SeqChange seqChange;
  	  		  bool inSeqChange = false;
  	    		FFOR (size_t, i, refSeq. size ())
  	    		  if (inSeqChange)
  	    		  {
  			        if (targetSeq [i] == refSeq [i])
  			        {
  			          seqChange. setSeq (targetSeq, refSeq);
  			          if (refSeq [seqChange. start] == '-')
  			          {
  			            QC_ASSERT (seqChange. start);
  			            seqChange. start--;
  			            seqChange. len++;
    			          seqChange. setSeq (targetSeq, refSeq);
  			          }
   			          ASSERT (! seqChange. reference. empty ());
   			          seqChange. setRefStart (refSeq, refStart);
   			          seqChange. setTargetStart (targetSeq, targetStart, targetStrand);
    		    	    seqChange. neighborhoodMismatch = getNeighborhoodMismatch (targetSeq, refSeq, seqChange);
    		    	    if (verbose ())
    		    	      seqChange. saveText (cout);
  			    			if (seqChange. neighborhoodMismatch <= 0.04)   // PAR  // PD-3191
   			            seqChanges << move (seqChange);
  			          seqChange = SeqChange ();
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
    		  while (j < pointMuts_all->size () && pointMuts_all->at (j). pos < refStart)
    		  {
    		  	if (point_mut_all. get ())
    		  	{
    		  		PointMut pm (pointMuts_all->at (j));
    		  		pm. name = "Non-called " + pmGene;
    		  		pm. additional = true;
    		  		pointMuts << move (pm);
    		  	}
    		    j++;
    		  }
    		  
      		size_t refStart_prev = NO_INDEX;
  			  for (const SeqChange& seqChange : seqChanges)
  			  {
  			    if (verbose ())
  			      seqChange. saveText (cout);
  			    QC_IMPLY (refStart_prev != NO_INDEX, refStart_prev < seqChange. refStart);
  			    while (j < pointMuts_all->size ())
  			    {
      			  const PointMut& pm = pointMuts_all->at (j);
      			  if (verbose ())
      			  {
      			    cout << "Processing PointMut: ";
      			    pm. print (cout);
      			  }
      			  if (pm. pos < seqChange. refStart)
  		      	{
  				    	if (point_mut_all. get ())
  				    		pointMuts << move (PointMut (pmGene, pm. pos, pm. reference, pm. reference));
  		      	}      			    
      			  if (pm. pos > seqChange. refStart)
      			    break;
      			  if (   pm. pos       == seqChange. refStart
      			      && pm. reference == seqChange. reference
      			      && pm. allele    == seqChange. allele
      			     )
      			  {
  	    		  	pointMuts << pm;
  	    		  	PointMut& pm_ = pointMuts. back ();
  	    		  	pm_. neighborhoodMismatch = seqChange. neighborhoodMismatch;
  	    		  	pm_. targetStart          = seqChange. targetStart;
  	    		  	if (verbose ())
  	    		  	{
  	    		  	  cout << "Found: ";
  	    		  	  pm_. print (cout);
  	    		  	}
  	    		  }
      		    j++;
  	    		}
    			  refStart_prev = seqChange. refStart;
      		}
  			    
   		  	if (point_mut_all. get ())
  	  		  while (j < pointMuts_all->size ())
  	  		  {
  			  		PointMut pm (pointMuts_all->at (j));
  			  		pm. name = "Non-called " + pmGene;
  			  		pm. additional = true;
  			  		pointMuts << move (pm);
  	  		    j++;
  	  		  }
  	  		  
  	  		pointMuts. sort ();
  	  		pointMuts. uniq ();

  	      
  	    #if 0
  	    	for (const PointMut& pm : *pointMuts_)
  	    	{
  	    		size_t pos = refStart;
  	    		// One pass for all *pointMuts_ ??
  	    		size_t targetPos = (targetStrand ? targetStart : (targetEnd - 1));
  	    		FFOR (size_t, i, refSeq. size ())
  	    		{
  	    		  if (targetSeq [i] != '-')
  	    		  {
  	    		    if (targetStrand)
  	    		      targetPos++;
  	    		    else
  	    		      targetPos--;
  	    		  }
  	    		  if (refSeq [i] != '-')
  		    	  {
  		    	    if (targetSeq [i] != '-')
  		    	    {
    		    	    const double neighborhoodMismatch = getNeighborhoodMismatch (targetSeq, refSeq, i);
    		    	  	if (verbose ())
    			    	  	if (targetSeq [i] != refSeq [i])  
    			    	  		cout        << i + 1 
    			    	  		     << ' ' << refSeq [i]  
    			    	  		     << ' ' << targetSeq [i] 
    			    	  		     << ' ' << pos + 1 
    			    	  		     << ' ' << pm. pos + 1
    			    	  		     << ' ' << pm. allele
    			    	  		     << ' ' << neighborhoodMismatch
    			    	  		     << endl;
    		    	  	if (pos == pm. pos)
    		    	  	{
    				    		if (toUpper (targetSeq [i]) == pm. allele [0])  // ??
    				    		{
    				    			ASSERT (targetSeq [i] != refSeq [i]);
    				    			if (neighborhoodMismatch <= 0.04)   // PAR  // PD-3191
    				    			{
    				    			  ASSERT (targetPos2pointMut [targetPos]. empty ());
    				    			  targetPos2pointMut [targetPos] = pm;
    				    			  targetPos2pointMut [targetPos]. neighborhoodMismatch = neighborhoodMismatch;
    				    			}
    				    		}
    		    	  		break;
    		    	  	}
    		    	  }
  	    		  	pos++;
  		    		}
  		    	}
  	      }
  	    #endif
  	    }
		  }
		  catch (...)
		  {
		  	cout << line << endl;
		  	throw;
		  }
    }
  void qc () const
    {
      if (! qc_on)
        return;
	    QC_ASSERT (targetStart < targetEnd);
	    QC_ASSERT (targetEnd <= targetLen);
      QC_ASSERT (refStart < refEnd);
	    QC_ASSERT (nident <= refEnd - refStart);
	    QC_ASSERT (refEnd <= refLen);
	    QC_ASSERT (refEnd - refStart <= length);	    
    }
  void saveText (ostream& os) const 
    { const string na ("NA");
    //for (const auto& it : targetPos2pointMut)
      for (const PointMut& pm : pointMuts)
      {
      //const PointMut& pm = it. second;
        if (pm. empty ())
          continue;
        TabDel td (2, false);
        td << na  // PD-2534
           << targetName 
           << targetStart + 1
           << targetEnd
           << (targetStrand ? '+' : '-');
        td << pm. geneMutation
           << pm. name
           << "core"  // PD-2825
           // PD-1856
           << "AMR"
           << "POINT"
           << nvl (pm. classS, na)
           << nvl (pm. subclass, na)
           //
           << "POINTN"  // PD-2088
           << targetLen;
        td << refLen
           << refCoverage () * 100  
           << pIdentity () * 100  
           << length
           << refAccessionFrag  // refName
           << product  // pm. gene
           ;
        // HMM
        td << na
           << na;
        os << td. str () << endl;
      }
    }
    

  double pIdentity () const
    { return (double) nident / (double) length; }
  double refCoverage () const
    { return (double) (refEnd - refStart) / (double) refLen; }
  bool good () const
    { return length >= min (refLen, 2 * flankingLen + 1); }
  bool operator< (const BlastAlignment &other) const
    { LESS_PART (*this, other, targetName);
      LESS_PART (other, *this, pIdentity ());
      LESS_PART (*this, other, targetStart);
      LESS_PART (*this, other, refName);
      return false;
    }
  double getNeighborhoodMismatch (const string &targetSeq, 
                                  const string &refSeq, 
                                  const SeqChange &seqChange) const
    { ASSERT (targetSeq. size () == refSeq. size ());
      ASSERT (seqChange. start < targetSeq. size ())
    //ASSERT (targetSeq [i] != '-');
    //ASSERT (refSeq    [i] != '-');
      ASSERT (targetEnd - targetStart <= targetSeq. size ());
      ASSERT (refEnd - refStart <= refSeq. size ());
      // PD-2001
      size_t len = 0;
      size_t mismatches = 0; 
      size_t j = 0;
      // Left flank
      if (seqChange. start)
      {
        j = seqChange. start - 1; 
        while (seqChange. start - j <= flankingLen)
        { len++;
          if (targetSeq [j] != refSeq [j])
                mismatches++;
          if (j == 0)
                break;
          j--;
        }
      }
      ASSERT (seqChange. start >= j);
      const size_t left = min (min (targetStart, refStart), flankingLen + 1 - (seqChange. start - j)); 
      len        += left;
      mismatches += left; 
      if (verbose ())
        cout << "left: " << seqChange. start  << ' ' << j << ' ' << left << endl;
      // Right flank
      const size_t stop = seqChange. getStop ();
      for (j = stop + 1; j < targetSeq. size () && j < refSeq. size () && j - stop <= flankingLen; j++)
      { len++;
        if (targetSeq [j] != refSeq [j])
          mismatches++;
      }
      ASSERT (j >= stop);
      const size_t right = min (min (targetLen - targetEnd, refLen - refEnd), flankingLen + 1 - (j - stop));
      len        += right;
      mismatches += right;     
      if (verbose ())
        cout << "right: " << j << ' ' << right << ' ' << mismatches << ' ' << len << endl;  
      //
      return (double) mismatches / (double) len;
    }
};




struct Batch
{
  vector<BlastAlignment> blastAls;   
  
  
  explicit Batch (const string &point_mut)
	  {
	    {
        LineInput f (point_mut);
  	  	string accession, geneMutation, classS, subclass, name;
  			int pos;
   	  	Istringstream iss;
    	  while (f. nextLine ())
    	  {
     	  	iss. reset (f. line);
    	  	iss >> accession >> pos >> geneMutation >> classS >> subclass >> name;
    	  	QC_ASSERT (pos > 0);
   	  		accession2pointMuts [accession]. push_back (move (PointMut ((size_t) pos, geneMutation, classS, subclass, name)));
    	  }	    
    	}
  	  for (auto& it : accession2pointMuts)
  	  {
  	  	it. second. sort ();
  	    if (! it. second. isUniq ())
  	  	  throw runtime_error ("Duplicate mutations for " + it. first);
  	  }
	  }
	  	  

	void report (ostream &os) const
	{
    {
    	// Cf. BlastAlignment::saveText()
	    TabDel td;
	    td << "Protein identifier"   // targetName  // PD-2534
         // Contig
         << "Contig id"
         << "Start"  // targetStart
         << "Stop"  // targetEnd
         << "Strand"   // targetStrand
	       << "Gene symbol"
	       << "Mutation name"
	       << "Scope"  // PD-2825
	       // PD-1856
	       << "Element type"
	       << "Element subtype"
	       << "class"
	       << "Subclass"
	       //
	       << "Method"
	       << "Target length" 
	       //
	       << "Reference gene length"         // refLen
	       << "% Coverage of reference gene"  // queryCoverage
	       << "% Identity to reference gene"  
	       << "Alignment length"                 // length
	       << "Accession of reference gene"    
	       << "Name of reference gene"
	       //
	       << "HMM id"
	       << "HMM description"
	       ;
	    os << td. str () << endl;
	  }

  	for (const auto& blastAl : blastAls)
  	{
   	  blastAl. qc ();
   	  blastAl. saveText (os);
    }
	}
};




// ThisApplication

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Find point mutations at DNA level and report in the format of amr_report.cpp")
    {
      addPositional ("blastn", "blastn output in the format: qseqid sseqid length nident qstart qend qlen sstart send slen sseq. sseqid is the 1st column of <point_mut> table");  
      addPositional ("point_mut", "Point mutation table");
	    version = SOFTWARE_VER;
    }



  void body () const final
  {
    const string blastnFName = getArg ("blastn");
    const string point_mut   = getArg ("point_mut");  
    
    

    Batch batch (point_mut);  
  
  
    // Input 
    {
      LineInput f (blastnFName);
  	  while (f. nextLine ())
  	  {
  	    { 
  	      Unverbose unv;
  	      if (verbose ())
  	        cout << f. line << endl;  
  	    }
  	    BlastAlignment al (f. line);
  	    al. qc ();  
  	    if (al. good ())
  	      batch. blastAls. push_back (move (al));
  	  }
  	}
  	if (verbose ())
  	  cout << "# Good Blasts: " << batch. blastAls. size () << endl;
  	
    
    // Output
    // Group by targetName and process each targetName separately for speed ??    
  //Common_sp::sort (batch. blastAls);
    if (verbose ())
    {
	    cout << "After process():" << endl;
		  for (const auto& blastAl : batch. blastAls)
		  {
		    blastAl. saveText (cout);
		    cout << ' ' << blastAl. pointMuts /*targetPos2pointMut*/. size () << endl;
		  }
		}
		
    
    for (const BlastAlignment& blastAl1 : batch. blastAls)
    //for (const auto& it1 : blastAl1. targetPos2pointMut)
      for (const PointMut& pm1 : blastAl1. pointMuts)
      {
      //const PointMut& pm1 = it1. second;
        if (pm1. empty ())
          continue;
        for (BlastAlignment& blastAl2 : batch. blastAls)
          if (   blastAl2. targetName == blastAl1. targetName
              && & blastAl2 != & blastAl1
             )  
          //for (auto& it2 : blastAl2. targetPos2pointMut)
            for (PointMut& pm2 : blastAl2. pointMuts)
            {
            //PointMut& pm2 = it2. second;
              if (pm2. empty ())
                continue;
              if (verbose ())
                cout        << blastAl2. targetName 
                   //<< ' ' << it1. first 
                   //<< ' ' << it2. first 
                     << ' ' << pm1. geneMutation 
                     << ' ' << pm2. geneMutation 
                     << ' ' << pm1. neighborhoodMismatch 
                     << ' ' << pm2. neighborhoodMismatch 
                     << ' ' << pm1. better (pm2)
                     << endl; 
              if (   pm1. targetStart == pm2. targetStart  // it2. first == it1. first
                  && pm1. better (pm2)
                 )
                pm2 = PointMut ();
            }
      }
		
		
    batch. report (cout);
  }
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



