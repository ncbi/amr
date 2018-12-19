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
*   Identification of point mutations at DNA level
*
*/
   
   
#undef NDEBUG 
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;



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
	  // Depends on the above
	string geneMutationGen;
	  // geneMutation generalized, may have a different pos
	string name;
	  // Species binomial + resistance

	
	PointMut (const string &gene_arg,
	          size_t pos_arg,
						char alleleChar_arg,
						const string &geneMutation_arg,
						const string &geneMutationGen_arg,
						const string &name_arg)
		: gene (gene_arg)
		, pos (pos_arg)
		, alleleChar (alleleChar_arg)
		, geneMutation (geneMutation_arg)
		, geneMutationGen (geneMutationGen_arg)
		, name (name_arg)
		{ 
			ASSERT (! gene. empty ());
			ASSERT (pos > 0);
			pos--;
			ASSERT (alleleChar != ' ');
			alleleChar = toUpper (alleleChar);
			ASSERT (! geneMutation. empty ());
			ASSERT (! geneMutationGen. empty ());
			ASSERT (! name. empty ());
      ASSERT (! contains (name, '\t'));
      replace (name, '_', ' ');
      ASSERT (geneMutation. back () == alleleChar);
      ASSERT (geneMutationGen. back () == alleleChar);
      ASSERT (geneMutation. front () != alleleChar);
      ASSERT (geneMutationGen. front () != alleleChar);
		}
	PointMut ()
	  {}


  bool empty () const
    { return gene. empty (); }
  string getResistance () const
    { size_t p = name. find (' ');
    	ASSERT (p != string::npos);
    	p = name. find (' ', p + 1);
    	ASSERT (p != string::npos);
    	return name. substr (p + 1);
    }
};


map <string/*accession*/, vector<PointMut>>  accession2pointMuts;



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
  bool targetStrand {true}; 
    // false <=> negative
  
  // Reference
  string refName; 

  map<size_t/*targetPos*/,PointMut> targetPos2pointMut;
  

  explicit BlastAlignment (const string &line)
    {
	    static Istringstream iss;
	    iss. reset (line);
      string targetSeq, refSeq;
	    iss >> targetName >> refName >> length >> nident >> targetStart >> targetEnd >> targetLen >> refStart >> refEnd >> refLen >> targetSeq >> refSeq;
	  // format:  qseqid      sseqid    length    nident         qstart         qend         qlen      sstart      send      slen         sseq
    // blastn:  ...         ...          733       733          62285        63017        88215         105       837       837          ...
	    ASSERT (! targetSeq. empty ());	
	    ASSERT (targetSeq. size () == refSeq. size ());    

    #if 0	    
	    string refAccession;
	    size_t refSegStart = 0;
	    size_t refSegEnd = 0;
	    {
		    string s (refName);
		    replace (s, ':', ' ');
		    replace (s, '-', ' ');
		    static Istringstream refNameIss;
		    refNameIss. reset (s);
		    refNameIss >> refAccession >> refSegStart >> refSegEnd;
		    ASSERT (refSegStart);
		    ASSERT (refSegStart < refSegEnd);
		    refSegStart--;
		  }
		#endif

	    ASSERT (refStart != refEnd);
	    bool refStrand = refStart < refEnd;  
	    if (! refStrand)
	      swap (refStart, refEnd);
	    
	    ASSERT (targetStart < targetEnd);
	      
	    ASSERT (refStart >= 1);
	    ASSERT (targetStart >= 1);
	    ASSERT (refStart < refEnd);
	    ASSERT (targetStart < targetEnd);
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

	    if (const vector<PointMut>* pointMuts_ = findPtr (accession2pointMuts, refName))
	    {
	    	if (verbose ())
	        cout << "PointMut DNA found: " << refName << endl;
	      ASSERT (! pointMuts_->empty ());
	    	for (const PointMut& pm : *pointMuts_)
	    	{
	    		size_t pos = refStart;
	    		// One pass for all *pointMuts_ ??
	    		FFOR (size_t, i, refSeq. size ())
	    		  if (refSeq [i] != '-')
		    	  {
		    	  	if (verbose ())
			    	  	if (targetSeq [i] != refSeq [i])  
			    	  		cout << i + 1 << ' ' << refSeq [i]  << ' ' << targetSeq [i] << ' ' << pos + 1 << endl;
		    	  	if (pos == pm. pos)
		    	  	{
				    		if (toUpper (targetSeq [i]) == pm. alleleChar)
				    		{
				    			ASSERT (targetSeq [i] != refSeq [i]);
				    			if (goodNeighborhood (targetSeq, refSeq, i))
				    			{
				    			  const size_t targetPos = (targetStrand ? targetStart + i : (targetEnd - 1 - i));
				    			  ASSERT (targetPos2pointMut [targetPos]. empty ());
				    			  targetPos2pointMut [targetPos] = pm;
				    			}
				    		}
		    	  		break;
		    	  	}
	    		  	pos++;
		    		}
	      }
	    }
    }
  void qc () const
    {
      if (! qc_on)
        return;
	    ASSERT (targetStart < targetEnd);
	    ASSERT (targetEnd <= targetLen);
      ASSERT (refStart < refEnd);
	    ASSERT (nident <= refEnd - refStart);
	    ASSERT (refEnd <= refLen);
	    ASSERT (refEnd - refStart <= length);	    
    }
  void saveText (ostream& os) const 
    { const string na ("NA");
      for (const auto& it : targetPos2pointMut)
      {
        const PointMut& pm = it. second;
        if (pm. empty ())
          continue;
        TabDel td (2, false);
        td << targetName;
        td << targetName 
           << targetStart + 1
           << targetEnd
           << (targetStrand ? '+' : '-');
        td << pm. geneMutation
           << pm. name
           << "POINTN"  // PD-2088
           << targetLen;
        td << refLen
           << refCoverage () * 100  
           << pIdentity () * 100  
           << length
           << refName
           << pm. gene
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
  bool goodNeighborhood (const string &targetSeq, 
                         const string &refSeq, 
                         size_t i) const
    { ASSERT (targetSeq. size () == refSeq. size ());
    	ASSERT (i < targetSeq. size ())
    	if (i == 0)
    		return false;
    	// PD-2001
      size_t mismatches = 0;
      size_t j = i - 1; 
      while (i - j <= flankingLen)
      { if (targetSeq [j] != refSeq [j])
        	mismatches++;
        if (j == 0)
        	break;
        j--;
      }
      if (i - j != flankingLen + 1)
      	return false;
      for (j = i + 1; j < targetSeq. size () && j < refSeq. size () && j - i <= flankingLen; j++)
        if (targetSeq [j] != refSeq [j])
        	mismatches++;
      if (j - i != flankingLen + 1)
      	return false;
      return (double) mismatches / (2.0 * flankingLen) <= 0.03;  // PAR
    }
  bool good () const
    { return length >= 2 * flankingLen + 1; }
  bool operator< (const BlastAlignment &other) const
    { LESS_PART (*this, other, targetName);
      LESS_PART (*this, other, targetStart);
      LESS_PART (*this, other, refName);
      return false;
    }
};




struct Batch
{
  vector<BlastAlignment> blastAls;   
  
  
  explicit Batch (const string &point_mut)
	  {
      LineInput f (point_mut);
 	  	Istringstream iss;
  	  while (f. nextLine ())
  	  {
  	  	string accession;
  	  	string gene;
				int pos;
				char alleleChar;
				string geneMutation;
				string geneMutationGen;
				string name;
   	  	iss. reset (f. line);
  	  	iss >> accession >> gene >> pos >> alleleChar >> geneMutation >> geneMutationGen >> name;
  	  	ASSERT (pos > 0);
 	  		accession2pointMuts [accession]. push_back (PointMut (gene, (size_t) pos, alleleChar, geneMutation, geneMutationGen, name));
  	  }	    
	  }
	  	  

	void report (ostream &os) const
	{
    {
    	// Cf. BlastAlignment::saveText()
	    TabDel td;
	    td << "Target identifier"   // targetName
         // Contig
         << "Contig id"
         << "Start"  // targetStart
         << "Stop"  // targetEnd
         << "Strand"   // targetStrand
	       << "Gene symbol"
	       << "Mutation name"
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
      // Input
      addPositional ("blastn", "blastn output in the format: qseqid sseqid length nident qstart qend qlen sstart send slen sseq. sseqid is the 1st column of <point_mut> table");  
      addPositional ("point_mut", "Point mutation table");
      // Output
    #ifdef SVN_REV
      version = SVN_REV;
    #endif
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
  	    const BlastAlignment al (f. line);
  	    al. qc ();  
  	    if (al. good ())
  	      batch. blastAls. push_back (al);
  	  }
  	}
  	if (verbose ())
  	  cout << "# Good Blasts: " << batch. blastAls. size () << endl;
  	
    
    // Output
    // Group by targetName and process each targetName separately for speed ??    
    Common_sp::sort (batch. blastAls);
    if (verbose ())
    {
	    cout << "After process():" << endl;
		  for (const auto& blastAl : batch. blastAls)
		  {
		    blastAl. saveText (cout);
		    cout << ' ' << blastAl. targetPos2pointMut. size () << endl;
		  }
		}
		
    
    FFOR (size_t, i, batch. blastAls. size ())
    {
      const BlastAlignment& blastAl1 = batch. blastAls [i];
      for (const auto& it1 : blastAl1. targetPos2pointMut)
      {
        const PointMut& pm1 = it1. second;
        FFOR_START (size_t, j, i + 1, batch. blastAls. size ())
        {
          BlastAlignment& blastAl2 = batch. blastAls [j];
          if (blastAl2. targetName == blastAl1. targetName)  
            for (auto& it2 : blastAl2. targetPos2pointMut)
            {
              PointMut& pm2 = it2. second;
              if (verbose ())
                cout        << blastAl2. targetName 
                     << ' ' << it1. first 
                     << ' ' << it2. first 
                     << ' ' << pm1. geneMutationGen 
                     << ' ' << pm2. geneMutationGen 
                     << endl; 
              if (   it1. first == it2. first
                  && pm1. geneMutationGen == pm2. geneMutationGen
                 )
                pm2 = PointMut ();
            }
        }
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



