// fasta_check.cpp

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
*   Extract sequences out of a FASTA file
*
*/
   
   
#undef NDEBUG 

#include "common.hpp"
using namespace Common_sp;

#include "common.inc"



namespace 
{
  
  
  
struct Segment
// not circular
{
  size_t start {0};
  size_t stop {0};
  bool strand {true};
    // false <=> negative
  string genesymbol;
  string name;
  
  
  bool isDna () const
    { return stop; }
  size_t size () const
    { return stop - start; }
  void saveText (ostream &os) const
    { os         << start 
         << '\t' << stop 
         << '\t' << strand 
         << '\t' << genesymbol 
         << '\t' << name
         << endl;
    }
};



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
    case '-': r = '-'; break;
    default: 
    	throw runtime_error ("Bad nucleotide " + to_string (wildNucleotide));
  }
  if (isupper (wildNucleotide))
    r = toUpper (r);

  return r;
}



bool process (const string &id, 
              string &seq, 
              const map<string/*id*/,Vector<Segment>> &id2segments)
{
  if (id. empty ())
    return false;
  const Vector<Segment>* segments = findPtr (id2segments, id);
  if (! segments)
    return false;
    
  replaceStr (seq, "-", "");
  QC_ASSERT (! seq. empty ());
  
  for (Segment& seg : var_cast (*segments))
  {
    cout << '>' << id;
    if (seg. isDna ())
    {
      QC_ASSERT (seg. start <= seq. size ());
      minimize (seg. stop, seq. size ());
      QC_ASSERT (seg. start < seg. stop);
      cout << ':' << seg. start + 1 << '-' << seg. stop << ' ' << "strand:" << (seg. strand ? '+' : '-');
    }
    cout << ' ' << seg. genesymbol << ' ' << seg. name << endl;
    string seq1 (seq);
    if (seg. isDna ())
    {
      ASSERT (seg. stop <= seq1. size ());
      seq1 = seq1. substr (seg. start, seg. size ());
      if (! seg. strand)
      {
        reverse (seq1);
        for (char &c : seq1)
          c = complementaryNucleotide (c);
      }
    //strLower (seq1);  // Letter case can indicate nucleotide quality
    }    
  //else
    //strUpper (seq1);
    constexpr size_t line_len = 60;  // PAR
    for (size_t i = 0; i < seq1. size (); i += line_len)
      cout << seq1. substr (i, line_len) << endl;
  }
  
  return true;
}



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Extract sequences out of a FASTA file")
    {
      addPositional ("fasta", "FASTA file");
      addPositional ("target", "Target identifiers in the FASTA file to extract.\n\
Line format for amino acid sequences : <id> <gene symbol> <product name>\n\
Line format for nucleotide sequences : <id> <start (>=1)> <stop (>= start)> <strand (+/-)> <gene symbol> <product name>\
");
      addFlag ("aa", "Amino acid sequenes, otherwise nucleotide");
	    version = SVN_REV;
    }



  void body () const final
  {
    const string fName       = getArg ("fasta");
    const string targetFName = getArg ("target");
    const bool aa            = getFlag ("aa");
    
    
    map<string/*id*/,Vector<Segment>> id2segments;
    {
      LineInput f (targetFName);
      string id;
      Istringstream iss;
      while (f. nextLine ())
      {
        iss. reset (f. line);
        Segment seg;
        iss >> id;
        if (! aa)
        {
          char strand = '\0';
          iss >> seg. start >> seg. stop >> strand;
          QC_ASSERT (seg. start);
          QC_ASSERT (seg. start <= seg. stop);
          seg. start--;
          QC_ASSERT (   strand == '+' 
                     || strand == '-'
                    );
          seg. strand = (strand == '+');
        }
        iss >> seg. genesymbol;
        seg. name = f. line. substr ((size_t) iss. tellg ());
        trim (seg. name);
        QC_ASSERT (aa == ! seg. isDna ());
        id2segments [id] << std::move (seg);
      }
    }
    if (verbose ())
      for (const auto& it : id2segments)
      {
        cout << it. first << ": " << endl;
        for (const Segment& seg : it. second)
        {
          cout << "  ";
          seg. saveText (cout);
        }
      }
    if (id2segments. empty ())
      return;
    

    size_t processed = 0;
    {
      LineInput f (fName); 
      string id;
      string seq;
      while (f. nextLine ())
      {
        trimTrailing (f. line);
        if (f. line. empty ())
        	continue;
      	if (f. line [0] == '>')
      	{
      	  processed += process (id, seq, id2segments);
      		size_t pos = 1;
      		while (pos < f. line. size () && ! isspace (f. line [pos]))
      		  pos++;
      		id = f. line. substr (1, pos - 1);
      		seq. clear ();
      	}
      	else 
      	  seq += f. line;
  	  }
   	  processed += process (id, seq, id2segments);
   	}
   	if (processed != id2segments. size ())  
   	  throw runtime_error ("Requested identifiers: " + to_string (id2segments. size ()) + ", but processed: " + to_string (processed));
   	  // Assumed: no duplicate identifiers in FASTA
  }
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



