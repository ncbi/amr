// disruption2genesymbol.cpp

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
*   Convert Disruption::genesymbol_raw() to a gene symbol
*
*/


#undef NDEBUG

#include "common.hpp"
using namespace Common_sp;
#include "seq.hpp"
using namespace Seq_sp;

#include "common.inc"



namespace 
{
  
  
constexpr char no_aa {'?'};

  
  
struct SymbolRaw final : Root
{
  static constexpr size_t del_size {10};  // PAR

  // Input
  string contig;
  string prot;
  Disruption::Type type {Disruption::eNone};
  // 0-based  
  size_t qstart {no_index}; 
  size_t qend {no_index};
    // aa
  size_t sstart {no_index};
  size_t send {no_index};
    // bp
  //
  Strand strand {0};
  bool stop {false};
  string rest;
  
  // Output
  string ref;
  string allele;


  explicit SymbolRaw (const string &line)
    {
      string s;
      {
        istringstream iss (line);
        iss >> contig >> prot >> s;
        QC_ASSERT (! s. empty ());
        ASSERT (! contig. empty ());
        ASSERT (! prot. empty ());

        constexpr size_t rest_size = 1024;
        char rest_ [rest_size];
        iss. getline (rest_, rest_size);
        rest = rest_;
      }

      if (trimSuffix (s, Disruption::stopSuf))
        stop = true;

      // strand
      {
        const string strandS = rfindSplit (s, '_');
        if (strandS == "0")
          strand = -1;
        else if (strandS == "1")
          strand = 1;
        else
          throw runtime_error ("Unknown strand: " + strQuote (strandS));
      }
          
      send   = str2<size_t> (rfindSplit (s, '_'));
      sstart = str2<size_t> (rfindSplit (s, '_'));      
      qend   = str2<size_t> (rfindSplit (s, '_'));
      qstart = str2<size_t> (rfindSplit (s, '_'));       
      QC_ASSERT (qstart <= qend);
      QC_ASSERT (sstart <= send);
      
      type = Disruption::name2type (s);
      QC_ASSERT (type != Disruption::eNone);
      QC_ASSERT (type != Disruption::eSmooth);
    }  
  void saveText (ostream &os) const final
    { 
      ASSERT (! ref. empty ());
      os         << contig  // 0
         << '\t' << prot    // 1 
         << '\t';
      if (verbose ())
        os << '\t' << Disruption::typeNames [type]
           << '\t' << qstart
           << '\t' << qend
           << '\t' << sstart
           << '\t' << send
           << '\t' << (int) strand
           << '\t' << stop
           << '\t' << ref
           << '\t' << allele
           << '\t';           
      ASSERT (! contains (ref, '*'));
      string allele_ (allele);
      const bool alleleStop = trimSuffix (allele_, "*");
      const size_t allele_size = allele_. size ();  // Without stop codon
    //QC_IMPLY (type != Disruption::eFrameshift, alleleStop == stop);
      QC_IMPLY (stop, alleleStop);
      constexpr size_t display_max = 1/*reference aa*/ + 5;  // PAR  // PD-5395
      if (allele_size > display_max)
        allele_ = "ins";
      if (alleleStop)
        allele_ += terminatorWord;
      ASSERT (! contains (allele_, '*'));
      // Standard gene symbol        
      // 2      
      if (ref. size () > display_max)
        os        << ref. front () << qstart + 1 
           << '_' << ref. back ()  << qstart + ref. size ();
      else
        os << ref << qstart + 1;  
      switch (type)
      {
        case Disruption::eFrameshift:
          ASSERT (ref. size () == 1)
          ASSERT (! allele. empty ());
          if (alleleStop && allele_size == 0)            
            os << terminatorWord;
          else
            os << allele [0];
          os << Disruption::typeNames [type];
          if (alleleStop)
            os << terminatorWord << allele_size;
          break;
        case Disruption::eDeletion:  // Or replacement
          if (allele_. empty ())
            os << Disruption::typeNames [type];
          else
          {
            os << allele_;
            if (allele_size > display_max)
              os << allele_size - 1/*reference aa*/;
          }
          break;
        case Disruption::eInsertion:
          ASSERT (ref. size () == 1);
          ASSERT (! allele_. empty ());
          os << allele_;
          if (allele_size > display_max)
            os << allele_size - 1/*reference aa*/;
          break;
        default:
          break;
      }
      // 3
      os << '\t' 
         // = <Disruption::genesymbol_raw()>
         // Opposite to SymbolRaw::SymbolRaw(line)
         << Disruption::typeNames [type] << '_' << qstart << '_' << qend << '_' << sstart << '_' << send << '_' << (strand == 1 ? 1 : 0);
      if (stop)
        os << Disruption::stopSuf;
         //
      os << '\t' << rest  // 4
         << '\n';
    }
    
    
  char contig2aa (const Dna &dna,
                  size_t offset,
                  Gencode gencode) const
    // Input: offset: from sstart/send
    // Return: no_aa <=> offset is outside dna
    {
      QC_ASSERT (send <= dna. seq. size ());
      
      if (strand == 1)
      {
        const size_t i = sstart + offset * 3;
        if (i + 3 > send)
          return no_aa;
        return codon2aa (& dna. seq [i], gencode, false);
      }
      
      ASSERT (strand == -1);
      if (send < (offset + 1) * 3)
        return no_aa;
      const size_t i = send - (offset + 1) * 3;
      ASSERT (i + 3 <= dna. seq. size ());
      if (i < sstart)
        return no_aa;
      string s (dna. seq. substr (i, 3));
      reverseDna (s);
      return codon2aa (s. c_str (), gencode, false);
    }
};

	
	
struct ThisApplication final : Application
{
  static constexpr char id_delim {'|'};
  
  
  ThisApplication ()
    : Application ("Convert Disruption::genesymbol_raw() to standard gene symbols according to https://hgvs-nomenclature.org/stable/recommendations/protein/frameshift/.\n\
A stop codon is '" + string (terminatorWord) + "'.\n\
Print: <tab row> where <genesymbol> is inserted before <Disruption::genesymbol_raw()>"
)
	  {
		  addPositional ("nucl", "Input nucleotide FASTA file");
		  addPositional ("prot", "Input protein FASTA file");
		  addPositional ("tab", "Table with lines: <contig identifier in <nucl>>  <protein identifier in <prot>>  <Disruption::genesymbol_raw()>");
		  addKey ("gencode", "NCBI genetic code for translated BLAST", "11");
	    addKey ("prot_id_pos", string ("Position of protein id in qseqid delimited by ") + id_delim + ", 1-based. 0 - use qseqid as a whole", "0");
	  }

  
  
	void body () const final
  {
	  const string nuclFName   = getArg ("nucl");
	  const string protFName   = getArg ("prot");
	  const string tabFName    = getArg ("tab");
    const Gencode gencode    = (Gencode) arg2uint ("gencode"); 
    const size_t prot_id_pos = str2<size_t> (getArg ("prot_id_pos"));
    
	  
	  Vector<SymbolRaw> symbolRaws;
	  {
	    LineInput f (tabFName);
	    while (f. nextLine ())
	      symbolRaws << std::move (SymbolRaw (f. line));
	  }
	  if (symbolRaws. empty ())
	    return;
	  
	  // SymbolRaw::allele
		{
		  Multifasta fa (nuclFName, false);
		  while (fa. next ())
		  {
	      const Dna dna (fa, 100000/*PAR*/, true);
		    dna. qc ();	
		    const string id (dna. getId ());
		    for (SymbolRaw& symbolRaw : symbolRaws)
		      if (symbolRaw. contig == id)
  		      for (size_t offset = 0; ; offset++)
  		      {
  		        const char aa = symbolRaw. contig2aa (dna, offset, gencode);
  		        if (aa == no_aa)
  		          break;
  		        symbolRaw. allele += aa;
  		        if (aa == '*')
  		          break;
  		      }
		  }
		}

    // SymbolRaw::{ref, allele for "del"}
		{
		  Multifasta fa (protFName, true);
		  while (fa. next ())
		  {
		    const Peptide pep (fa, 1000/*PAR*/, true);
		    pep. qc ();	
		    
		    string id;
		    const string id_whole (pep. getId ());
		    if (prot_id_pos)
		    {
    		  const StringVector vec (id_whole, id_delim, true);
    		  if (prot_id_pos - 1 >= vec. size ())
    		    throw runtime_error ("Protein identifier position " + to_string (prot_id_pos) + " is outside of the list of identifiers: " + strQuote (id_whole));
    		  id = vec [prot_id_pos - 1];
    		}
		    else
		      id = id_whole;
    		
		    for (SymbolRaw& symbolRaw : symbolRaws)
		      if (symbolRaw. prot == id)
		        symbolRaw. ref = pep. seq. substr (symbolRaw. qstart, symbolRaw. qend - symbolRaw. qstart);
		  }
		}
		
		// symbolRaw's
		for (const SymbolRaw& symbolRaw : symbolRaws)
		  symbolRaw. saveText (cout);
  }
};


}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



