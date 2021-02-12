// amrfinder_customize.cpp

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
*   Customize AMRFinder database
*
* Dependencies: cp, makeblastdb
*
* Release changes: 
*
*/


#ifdef _MSC_VER
  #error "UNIX is required"
#endif
   
#undef NDEBUG 
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;



namespace 
{



struct ThisApplication : ShellApplication
{
  ThisApplication ()
    : ShellApplication ("Create a custom database for AMRFinder", true, true, false)
    {
    	addKey ("database_in", "Directory of the original AMRFinder database", "$BASE/data", 'd', "DATABASE_IN_DIR");
      addKey ("database_out", "Customized AMRFinder database directory for output", "", 'o', "DATABASE_OUT_DIR");
      addKey ("prot", "\
Protein FASTA file with header lines: >protein_identifier gene_symbol\n\
    where gene_symbol is a new family id or an existing family id in the column \"fam_id\" of DATABASE_IN_DIR/fam.tab", "", 'p', "PROT");
      addKey ("metadata", "\
Protein sequence metadata tab-delimited file for new non-mutation gene symbols with the header:\n\
        #gene_symbol\treportable\ttype\tsubtype\tclass\tsubclass\tprotein_name\n"
    //   0            1           2     3        4      5         6
"    where reportable: 1 - plus, 2 - core", "", 'm', "METADATA");
      addKey ("mutation", "\
Protein mutations tab-delimited file with the header:\n\
        #protein_identifier\tgene_symbol\torganism\tmutation\tclass\tsubclass\tmutated_protein_name\n"
    //   0                   1            2         3         4      5         6
"    where organism is either in the column \"taxgroup\" of the file DATABASE_IN_DIR/taxgroup.tab or is a new name\n\
    mutation = {<reference sequence>|INS}<position>{<allele sequence>|DEL}", "", 't', "MUTATION");
      addFlag ("quiet", "Suppress messages to STDERR", 'q');
	    version = SVN_REV;
    }



  void shellBody () const final
  {
    const string dbDirIn   = getArg ("database_in");
    const string dbDirOut  = getArg ("database_out");
    const string protFName = getArg ("prot");
    const string metaFName = getArg ("metadata");
    const string pmFName   = getArg ("mutation");
    const bool   quiet     = getFlag ("quiet");
    
        
    Stderr stderr (quiet);
    stderr << "Running: "<< getCommandLine () << '\n';
    const Verbose vrb (qc_on);
    const string qcS (qc_on ? " -qc" : "");

    if (dbDirIn. empty ())
      throw runtime_error ("Parameter -database_in is not specified");
    if (dbDirOut. empty ())
      throw runtime_error ("Parameter -database_out is not specified");
    if (protFName. empty ())
      throw runtime_error ("Parameter -prot is not specified");
    if (metaFName. empty ())
      throw runtime_error ("Parameter -metadata is not specified");
    if (pmFName. empty ())
      throw runtime_error ("Parameter -mutation is not specified");
    
    if (! directoryExists (dbDirIn))
      throw runtime_error ("Input database directory " + dbDirIn + " does not exist");
    if (directoryExists (dbDirOut))
      throw runtime_error ("Output database directory " + dbDirOut + " exists");
    if (dbDirIn == dbDirOut)
      throw runtime_error ("Input database directory = output database directory");

    prog2dir ["fasta_check"] = execDir;
    findProg ("makeblastdb");


    stderr << "Copying database to " << dbDirOut << " ";
    createDirectory (dbDirOut, false);
    {
      FileItemGenerator fileGen (0, true, dbDirIn);
      string fName;
      while (fileGen. next (fName))
      {
        stderr << ".";
        exec ("cp " + dbDirIn + "/" + fName + " " + dbDirOut + "/");
      }
      stderr << "\n";
    }
    

    // Old files
    stderr << "Reading " << dbDirOut << "/fam.tab ...\n";
    const TextTable fam (dbDirOut + "/fam.tab");
    fam. qc ();
    const TextTable::Key fam_key (fam, StringVector {"fam_id"});
    
    stderr << "Reading " << dbDirOut << "/taxgroup.tab ...\n";
    const TextTable taxgroup (dbDirOut + "/taxgroup.tab");
    taxgroup. qc ();
    const TextTable::Key taxgroup_key (taxgroup, StringVector {"taxgroup"});
      
    
    // Customization files
    stderr << "Reading " << metaFName << " ...\n";
    const TextTable meta (metaFName);
    meta. qc ();
    const TextTable::Key meta_key (meta, StringVector {"gene_symbol"});
      
    stderr << "Reading " << pmFName << " ...\n";
    const TextTable pm (pmFName);
    pm. qc ();

      
    Set<string> newGenesymbols;
      
      
    // AMRProt
    stderr << "Customizing " + dbDirOut + "/AMRProt ...\n";
	  exec (fullProg ("fasta_check") + protFName + " -aa" + qcS);  
    {
      ofstream fOut (dbDirOut + "/AMRProt", ios_base::app);
      LineInput fIn (protFName);
      Istringstream iss;
      StringVector fam_value (1);
      const Vector<size_t> pm_index (pm. columns2indexes (StringVector {"protein_identifier"})); 
      StringVector pm_value (1);
      while (fIn. nextLine ())
        if (   ! fIn. line. empty () 
            && fIn. line [0] == '>'
           )
        try
        {
          iss. reset (fIn. line. substr (1));
          string accession, genesymbol;
          iss >> accession >> genesymbol;
          if (accession. empty ())
            throw runtime_error ("Empty protein accession");
          if (genesymbol. empty ())
            throw runtime_error ("Empty gene symbol");
          string resistance, reportable, classS, subclass, product;
          bool isMutation = false;
          {
            fam_value [0] = genesymbol;
            pm_value  [0] = accession;
            const TextTable::RowNum fam_row_num  = fam_key.  find (fam_value);
            const TextTable::RowNum meta_row_num = meta_key. find (fam_value);
            if (pm. find (pm_index, pm_value, 0) != no_index)
              isMutation = true;
            if (fam_row_num != no_index)
            {
              if (meta_row_num != no_index)
                throw runtime_error ("fam.tab contains gene symbol " + genesymbol + " from " + metaFName);
              if (isMutation)
                throw runtime_error ("fam.tab contains a mutation gene symbol " + genesymbol);
              const StringVector& row = fam. rows [fam_row_num];
              reportable = row [fam. col2index ("reportable")];
              subclass   = row [fam. col2index ("subclass")];
              classS     = row [fam. col2index ("class")];
              product    = row [fam. col2index ("family_name")];
            }
            else if (meta_row_num != no_index)
            {
              if (isMutation)
                throw runtime_error (metaFName + " contains a mutation gene symbol " + genesymbol);
              newGenesymbols << genesymbol;
              const StringVector& row = meta. rows [meta_row_num];
              reportable = row [meta. col2index ("reportable")];
              subclass   = row [meta. col2index ("subclass")];
              classS     = row [meta. col2index ("class")];
              product    = row [meta. col2index ("protein_name")];
            }
            else if (isMutation)
            {
              resistance = "mutation";
              reportable = "2";
              product    = "WILDTYPE NAME";              
            }
            else
              throw runtime_error ("Unknown gene symbol: " + genesymbol);
          }
          if (   reportable != "1"
              && reportable != "2"
             )
            throw runtime_error ("Wrong reportable value: " + reportable);
          if (product. empty ())
            throw runtime_error ("Empty product name");
          replace (product, ' ', '_');
          fOut << ">0" 
               << '|' << accession 
               << "|1"
               << "|1"
               << '|' << genesymbol 
               << '|' << genesymbol 
               << '|' << resistance
               << '|' << reportable
               << '|' << subclass
               << '|' << classS
               << '|' << product
               << endl;
        }
        catch (const exception &e)
        {
          throw runtime_error (string (e. what ()) + "\nLine number " + to_string (fIn. lineNum) + " of " + protFName);
        }
        else
          fOut << fIn. line << endl;
    }
	  exec (fullProg ("makeblastdb") + " -in " + dbDirOut + "/AMRProt" + "  -dbtype prot  -logfile /dev/null");  
    
   
    stderr << "Customizing " + dbDirOut + "/fam.tab ...\n";
    {
      ofstream famF (dbDirOut + "/fam.tab", ios_base::app);
      StringVector fam_value (1);
      for (const string& genesymbol : newGenesymbols)
      {
        fam_value [0] = genesymbol;
        const TextTable::RowNum row_num = meta_key. find (fam_value);
        ASSERT (row_num != no_index);
        const StringVector& row = meta. rows [row_num];
        famF         << genesymbol
             << '\t' << "ALL"
             << '\t' << genesymbol
             << "\t-\t0\t0\t0\t0\t0\t0\t0\t0"
             << '\t' << row [meta. col2index ("reportable")]
             << '\t' << row [meta. col2index ("type")]
             << '\t' << row [meta. col2index ("subtype")]
             << '\t' << row [meta. col2index ("class")]
             << '\t' << row [meta. col2index ("subclass")]
             << '\t' << row [meta. col2index ("protein_name")]
             << endl;
      }
    }


    stderr << "Customizing " + dbDirOut + "/AMRProt-mutation.tab ...\n";
    {
      ofstream mutF (dbDirOut + "/AMRProt-mutation.tab", ios_base::app);
      ofstream taxgroupF (dbDirOut + "/taxgroup.tab", ios_base::app);
      StringVector values (1);
      Set<string> newOrganisms;
      for (const StringVector& row : pm. rows)
      {
        string organism = row [pm. col2index ("organism")];
        replace (organism, ' ', '_');
        try
        {
          string pos = row [pm. col2index ("mutation")];
          while (   ! pos. empty ()
                 && ! isDigit (pos. front ())
                )
            pos. erase (0, 1);
          while (   ! pos. empty ()
                 && ! isDigit (pos. back ())
                )
            pos. erase (pos. size () - 1, 1);
          if (pos. empty ())
            throw runtime_error ("No mutation position");
          string product = row [pm. col2index ("mutated_protein_name")];
          replace (product, ' ' , '_');
          if (product. empty ())
            throw runtime_error ("Empty product name");
          mutF         << organism 
               << '\t' << row [pm. col2index ("protein_identifier")]
               << '\t' << pos
               << '\t' << row [pm. col2index ("gene_symbol")] << '_' << row [pm. col2index ("mutation")]
               << '\t' << row [pm. col2index ("class")]
               << '\t' << row [pm. col2index ("subclass")]
               << '\t' << product
               << endl;
        }
        catch (const exception &e)
        {
          throw TextTable::Error (pm, string (e. what ()) + "\nIn row: " + row. toString (","));
        }
        values [0] = organism;        
        if (   taxgroup_key. find (values) == no_index
            && ! newOrganisms. contains (organism)
           )
        {
          taxgroupF << organism << '\t' << organism << "\t0" << endl;
          newOrganisms << organism;
        }
      }
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



