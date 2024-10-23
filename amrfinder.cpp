// amrfinder.cpp

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
*  Please cite NCBI in any work or product based on this material.   
*
* ===========================================================================
*
* Author: Vyacheslav Brover
*
* File Description:
*   AMRFinderPlus
*   https://github.com/ncbi/amr
*
* Dependencies: NCBI BLAST, HMMer, libcurl, gunzip (optional)
*
* Release changes:
*   4.0.2   10/23/2024 PD-5155  StxTyper version 1.0.27
*   4.0.1   10/22/2024 PD-5155  "::" is a fusion infix for the column "Hierarchy node"
*                      PD-5085  Change column "Element length" to "Target length"
*                      PD-5085  StxTyper version 1.0.26
*   4.0.0   10/22/2024 PD-5156  StxTyper etc.
*   3.13.3  09/11/2024          Instruction for -gff: "Locations are in the --nucleotide file."
*                               StxTyper version 1.0.25
*   3.13.2  08/16/2024 PD-5085  column names to match MicroBIGG-E
*   3.13.1  08/14/2024 PD-5084  AmrMutation(geneMutation_std,geneMutation_reported)
*                               point mutation files have columns: standard_mutation_symbol, reported_mutation_symbol
*                      PD-5091  files: .fa, .tsv
*   3.12.25 08/05/2024 PD-5076  StxTyper version 1.0.24
*                               colorize() is suppressed for redirected output
*   3.12.24            PD-5064  StxTyper version 1.0.23
*                               colorizeUrl()
*   3.12.23            PD-5054  Check that the file "AMR_DNA-" + organism1 + ".ndb" exists
*                      PD-5038  StxTyper version 1.0.21
*                      PD-4969  redundant QC check is removed
*                               -db_gencode parameter is added to tblastn
*   3.12.22 06/03/2024 PD-5015  HMMs are not applied to mutation proteins
*   3.12.21 05/31/2024 PD-5014  duplicate susceptible proteins are not reported
*   3.12.20 05/22/2024 PD-4078  a regular reference protein can have point mutations
*   3.12.19 05/21/2024 PD-5002  StxTyper 1.0.20
*           05/06/2024          BlastAlignment: isMutation() => !seqChanges.empty()
*                                               good() is not using refMutation
*                                               betterEq(): only within the same category: inFam(), isMutationProt(), isSusceptibleProt(), refMutation
*           04/30/2024 PD-4981  !isMutation() --> inFam() 
*           04/26/2024          StxTyper does not use --log parameter
*   3.12.18 03/26/2024          StxTyper ver. 1.0.19
*           03/13/2024 PD-4926  StxTyper ver. 1.0.16; amr_report.cpp reports all stx genes
*           03/11/2024 PD-4924  StxTyper 1.0.15: dead stxA2j EFK5907329.1 is replaced by EMA1832120.1
*   3.12.17 03/08/2024 PD-4925  more detailed explanation for the message "The BLAST database for AMRProt was not found."
*   3.12.16 03/08/2024 PD-4923  verify the version of StxTyper
*   3.12.15 03/05/2024 PD-4918  stxtyper ver. 1.0.14: --print_node 
*                               stxtyper -q
*   3.12.14 03/05/2024 PD-4910  -n <nucl> --plus -O Escherichia: stxtyper 1.0.13 is invoked, STX* classes are suppressed in amr_report.cpp
*   3.12.13 02/27/2024 PD-4888  mutation.cpp is added to the distribution (used only in testing)
*   3.12.12 02/21/2024 PD-4906  --print_node for fusion proteins
*                               move DOCUMENTATION section to the end of the usage message
*   3.12.11 02/13/2024 PD-4893  move DOCUMENTATION and UPDATES sections to the end of the --help message
*   3.12.10 02/12/2024 PD-4874  -v == --version
*   3.12.9  02/02/2024          moving some functions and variables to common.{hpp,cpp}
*   3.12.8  02/01/2024 PD-4872  ALLELEP = EXACTP for alleles; exact matches with possibly hanging tailes are preferred
*   3.12.7  02/01/2024 PD-4872  hanging tails of target protein are allowed for EXACTP matches, but ALLELEP requires EXACTP with no hanging tails
*           01/27/2024          replacing getFam() by getMatchFam() in amr_report.cpp
*   3.12.6  01/26/2024          memory leaks in amr_report.cpp
*   3.12.5  01/19/2024          BlastAlignment::BlastRule's are valid iff !fromHmm and !inFam()
*   3.12.4  01/18/2024 PD-4856  allow multiple Blast Rules for the same protein
*   3.12.3  01/12/2024          improved error message reporting for GFF files (https://github.com/ncbi/amr/issues/135)
*   3.12.2  12/21/2023 PD-4843  --mutation_all should report only point mutations
*   3.12.1  12/15/2023 PD-4838  stop codons are added to the reference proteins in AMRFinderPlus
*                               input proteins may miss '*' at the ends
*                               target hits have a new three-valued flag targetStopCodon: detected, missing, unknown
*                               procsssing of Prodigal GFF format is restored
*   3.11.26 10/16/2023 PD-4772  remove Prodigal GFF format from AMRFinderPlus
*   3.11.25 10/13/2023 PD-4771  revert removing '*' from Prodigal output to ensure ALLELEP and EXACTP matches
*   3.11.24 10/12/2023 PD-4769  --print_node prints FAM.id replaced by FAM.parent for non-exact allele matches
*   3.11.23 10/06/2023 PD-4764  remove '*' from Prodigal output to ensure ALLELEP and EXACTP matches
*           10/05/2023 PD-4761  remove protein sequences with >= 20 Xs
*   3.11.22 10/05/2023 PD-4754  Prodigal GFF
*   3.11.21 10/02/2023 PD-4755  bug: calling fusion2geneSymbols() for a mutation protein
*   3.11.20 09/06/2023 PD-4722  bug: calling fusion2geneSymbols() for a mutation protein
*                               color codes are printed only when output is to screen
*   3.11.19 08/09/2023 PD-4698  if a pseudogene is overlapped by a different gene on the length >= 20 aa with the same gene symbol then the pseudogene is not reported
*           08/04/2023 PD-4706  protein match overrides a nucleotide match for point mutations
*   3.11.18 07/25/2023          parameter order in instruction; "can be gzipped" is added to help
*   3.11.17 07/19/2023 PD-4687  distinct overlapping hits are not reported separately for protein targets for the same alleles or gene symbols
*   3.11.16 07/18/2023 PD-4687  distinct overlapping hits are not reported separately for protein targets (because the start/stop are not reported)
*!  3.11.15 05/23/2023 PD-4629  "amrfinder_update -d DIR" will create DIR if DIR is missing
*   3.11.14 05/06/2023 PD-4598  error messages in curl_easy.cpp
*   3.11.13 05/04/2023 PD-4596  prohibit ASCII characters only between 0x00 and 0x1F in GFF files
*           04/24/2023 PD-4583  process files ending with ".gz", see https://github.com/ncbi/amr/issues/61, dependence on gunzip (optional)
*           04/19/2023          On failure no empty output file (-o) is created
*   3.11.12 04/13/2023          application::makeKey()
*                      PD-4548  fasta_check.cpp prohibits '\t' (not any '\'), and all restrictions are only for nucleotide sequences
*   3.11.11 04/13/2023 PD-4566  --hmmer_bin
*   3.11.10 04/12/2023 PD-4548  fasta_check.cpp prohibits ';', '.', '~' in the last position of a sequence identifier
*                      PD-4548  fasta_check.cpp prohibits: ',,' and '\\' in all positions, '?' in initial position, and ',' in the last position of a sequence identifier
*   3.11.9  04/11/2023 PD-4560  BLAST -mt_mode is used on Mac only for BLAST version >= 2.13.0
*           04/05/2023 PD-4522  blastp -task blastp-fast
*           04/05/2023 PD-4548  "-a standard" is added
*   3.11.8  04/01/2023          fasta_extract.cpp checks whether all requested identifiers are found in FASTA
*   3.11.7  03/30/2023 PD-4548  GFF parsing processes '%<hex>' characters
*   3.11.6  03/21/2023 PD-4522  tblastn: -word_size 3 --> -task tblastn-fast -threshold 100  -window_size 15
*           03/21/2023 PD-4533  '_' are incorrectly trimmed from contig names
*           03/13/2023          -mt_mode is restored for __APPLE__ 
*   3.11.5  03/10/2023          directories --blast_bin, directory DATABASE in amrfinder_index.cpp do not need to end with '/'
*           03/01/2023 PD-3597  amrfinder_index
*           02/27/2023          section()
*   3.11.4  01/24/2023          GPipe organism string in taxgroup.tab is a comma-separated list of GPipe organisms
*   3.11.3  12/27/2022          "No valid AMRFinder database is found.\nThis directory (or symbolic link to directory) is not found: " + db
*   3.11.2  12/13/2022 PD-4427  a database of the older software minor is loaded for a new software minor version
*           12/05/2022          detect reference frameshited proteins
*   3.11.1  11/23/2022 PD-4414  modified reference proteins can have unequal lengths of reference and allele sequences
*           11/04/2022 PD-4394  --print_node
*           10/25/2022          setSymlink(,,bool); g++ -std=gnu++17
*   3.10.47 10/22/2022          exec(,logFile) is added
*           10/21/2022          simplified: "No valid AMRFinder database found.\nSymbolic link is not found: " + db
*   3.10.46 10/21/2022          blast error messages are printed if blast fails
*   3.10.45 10/21/2022 PD-4348  path2canonical() is sometimes needed for setSymlink(); exec() throw prints command and exit code
*   3.10.44 10/12/2022 PD-4348  "latest" symbolic link is relative, and it is updated by "amrfinder -u"
*   3.10.43 10/05/2022 PD-4341  non-existant GPipe taxgroup
*   3.10.42 10/03/2022 PD-4333  https://github.com/ncbi/amr/issues/99
*   3.10.41 08/19/2022          "-max_target_seqs 10000" is used in all BLAST commands
*   3.10.40 08/12/2022 PD-4297  duplicated rows in amr_report output
*   3.10.39 08/10/2022 PD-4290  wrong QC in amr_report.cpp is removed
*   3.10.38 08/04/2022 PD-4264  prohibition of the same contig and protein names is removed
*   3.10.37 08/02/2022 PD-4277  --annotation_format is restored, amr_report is faster, -a patric allows "accn|" in nucleotide FASTA, -a pseudomonasdb
*   3.10.36 08/02/2022 PD-4277  --annotation_format is removed 
*   3.10.35 08/01/2022 PD-4277  repeated allele is not reported correctly
*   3.10.34 08/01/2022 PD-4264  prokka and bakta have nucleotide FASTA in GFF files
*   3.10.33 07/29/2022 PD-4264  --annotation_format 
*   3.10.32 07/29.2022 PD-4274  "amrfinder -u" crashes if the directiry "data/" already exists
*                               GFF code refactoring
*           07/28/2022 PD-4274  "amrfinder -u" can use --blast_bin
*           07/28/2022 PD-3292  dependence on ls, rm, tail, cat is removed
*   3.10.31 07/27/2022          AMRProt.phr is checked earlier; locus_tagP is checked faster; fasta2parts is not using tmp.db
*   3.10.30 05/28/2022 PD-4217
*   3.10.29 05/27/2022 PD-4217  multi-domain tccP BLASTP result is confused with a fusion protein
*   3.10.28 05/11/2022 PD-4169  CDSs are the same if CDS difference is shorter than 60 aa
*   3.10.27 05/06/2022 PD-4119  --database_version
*   3.10.26 05/02/2022 PD-3292  dependence on "ln" is removed
*   3.10.25 04/29/2022 PD-3292  dependence on "mv", "cp", "cut", "head" and "sort" is removed
*   3.10.24 03/25/2022 PD-4132  pmrB_RPISLR6del shadows pmrB_L10P
*   3.10.23 02/11/2022 PD-4098  HTTPS connection for database downloading is restored (only this connection is guaranteed to exist for a user because it is needed for software installation)
*   3.10.22 02/10/2022 PD-4098  for FTP: FTP EPSV mode is turned off (PASV is turned on)
*   3.10.21 01/31/2022 PD-4071  https --> ftp
*   3.10.20 01/13/2022 PD-4069  "Vulnerability Disclosure" at the AMRFinder download web page
*   3.10.19 12/08/2021          $PATH is printed is case of error
*   3.10.18 11/01/2021 PD-4013  multi-domain protein hits
*   3.10.17 09/27/2021          program "dna_mutation" is found if the DNA file is empty
*   3.10.16 09/22/2021 PD-3958  point mutations oiverride HMM hits
*   3.10.15 08/23/2021          -mt_mode 1 is used only for multi-FASTA files: SB-3162
*   3.10.14 08/19/2021 PD-3918  BLAST output to stderr is not checked
*   3.10.13 08/19/2021 PD-3918  BLAST -mt_mode 1 is turned off for Mac, see SB-3163
*   3.10.12 08/19/2021 PD-3918  BLAST output to stderr is reported as error, except for BLASTN due to SB-3162
*   3.10.11 08/18/2021 PD-3826  dashes in a protein FASTA file are removed with a warning (crashes with some versions of HMMer)
*   3.10.10 08/16/2021 PD-3910  alien organism's proteins are removed from processing in amr_report.cpp (point mutations, susceptible)
*   3.10.9  08/13/2021 PD-3888  temporary files are named "amrfinder.XXXXXX"
*   3.10.8  07/06/2021 PD-3865  creating a directory does not break if upper-level directories are not readable
*                      PD-3867  "BLAST -mt_mode 1" is used
*                               makeblastdb and tblastn are used for large sequences 
*                               hmmsearch and tbasltn: query file is split into threads_max files
*                               "blastp -seg no" is used
*   3.10.7  05/18/2021 PD-3820  message for missing AMRProt blast database; https://github.com/ncbi/amr/issues/53
*   3.10.6  05/07/2021 PD-3796  for POINTN reported "target length" column = targetEnd - targetStart 
*   3.10.5  04/12/2021 PD-3772  --report_equidistant
*   3.10.4  03/24/2021 PD-3761  amrfinder --help will not break if /tmp is full
*   3.10.3  03/15/2021 PD-3749  --nucleotide_flank5_output, --nucleotide_flank5_size
*   3.10.2  03/03/2021 PD-3729  neighboring point mutations are reported
*   3.10.1  02/17/2021 PD-3679  AMRProt-susceptible.tab
*   3.9.10  02/16/2021 PD-3694  message about missing "latest/" symbolic link; amrfinder_update.cpp: createLatestLink()
*   3.9.9   01/27/2021 PD-3674  crash for a custom database
*   3.9.8   01/04/2021 PD-3613  --dir is removed
*   3.9.7   12/03/2020 PD-3292  dependence on "mkdir" is removed
*   3.9.6   11/20/2020 PD-3613  --dir
*                               prepare_fasta_extract()
*   3.9.5   11/18/2020 PD-3292  dependence on awk is removed
*                               --help prints instruction on $TMPDIR
*   3.9.4   11/16/2020 PD-3609  ($TMPDIR or "/tmp") + "/XXXXXX"
*   3.9.3   11/05/2020 PD-3577  merge lines for bifunctional proteins
*   3.9.2   11/04/2020 PD-3590  AMRProt has new fields #9 and #10: "subclass" and "class"
*   3.9.1   10/27/2020 PD-3583  AMRProt has a new field #8 "reportable"
*           09/30/2020 PD-2407  option --type is removed
*   3.8.28  09/29/2020 PD-3292  dependence on "uniq" is removed
*   3.8.27  09/28/2020 PD-2381  non-standard start codons are not changed in fusion proteins
*   3.8.26  09/25/2020 PD-2381  proteins with non-standard start codons that are extended in the N-terminal direction are EXACTP
*   3.8.25  09/25/2020 PD-3547  identification of frameshifts is disabled
*                               POINTX method with more SNPs is preferred over POINTP method
*   3.8.24  09/21/2020 PD-3536  --pointmut_all reports all SNPs in a reference gene repetition
*   3.8.23  09/16/2020 PD-3536  simplifying point mutations preference
*   3.8.22  09/15/2020 PD-3470  frameshift detection bug; preference of point mutation reference proteins 
*   3.8.21  09/14/2020 PD-3536  point mutations merging bug
*                      PD-3469  --force_update implies --update; -U
*   3.8.20  09/14/2020 PD-3531  "--parm -print_fam" bug
*   3.8.19  09/04/2020 PD-3292  removed the dependence on "grep"
*   3.8.18  09/03/2020 PD-3292  removed the dependence on "which"
*   3.8.17  09/02/2020 PD-3528  ordering of rows in the report is broken with parameter --name
*   3.8.16  09/01/2020 PD-2322  a complete nucleotide hit is not preferred to a partial protein hit; stopCodon field is borrowed from BLASTX to BPASTP
*   3.8.15  08/28/2020 PD-3475  return BLAST alignment parameters for HMM-only hits where available
*   3.8.14  08/27/2020 PD-3470  method FRAME_SHIFT, amr_report is faster
*   3.8.13  08/25/2020 PD-2322  a complete nucleotide hit is preferred to a partial protein hit
*   3.8.12  08/24/2020 PD-2394  fusion genes are reported to include both gene symbols on each line
*   3.8.11  08/21/2020 PD-2407  --type
*   3.8.10  08/20/2020 PD-3469  --force_update
*   3.8.9   08/13/2020          BLAST -show_gis parameter is removed, more mutations are reported for --mutation_all
*   3.8.8   08/04/2020          bug in fasta_extract.cpp, more output in --nucleotide_output
*   3.8.7   08/03/2020 PD-3504  --protein_output, --nucleotide_output options by fasta_extract.cpp
*   3.8.6   07/29/2020 PD-3468  --name option
*           07/13/2020 PD-3484  -l for old database versions
*   3.8.5   07/10/2020 PD-3482  --ident_min instruction
*   3.8.4   05/13/2020 PD-3447  custom point mutation does not match the reference sequence
*                               text "*** ERROR ***" is not repeated twice
*   3.8.3   05/01/2020          WILDTYPE mutations were reported as 0-based
*   3.8.2   05/01/2020 PD-3419  taxgroup is removed from the DNA files, dna_mutation parameter: organism
*                      PD-3437  --mutation_all requires --organism
*                               all warnings are printed to stderr
*                               warnings are printed in bright yellow color; ERROR is printed in bright red color
*                      PD-3363  WILDTYPE mutations map on the reference gene with offset
*                               NOVEL is changed to UNKNOWN
*   3.8.1   04/30/2020 PD-3419  dna_mutation: reporting gene symbol for novel mutations; taxgroup and genesymbol are added to the DNA files
*   3.7.6   04/29/2020 PD-3419  dna_mutation: reporting gene symbol for novel mutations
*   3.7.5   04/22/2020 PD-3427  -h prints the help message
*   3.7.4   04/14/2020 PD-3391  Mac Conda installation
*   3.7.3   04/09/2020 PD-3416  redundant QC check in alignment.cpp
*   3.7.2   04/08/2020 PD-3363  "WILDTYPE" was not reported
*   3.7.1   04/02/2020 PD-3154  GIs may be 0, accessions are main identifiers; file "AMRProt-suppress" is added accessions; dataVer_min is "2020-04-02.1"
*   3.6.19  03/24/2020          Check of ">lcl|" is done only for the first sequence in FASTA
*           03/24/2020 PD-3347  -lcl parameter in gff_check and amr_report
*   3.6.18  03/17/2020 PD-3396  amr_report.cpp prints a better error message on missing sublcass in data
*   3.6.17  03/12/2020          software version is printed after software directory
*   3.6.16  03/06/2020 PD-3363  --mutation_all: UNKNOWN are not reported
*                      PD-2328  last 2 columns of report are real HMM hits
*   3.6.15  02/24/2020          "database" is printed to stderr in one line in a canonical form (without links)
*   3.6.14  02/19/2020 PD-3363   --mutation_all: 4 types of mutations, adding DNA mutations
*   3.6.13  02/13/2020 PD-3359,issue#23   ln -s <db>: uses path2canonical()
*   3.6.12  02/13/2020 PD-3359,issue#23   AMRFinder database directory may contain spaces
*   3.6.11  02/13/2020 PD-3359,issue#23   AMRFinder code directory may contain spaces
*   3.6.10  02/06/2020 PD-3357,issue#21  --mutation_all bug
*           01/24/2020 PD-3345   improved error message for "GFF file mismatch"
*   3.6.9   01/13/2020           "Database directory" is printed to stederr
*           01/10/2020 PD-3329   ln -s .../amrfinder abc: abc calls the right executables
*           01/20/2020           'rm" dependence is removed
*   3.6.8   01/10/2020           'gnl|' processing is simplified
*           01/09/2020 PD-3327   allow empty input files
*   3.6.7   01/09/2020           do not remove 'lcl|' from DNA FASTA
*   3.6.6   01/09/2020 PD-3326   'gnl|' is added only for gnl|PROJECT|ACC accessions if --pgap
*           01/09/2020 PD-3324   pipefail requires bash
*           01/08/2020 GP-28123  'gnl|' is added to report if --pgap
*   3.6.5                        sorting of reported rows: gene symbol is used as the last sorting column if contig is available
*   3.6.4   01/03/2020 PD-3230   sorting of reported rows: protein accession is ignored if contig is available
*   3.6.3   01/03/2020 PD-3230   sorting of reported rows
*           12/28/2019           QC in dna_mutation
*   3.6.2   12/27/2019 PD-3230   redundant reported lines are removed for mutated reference proteins
*                                reports are sorted by sort
*   3.6.1   12/27/2019 PD-3230   mutated proteins are added to AMRProt
*   3.5.10  12/20/2019           --log
*   3.5.9   12/19/2019 PD-3294   blastx parameters: space added
*   3.5.8   12/18/2019 issues/19 changed message if db path is bad
*   3.5.7   12/18/2019 PD-3289   improved message for gff_check failure
*   3.5.6   12/18/2019 PD-3269   --gpipe is removed, --pgapx is replaced by --pgap
*   3.5.5   12/17/2019 PD-3287   short proteins at an end of a contig are reported
*   3.5.4   12/17/2019 PD-3287   truncated short proteins are not reported
*   3.5.3   12/16/2019 PD-3279   GPipe-GenColl assemblies, --gpipe_org
*                      GP-28025
*   3.5.2   12/13/2019 PD-3269   new flag --pgapx
*   3.5.1   12/12/2019 PD-3277   files AMRProt-mutation.tab, AMRProt-suppress, AMR_DNA-<TAXGROUP>.tab and taxgroup.tab have headers
*   3.4.3   12/11/2019 PD-2171   --mutation_all bug
*                                --debug does not imply "-verbose 1"
*   3.4.2   12/10/2019 PD-3209   alignment correction for mutations
*                                point_mut.{hpp,cpp} -> alignment.{hpp,cpp}
*                                dna_point_mut.cpp -> dna_mutation.cpp
*                                AMRProt-point_mut.tab -> AMRProt-mutation.tab
*                                protein resistance: "point_mutation" -> "mutation"
*                                amrfinder: --point_mut_all -> --mutation_all
*                      PD-3232   mutation detection redesign
*                      PD-3267   mutation in a mutated context
*   3.4.1   12/03/2019 PD-3193   AMR_DNA-*.tab: column "genesymbol" is removed
*                                product name is fixed for point mutations
*                                point_mut.cpp -> dna_point_mut.cpp
*   3.3.2   11/26/2019 PD-3193   indel mutations: partially implemented
*                                bug fixed: protein point mutations were reported incorrectly if there was an offset w.r.t. the reference sequence
*                                files AMRProt-point_mut.tab and AMR_DNA-<taxgroup>.tab: columns allele, symbol are removed
*                                files taxgroup.list and gpipe.tab are replaced by taxgroup.tab
*   3.3.1   11/22/2019 PD-3206   new files: taxgroup.list, gpipe.tab; new option --list_organisms
*   3.2.4   11/15/2019 PD-3191   dna_mutation.cpp: neighborhoodMismatch <= 0.04; good(): length >= min (refLen, 2 * flankingLen + 1)
*   3.2.3   11/14/2019 PD-3192   fixed error made by PD-3190
*   3.2.3   11/13/2019 PD-3190   organisms for --gpipe
*   3.2.3   11/12/2019 PD-3187   sequence name is always from AMRProt, not from fam.tab
*   3.2.2   11/06/2019 PD-2244   added "LANG=C" before "sort"
*
*/


#ifdef _MSC_VER
  #error "UNIX is required"
#endif
   
#undef NDEBUG 

#include <unistd.h>

#include "common.hpp"
#include "tsv.hpp"
using namespace Common_sp;
#include "gff.hpp"
using namespace GFF_sp;
#include "columns.hpp"

#include "common.inc"



#undef DIR  // PD-3613


// PAR!
// PD-3051
const string dataVer_min ("2024-08-14.2");
  // 3.12: "2023-12-15.2"
  // 3.11: "2021-02-18.1"  
const string stxTyperVersion ("1.0.27");  



namespace 
{


// PAR
//constexpr size_t threads_max_min = 1;  
constexpr size_t threads_def = 4;
// Cf. amr_report.cpp
constexpr double ident_min_def = 0.9;
constexpr double partial_coverage_min_def = 0.5;
const string ambigS ("20");  


struct ThisApplication : ShellApplication
{
  ThisApplication ()
    : ShellApplication ("Identify AMR and virulence genes in proteins and/or contigs and print a report", true, true, true, true)
    {
    	addFlag ("update", "Update the AMRFinder database", 'u');  // PD-2379
    	addFlag ("force_update", "Force updating the AMRFinder database", 'U');  // PD-3469

    #ifdef DIR
      addKey ("dir", "Common directory of the --protein, --nucleotide and --gff files", "", '\0', "DIRECTORY");
    #endif
    	addKey ("protein", "Input protein FASTA file (can be gzipped)", "", 'p', "PROT_FASTA");
    	addKey ("nucleotide", "Input nucleotide FASTA file (can be gzipped)", "", 'n', "NUC_FASTA");
    	addKey ("gff", "GFF file for protein locations (can be gzipped). Locations are in the --nucleotide file. Protein ids should be in the attribute 'Name=<id>' (9th field) of the rows with type 'CDS' or 'gene' (3rd field).", "", 'g', "GFF_FILE");

      {    	
      	const string annots (Gff::names. toString (", "));
        addKey ("annotation_format", "Type of GFF file: " + annots, "genbank", 'a', "ANNOTATION_FORMAT");  
      }

    	addKey ("database", "Alternative directory with AMRFinder database. Default: $AMRFINDER_DB", "", 'd', "DATABASE_DIR");
    	addFlag ("database_version", "Print database version", 'V');

    	addKey ("ident_min", "Minimum proportion of identical amino acids in alignment for hit (0..1). -1 means use a curated threshold if it exists and " + toString (ident_min_def) + " otherwise", "-1", 'i', "MIN_IDENT");
    	  // "PD-3482
    	addKey ("coverage_min", "Minimum coverage of the reference protein (0..1)", toString (partial_coverage_min_def), 'c', "MIN_COV");

      addKey ("organism", "Taxonomy group. To see all possible taxonomy groups use the --list_organisms flag", "", 'O', "ORGANISM");
      addFlag ("list_organisms", "Print the list of all possible taxonomy groups for mutations identification and exit", 'l');
    	addKey ("translation_table", "NCBI genetic code for translated BLAST", "11", 't', "TRANSLATION_TABLE");

    	addFlag ("plus", "Add the plus genes to the report");  // PD-2789
      addFlag ("report_common", "Report proteins common to a taxonomy group");  // PD-2756
      addFlag ("report_all_equal", "Report all equally-scoring BLAST and HMM matches");  // PD-3772
      addKey ("name", "Text to be added as the first column \"name\" to all rows of the report, for example it can be an assembly name", "", '\0', "NAME");
      addFlag ("print_node", "Print hierarchy node (family)");  // PD-4394
    	addKey ("mutation_all", "File to report all mutations", "", '\0', "MUT_ALL_FILE");

      addKey ("output", "Write output to OUTPUT_FILE instead of STDOUT", "", 'o', "OUTPUT_FILE");
      addKey ("protein_output", "Output protein FASTA file of reported proteins", "", '\0', "PROT_FASTA_OUT");
      addKey ("nucleotide_output", "Output nucleotide FASTA file of reported nucleotide sequences", "", '\0', "NUC_FASTA_OUT");
      addKey ("nucleotide_flank5_output", "Output nucleotide FASTA file of reported nucleotide sequences with 5' flanking sequences", "", '\0', "NUC_FLANK5_FASTA_OUT");
      addKey ("nucleotide_flank5_size", "5' flanking sequence size for NUC_FLANK5_FASTA_OUT", "0", '\0', "NUC_FLANK5_SIZE");

    	addKey ("blast_bin", "Directory for BLAST. Deafult: $BLAST_BIN", "", '\0', "BLAST_DIR");
    	addKey ("hmmer_bin", "Directory for HMMer", "", '\0', "HMMER_DIR");

      addFlag ("pgap", "Input files PROT_FASTA, NUC_FASTA and GFF_FILE are created by the NCBI PGAP");  // = --annotation_format pgap 
      addFlag ("gpipe_org", "NCBI internal GPipe organism names");

    	addKey ("parm", "amr_report parameters for testing: -nosame -noblast -skip_hmm_check -bed", "", '\0', "PARM");

	    version = SVN_REV;  
	    documentationUrl = "https://github.com/ncbi/amr/wiki";
	    updatesUrl = "https://www.ncbi.nlm.nih.gov/mailman/listinfo/amrfinder-announce";
	    updatesDoc = "\
    Subscribe to the amrfinder-announce mailing list for database and software update notifications";
    }



  void initEnvironment () final
  {
    ShellApplication::initEnvironment ();
    var_cast (name2arg ["threads"] -> asKey ()) -> defaultValue = to_string (threads_def);  
  }
  
  
  
  void fastaCheck (const string &fName, 
                   bool prot, 
                   const string &qcS, 
                   const string &logFName, 
                   size_t &nSeq, 
                   size_t &len_max,
                   size_t &len_total,
                   const string &outFName) const
  // Input: fName, outFName: quoted
  {
    ASSERT (fName != logFName);
    if (! outFName. empty ())
    {
      ASSERT (outFName != fName);
      ASSERT (outFName != logFName);
    }
    exec (fullProg ("fasta_check") + fName + "  " + (prot ? "-aa  -stop_codon  -ambig_max " + ambigS + prependS (outFName, "  -out ") : "-len " + tmp + "/len  -hyphen  -ambig") + qcS + "  -log " + logFName + " > " + tmp + "/nseq", logFName); 
      // "-stop_codon" PD-4771 

  	const StringVector vec (tmp + "/nseq", (size_t) 10, true); 
  	if (vec. size () != 3)
      throw runtime_error (string (prot ? "Protein" : "DNA") + " fasta_check failed: " + vec. toString ("\n"));
    nSeq      = str2<size_t> (vec [0]);
    len_max   = str2<size_t> (vec [1]);
    len_total = str2<size_t> (vec [2]);
    QC_ASSERT (nSeq);
    QC_ASSERT (len_max);
    QC_ASSERT (len_total);
  }

  
  
  StringVector db2organisms () const
  {
    const TextTable taxgroup            (tmp + "/db/taxgroup.tsv");
    const TextTable AMRProt_mutation    (tmp + "/db/AMRProt-mutation.tsv");
    const TextTable AMRProt_susceptible (tmp + "/db/AMRProt-susceptible.tsv");
    taxgroup. qc ();
    AMRProt_mutation. qc ();
    AMRProt_susceptible. qc ();
    StringVector vec (taxgroup. col2values (0));
    vec << AMRProt_mutation. col2values (0);
    vec << AMRProt_susceptible. col2values (0);
    vec. sort ();
    vec. uniq ();
    return vec;
  }
  
  
  
  void prepare_fasta_extract (const StringVector &columns,
                              const string &tmpSuf,
                              bool saveHeader) const
  // Input: tmp + "/amr"
  {
    TextTable t (tmp + "/amr");
    t. qc ();
    t. filterColumns (columns);
    t. rows. filterValue ([] (const StringVector& row) { return row [0] == "NA"; });
    t. rows. sort ();
    t. rows. uniq ();
    t. saveHeader = saveHeader;
    t. qc ();
    t. saveFile (tmp + "/" + tmpSuf);
  }
  
  
  
  void shellBody () const final
  {
    const bool    force_update    =             getFlag ("force_update");
    const bool    update          =             getFlag ("update") || force_update;
    const string  dir             =    
  #ifdef DIR
                                    appendS (getArg ("dir"), "/"); 
  #else
                                    "";
  #endif
    const string  prot             = shellQuote (prependS (getArg ("protein"),    dir));
    const string  dna              = shellQuote (prependS (getArg ("nucleotide"), dir));
    const string  gff              = shellQuote (prependS (getArg ("gff"),        dir));
          string  db               =             getArg ("database");
    const bool    pgap             =             getFlag ("pgap");
          Gff::Type gffType        = Gff::name2type (getArg ("annotation_format"));  
    const double  ident            =             arg2double ("ident_min");
    const double  cov              =             arg2double ("coverage_min");
    const string  organism         = shellQuote (getArg ("organism"));   
    const bool    list_organisms   =             getFlag ("list_organisms");
    const uint    gencode          =             arg2uint ("translation_table"); 
    const bool    add_plus         =             getFlag ("plus");
    const bool    report_common    =             getFlag ("report_common");
    const string  mutation_all     =             getArg ("mutation_all");  
          string  blast_bin        =             getArg ("blast_bin");
          string  hmmer_bin        =             getArg ("hmmer_bin");
    const bool    equidistant      =             getFlag ("report_all_equal");
    const bool    print_node       =             getFlag ("print_node");
    const string  input_name       = shellQuote (getArg ("name"));
    const string  parm             =             getArg ("parm");  
    const string  output           =             getArg ("output");
    const string  prot_out         = shellQuote (getArg ("protein_output"));
    const string  dna_out          = shellQuote (getArg ("nucleotide_output"));
    const string  dnaFlank5_out    = shellQuote (getArg ("nucleotide_flank5_output"));
    const uint    dnaFlank5_size   =             arg2uint ("nucleotide_flank5_size");
    const bool    gpipe_org        =             getFlag ("gpipe_org");
    const bool    database_version =             getFlag ("database_version");
    
    
		const string logFName (tmp + "/log");  // Command-local log file


    if (database_version)
      cout   << "Software directory: " << shellQuote (execDir) << endl;
    else
      stderr << "Software directory: " << shellQuote (execDir) << '\n';
    if (database_version)
      cout   << "Software version: " << version << endl; 
    else
	    stderr << "Software version: " << version << '\n'; 
    
    if (contains (input_name, '\t'))
      throw runtime_error ("NAME cannot contain a tab character");

  #if 0
    if (threads_max < threads_max_min)
      throw runtime_error ("Number of threads cannot be less than " + to_string (threads_max_min));
  #endif
    
		if (ident != -1.0 && (ident < 0.0 || ident > 1.0))
		  throw runtime_error ("ident_min must be between 0 and 1");
		
		if (cov < 0.0 || cov > 1.0)
		  throw runtime_error ("coverage_min must be between 0 and 1");
		  
	  if (report_common && emptyArg (organism))
		  throw runtime_error ("--report_common requires --organism");
	  if (report_common && ! add_plus)
		  throw runtime_error ("--report_common requires --plus");
		  
		  		  
		// PD-3437
	  if (! mutation_all. empty () && emptyArg (organism))
	  {
	    const Warning warning (stderr);
		  stderr << "--mutation_all option used without -O/--organism option. No point mutations will be screened";
		}

    OFStream::prepare (output);

    
    string defaultDb;
    #ifdef CONDA_DB_DIR
    // we're in condaland
      if (const char* s = getenv("CONDA_PREFIX")) {
        defaultDb = string (s) + "/share/amrfinderplus/data/latest";
      } else if (const char* s = getenv("PREFIX")) {
        const Warning warning (stderr);
        stderr << "This was compiled for running under bioconda, but $CONDA_PREFIX was not found" << '\n';
        defaultDb = string (s) + "/share/amrfinderplus/data/latest";
        stderr << "Reverting to $PREFIX: " << defaultDb;
      } else {
        const Warning warning (stderr);
        stderr << "This was compiled for running under bioconda, but $CONDA_PREFIX was not found" << '\n';
        stderr << "Reverting to hard coded directory: " << CONDA_DB_DIR "/latest";
        defaultDb = CONDA_DB_DIR "/latest";
      }
    #else
    // not in condaland
      defaultDb = execDir + "data/latest";
    #endif
    ASSERT (isRight (defaultDb, "/latest"));
        
		// db
		if (db. empty ())
		{
    	if (const char* s = getenv ("AMRFINDER_DB"))
    		db = s;
    	else
			  db = defaultDb;
		}
		ASSERT (! db. empty ());		  


    // blast_bin
    if (blast_bin. empty ())
    	if (const char* s = getenv ("BLAST_BIN"))
    		blast_bin = string (s);
    if (! blast_bin. empty ())
    {
	    addDirSlash (blast_bin);
	    prog2dir ["blastp"]      = blast_bin;
	    prog2dir ["blastx"]      = blast_bin;
	    prog2dir ["tblastn"]     = blast_bin;
	    prog2dir ["blastn"]      = blast_bin;
      prog2dir ["makeblastdb"] = blast_bin;
	  }

    if (! hmmer_bin. empty ())
    {
	    addDirSlash (hmmer_bin);
	    prog2dir ["hmmsearch"] = hmmer_bin;
	  }
	  

		if (update)
    {
      // PD-2447
      if (! emptyArg (prot) || ! emptyArg (dna))
        throw runtime_error ("AMRFinder -u/--update option cannot be run with -n/--nucleotide or -p/--protein options");
      if (! getArg ("database"). empty ())
        throw runtime_error ("AMRFinder update option (-u/--update) only operates on the default database directory. The -d/--database option is not permitted");
      if (getenv ("AMRFINDER_DB"))
      {
        const Warning warning (stderr);
        stderr << "AMRFINDER_DB is set, but AMRFinder auto-update only downloads to the default database directory";
        db = defaultDb;
      }
  		const Dir dbDir (db);
      if (! dbDir. items. empty () && dbDir. items. back () == "latest")
      {
        prog2dir ["amrfinder_update"] = execDir;
  		  exec (fullProg ("amrfinder_update") + " -d " + shellQuote (dbDir. getParent ()) + ifS (force_update, " --force_update") 
  		          + makeKey ("blast_bin", blast_bin)  
  		          + makeKey ("hmmer_bin", hmmer_bin)  
  		          + ifS (getQuiet (), " -q") + ifS (qc_on, " --debug") + " > " + logFName, logFName);
      }
      else
      {
        const Warning warning (stderr);
        stderr << "Updating database directory works only for databases with the default data directory format." << '\n'
               << "         Please see " + colorizeUrl ("https://github.com/ncbi/amr/wiki", ! isRedirected (cerr)) + " for details." << "\n"
               << "         Current database directory is: " << dbDir. get () << "\n"
               << "         New database directories will be created as subdirectories of " << dbDir. getParent ();
      }
		}


    const string downloadLatestInstr ("\nTo download the latest version to the default directory run: amrfinder -u");
    
		if (! directoryExists (db))  // PD-2447
		//throw runtime_error ("No valid AMRFinder database found: " + db + ifS (! update, downloadLatestInstr));
		  throw runtime_error ("No valid AMRFinder database is found.\nThis directory (or symbolic link to directory) is not found: " + db + ifS (! update, downloadLatestInstr));
    if (database_version)
      cout   << "Database directory: " << shellQuote (path2canonical (db)) << endl;
    else
		  stderr << "Database directory: " << shellQuote (path2canonical (db)) << '\n';
    setSymlink (db, tmp + "/db", true);

    {
      // PD-4925
      const string dbTest (db + "/AMRProt.fa.phr");
  		if (! fileExists (dbTest))
  			throw runtime_error ("The BLAST database for AMRProt.fa (" + dbTest + ") was not found.\nUse amrfinder -u or amrfinder --force_update to download and prepare database for AMRFinderPlus");
    }


		// PD-3051
		{
  	  istringstream versionIss (version);
  		const SoftwareVersion softwareVersion (versionIss);
  		const SoftwareVersion softwareVersion_min (db + "/database_format_version.txt"); 
  	//stderr << "Software version: " << softwareVersion. str () << '\n'; 
  		const DataVersion dataVersion (db + "/version.txt");
  		istringstream dataVersionIss (dataVer_min); 
  		const DataVersion dataVersion_min (dataVersionIss);  
      if (database_version)
        cout   << "Database version: " << dataVersion. str () << endl;
      else
        stderr << "Database version: " << dataVersion. str () << '\n';
      if (softwareVersion < softwareVersion_min)
        throw runtime_error ("Database requires software version at least " + softwareVersion_min. str ());
      if (dataVersion < dataVersion_min)
        throw runtime_error ("Software requires database version at least " + dataVersion_min. str () + downloadLatestInstr);
      if (database_version)
        return;
    }
    
    
    if (list_organisms)
    {
      const StringVector organisms (db2organisms ());
      cout << endl << "Available --organism options: " + organisms. toString (", ") << endl;
      return;
    }    		  

		  
    {
      string searchMode;
      StringVector includes;
      if (emptyArg (prot))
      {
        if (emptyArg (dna))
        {
          if (update)
            return;
  	  	  throw runtime_error ("Parameter --protein or --nucleotide must be present");
    		}
        else
        {
      		if (! emptyArg (gff))
            throw runtime_error ("Parameter --gff is redundant");
          searchMode = "translated nucleotide";
        }
      }
      else
      {
        searchMode = "protein";
        if (emptyArg (dna))
        {
          searchMode += "-only";
          includes << key2shortHelp ("nucleotide") + " and " + key2shortHelp ("gff") + " options to add translated searches";
        }
        else
        {
      		if (emptyArg (gff))
            throw runtime_error ("If parameters --protein and --nucleotide are present then parameter --gff must be present");
          searchMode = "combined translated and protein";
        }
      }
      if (emptyArg (prot) && ! emptyArg (prot_out))
        throw runtime_error ("Parameter --protein must be present for --protein_output");
      if (emptyArg (dna) && ! emptyArg (dna_out))
        throw runtime_error ("Parameter --nucleotide must be present for --nucleotide_output");
      if (emptyArg (dna) && ! emptyArg (dnaFlank5_out))
        throw runtime_error ("Parameter --nucleotide must be present for --nucleotide_flank5_output");
      if (emptyArg (dnaFlank5_out) && dnaFlank5_size > 0)
        throw runtime_error ("Parameter --nucleotide_flank5_output must be present for --nucleotide_flank5_size");
      if (! emptyArg (dnaFlank5_out) && dnaFlank5_size == 0)
        throw runtime_error ("Parameter --nucleotide_flank5_size must be present with a positive value for --nucleotide_flank5_output");
      ASSERT (! searchMode. empty ());
      if (emptyArg (organism))
        includes << key2shortHelp ("organism") + " option to add mutation searches and suppress common proteins";
      else
        searchMode += " and mutation";
      stderr << "AMRFinder " << searchMode << " search\n";

      for (const string& include : includes)
        stderr << "  - include " << include << '\n';
    }
    
    
    // Quoted names
    const string prot_flat = uncompress (prot, "prot_flat");
    const string dna_flat  = uncompress (dna,  "dna_flat");
    const string gff_flat  = uncompress (gff,  "gff_flat");
      

    {
      StringVector emptyFiles;
      if (! emptyArg (prot) && ! getFileSize (unQuote (prot_flat)))  emptyFiles << prot;
      if (! emptyArg (dna)  && ! getFileSize (unQuote (dna_flat)))   emptyFiles << dna;
      if (! emptyArg (gff)  && ! getFileSize (unQuote (gff_flat)))   emptyFiles << gff;      
      for (const string& emptyFile : emptyFiles)
      {
        const Warning warning (stderr);
        stderr << "Empty file: " << emptyFile;
      }
    }


	  // organism --> organism1
	  string organism1;
	  bool suppress_common = false;	  
	  if (! emptyArg (organism))
	  {
	  	organism1 = unQuote (organism);
 	  	replace (organism1, ' ', '_');
 	  	ASSERT (! organism1. empty ());
      if (gpipe_org)
      {
        LineInput f (db + "/taxgroup.tsv");
        Istringstream iss;
        bool found = false;
        while (f. nextLine ())
        {
	  	    if (isLeft (f. line, "#"))
	  	      continue;
          iss. reset (f. line);
          string org, gpipeOrgs;
          int num = -1;
          iss >> org >> gpipeOrgs >> num;
          QC_ASSERT (! org. empty ());
          QC_ASSERT (num >= 0);
          QC_ASSERT (iss. eof ());
          const StringVector gpipeOrgVec (gpipeOrgs, ',', true);
          QC_ASSERT (gpipeOrgVec. size () >= 1);
          if (gpipeOrgVec. contains (organism1))
          {
            organism1 = org;
            found = true;
            break;
          }
        }
        if (! found)  // PD-4341
        #if 0
          throw runtime_error ("Non-existant GPipe taxgroup: " + organism);  
        #else
        {
    	    const Warning warning (stderr);
    		  stderr << "Non-existant GPipe taxgroup: " << organism;
    		  organism1. clear ();
        }
        #endif
      }
 	  }
 	  ASSERT (! contains (organism1, ' '));

	  if (! organism1. empty ())
	  {
      const StringVector organisms (db2organisms ());
      if (! organisms. contains (organism1))
        throw runtime_error ("Possible organisms: " + organisms. toString (", "));  
 	  	if (! report_common)
 	  	  suppress_common = true;
 	  }


    const string qcS (qc_on ? " -qc" : "");
		const string force_cds_report (! emptyArg (dna) && ! organism1. empty () ? "-force_cds_report" : "");  // Needed for dna_mutation
		
								  
    prog2dir ["fasta_check"]   = execDir;
    prog2dir ["fasta2parts"]   = execDir;
    prog2dir ["amr_report"]    = execDir;	
		prog2dir ["dna_mutation"]  = execDir;
    prog2dir ["fasta_extract"] = execDir;
		prog2dir ["stxtyper"]      = execDir + "stxtyper/";
    

    if (pgap)
    {
      switch (gffType)
      {
        case Gff::genbank: gffType = Gff::pgap; break;
        case Gff::pgap: break;
        default: throw runtime_error ("--pgap conflicts with GFF type " + strQuote (Gff::names [(size_t) gffType]));
      }
    }
    
    bool lcl = false;
    if (gffType == Gff::pgap && ! emptyArg (dna))  // PD-3347
    {
      LineInput f (unQuote (dna_flat));
      while (f. nextLine ())
        if (isLeft (f. line, ">"))
        {
          lcl = isLeft (f. line, ">lcl|");
          break;
        }
    }
    

	  const bool blastn = ! emptyArg (dna) && ! organism1. empty () && fileExists (db + "/AMR_DNA-" + organism1 + ".fa");
		const bool stxTyper = blastn && organism1 == "Escherichia" && add_plus;
		
		
		if (blastn)
    {
      // PD-5054
      const string dbTest (db + "/AMR_DNA-" + organism1 + ".fa.ndb");
  		if (! fileExists (dbTest))
  			throw runtime_error ("The BLAST database for AMR_DNA-" + organism1 + ".fa was not found.\nUse amrfinder -u or amrfinder --force_update to download and prepare database for AMRFinderPlus");
    }

		if (stxTyper)
		{
		  const string verFName (tmp + "/stxtyper-ver");
			exec (fullProg ("stxtyper") + " -v > " + verFName, logFName);
	    LineInput f (verFName);
	    if (! f. nextLine ())
	      throw runtime_error ("Cannot get the version of StxTyper");
	    if (f. line != stxTyperVersion)
	      throw runtime_error ("AMRFinder invokes StxTyper version " + f. line + ". Expected StxTyper version is " + stxTyperVersion);		  
		}


    // Create files for amr_report    
    string amr_report_blastp;	
 		string amr_report_blastx;
 		bool hmmChunks = false;
 		bool tblastnChunks = false;
	  const string annotS (" -gfftype " + Gff::names [(size_t) gffType] + ifS (lcl, " -lcl"));
    {
      //                               target ref    
  		#define BLAST_FMT    "-outfmt '6 qseqid sseqid qstart qend qlen sstart send slen qseq sseq'"
  		#define TBLASTN_FMT  "-outfmt '6 sseqid qseqid sstart send slen qstart qend qlen sseq qseq'"

  		
  		// PD-2967
  		const string blastp_par ("-comp_based_stats 0  -evalue 1e-10  -seg no  -max_target_seqs 10000");  // PAR
  		  // was: -culling_limit 20  // PD-2967
  		if (! emptyArg (prot))
  		{
  			string gff_prot_match;
 			  string gff_dna_match;
  			if (getFileSize (unQuote (prot_flat)))
  			{
    			findProg ("blastp");  			
    			findProg ("hmmsearch");
    			
    			string prot1 (prot_flat);  // Protein FASTA with no dashes in the sequences
          size_t nProt = 0;
          size_t protLen_max = 0;
          size_t protLen_total = 0;

          try
          {
      	    fastaCheck (prot_flat, true, qcS, logFName, nProt, protLen_max, protLen_total, noString);
      	  }
      	  catch (...)
      	  {
       	    bool fixable = false;
        	  {
        	    LineInput f (logFName);
        	    while (f. nextLine ())
        	      // Cf. fasta_check.cpp
          	    if (contains (f. line, "Hyphen in the sequence"))  
                {
            	    const Warning warning (stderr);
            		  stderr << "Ignoring dash '-' characters in the sequences of the protein file " << prot;
            		  fixable = true;
            		  break;
            		}
          	    else if (contains (f. line, "Too many ambiguities"))  
                {
            	    const Warning warning (stderr);
            		  stderr << "Removing sequences with >= " << ambigS << " Xs from the protein file " << prot;
            		  fixable = true;
            		  break;
            		}
              #if 0
          	    else if (contains (f. line, "'*' at the sequence end"))  
                {
            	    const Warning warning (stderr);
            		  stderr << "Removing '*' from the ends of protein sequences in " << prot;
            		  fixable = true;
            		  break;
            		}
            	#endif
          	}
          	if (fixable)
          	  removeFile (logFName);
          	else
          	  throw;
            prot1 = shellQuote (tmp + "/prot");
          	fastaCheck (prot_flat, true, qcS, logFName, nProt, protLen_max, protLen_total, prot1);
      	  }
    			
   			  // gff_check
    			if (! emptyArg (gff) && ! contains (parm, "-bed"))
    			{
    			  {
      			  bool gffProtMatchP = false;
      			  switch (gffType)
      			  {
      			    case Gff::genbank:
          			  {
            			  LineInput f (unQuote (prot1));
            			  while (f. nextLine ())
            			    if (   ! f. line. empty () 
            			        && f. line [0] == '>'
            			       )
            			    {
            			      gffProtMatchP = contains (f. line, "[locus_tag=");  
            			      break;
            			    }
            			}
            			break;
            	  case Gff::microscope: gffProtMatchP = true; break;
            	  case Gff::prodigal:   gffProtMatchP = true; break;
            	  default: break;
            	}
      			  if (gffProtMatchP)
      			    gff_prot_match = " -gff_prot_match " + tmp + "/prot_match";
      			}
    			  prog2dir ["gff_check"] = execDir;		
    			  string dnaPar;
    			  if (! emptyArg (dna))
    			  {
    			    dnaPar = " -dna " + dna_flat;
    			    if (gffType == Gff::pseudomonasdb)
    			      gff_dna_match = " -gff_dna_match " + tmp + "/dna_match";
    			  }
    			  try 
    			  {
    			    exec (fullProg ("gff_check") + gff_flat + annotS + " -prot " + prot1 + dnaPar + gff_prot_match + gff_dna_match + qcS + " -log " + logFName, logFName);
    			  }
    			  catch (...)
    			  {
    			    StringVector vec (logFName, (size_t) 10, false);  // PAR
    			    if (! vec. empty ())
    			      if (vec [0]. empty ())
    			        vec. eraseAt (0);
    			    if (! vec. empty ())
    			      if (vec [0] == error_caption)
    			        vec. eraseAt (0);
    			    throw runtime_error ("GFF file mismatch.\n" + vec. toString ("\n"));  // PD-3289, PD-3345
    			  } 
    			}
    			    			
    			stderr. section ("Running blastp");
    			{
      			const Chronometer_OnePass_cerr cop ("blastp");
      			// " -task blastp-fast -word_size 6  -threshold 21 "  // PD-2303
      			exec (fullProg ("blastp") + " -query " + prot1 + " -db " + tmp + "/db/AMRProt.fa" + "  " 
      			      + blastp_par + " -task blastp-fast"  // "-threshold 100 -window_size 15" are faster, but may miss hits, see SB-3643
      			      + getBlastThreadsParam ("blastp", min (nProt, protLen_total / 10000)) + " " BLAST_FMT " -out " + tmp + "/blastp > /dev/null 2> " + tmp + "/blastp-err", tmp + "/blastp-err");
      		}
    			  
    			stderr. section ("Running hmmsearch");
    			{
       			const Chronometer_OnePass_cerr cop ("hmmsearch");
      			ASSERT (threads_max >= 1);
      			if (threads_max > 1 && nProt > threads_max / 2)  // PAR
      			{
        		  createDirectory (tmp + "/hmm_chunk");
        		  exec (fullProg ("fasta2parts") + prot1 + " " + to_string (threads_max) + " " + tmp + "/hmm_chunk" + qcS + " -log " + logFName, logFName);
        		  createDirectory (tmp + "/hmmsearch_dir");
        		  createDirectory (tmp + "/dom_dir");
              Threads th (threads_max - 1, true);  
        		  DirItemGenerator dig (0, tmp + "/hmm_chunk", false);
        		  string item;
        		  while (dig. next (item))
          			th. exec (fullProg ("hmmsearch") 
      			              + "  --tblout "    + tmp + "/hmmsearch_dir/" + item + "  --noali"
      			              + "  --domtblout " + tmp + "/dom_dir/"       + item + "  --cut_tc  -Z 10000  --cpu 0  " + tmp + "/db/AMR.LIB" + " " + tmp + "/hmm_chunk/" + item + " > /dev/null 2> /dev/null"
      			             );
        		  hmmChunks = true;
        	  }
        	  else
        			exec (fullProg ("hmmsearch") + " --tblout " + tmp + "/hmmsearch  --noali  --domtblout " + tmp + "/dom  --cut_tc  -Z 10000  --cpu " + to_string (threads_max - 1) + "  " + tmp + "/db/AMR.LIB" + " " + prot1 + " > /dev/null 2> /dev/null", logFName);
        	}
  		  }
  		  else
  		  {
  		    OFStream::create (tmp + "/blastp");
  		    OFStream::create (tmp + "/hmmsearch");
  		    OFStream::create (tmp + "/dom");
  		  }  

  		  amr_report_blastp = "-blastp " + tmp + "/blastp  -hmmsearch " + tmp + "/hmmsearch  -hmmdom " + tmp + "/dom";
  			if (! emptyArg (gff))
  			  amr_report_blastp += "  -gff " + gff_flat + gff_prot_match + gff_dna_match + annotS;
  		}  		

  		
  		if (! emptyArg (dna))
  		{
  		  if (getFileSize (unQuote (dna_flat)))
    		{
          size_t nDna = 0;
          size_t dnaLen_max = 0;
          size_t dnaLen_total = 0;
          fastaCheck (dna_flat, false, qcS, logFName, nDna, dnaLen_max, dnaLen_total, noString); 
          const string blastx (dnaLen_max > 100000 ? "tblastn" : "blastx");  // PAR  // SB-3643

    			stderr. section ("Running " + blastx);
    			findProg (blastx);
          {
       			const Chronometer_OnePass_cerr cop (blastx);
            const string tblastn_par (blastp_par + "  -task tblastn-fast  -threshold 100  -window_size 15  -db_gencode " + to_string (gencode));  // SB-3643, PD-4522
        		const string blastx_par  (blastp_par + "  -word_size 3  -query_gencode " + to_string (gencode));
      			ASSERT (threads_max >= 1);
      			if (blastx == "blastx")
        			exec (fullProg ("blastx") + "  -query " + dna_flat + " -db " + tmp + "/db/AMRProt.fa" + "  "
            			  + blastx_par + " " BLAST_FMT " " + getBlastThreadsParam ("blastx", min (nDna, dnaLen_total / 10002))
            			  + " -out " + tmp + "/blastx > /dev/null 2> " + tmp + "/blastx-err", tmp + "/blastx-err");
            else
            {
              ASSERT (blastx == "tblastn");
        			findProg ("makeblastdb");
           	  exec (fullProg ("makeblastdb") + " -in " + dna_flat + " -out " + tmp + "/nucl" + "  -dbtype nucl  -logfile " + tmp + "/makeblastdb.log", tmp + "/makeblastdb.log");  
        			if (threads_max > 1)
        			{
          		  createDirectory (tmp + "/AMRProt_chunk");
          		  exec (fullProg ("fasta2parts") + " " + shellQuote (db + "/AMRProt.fa") + " " + to_string (threads_max) + " " + tmp + "/AMRProt_chunk" + qcS + " -log " + logFName, logFName);
          		  createDirectory (tmp + "/tblastn_dir");
          		  createDirectory (tmp + "/tblastn_dir.err");
                Threads th (threads_max - 1, true);  
          		  DirItemGenerator dig (0, tmp + "/AMRProt_chunk", false);
          		  string item;
          		  while (dig. next (item))
            			th. exec (fullProg ("tblastn") + "  -db " + tmp + "/nucl  -query " + tmp + "/AMRProt_chunk/" + item + "  "
              			        + tblastn_par + " " TBLASTN_FMT "  -out " + tmp + "/tblastn_dir/" + item + " > /dev/null 2> " + tmp + "/tblastn_dir.err/" + item);
          		  tblastnChunks = true;
          	  }
          	  else
          			exec (fullProg ("tblastn") + "  -db " + tmp + "/nucl  -query " + tmp + "/db/AMRProt.fa  "
              			  + tblastn_par + " " TBLASTN_FMT "  -out " + tmp + "/blastx > /dev/null 2> " + tmp + "/tblastn-err", tmp + "/tblastn-err");
            }
          }

          if (blastn)
      		{
      			findProg ("blastn");
      			stderr. section ("Running blastn");
       			const Chronometer_OnePass_cerr cop ("blastn");
      			exec (fullProg ("blastn") + " -query " + dna_flat + " -db " + tmp + "/db/AMR_DNA-" + organism1 + ".fa  -evalue 1e-20  -dust no  -max_target_seqs 10000  " 
      			      + getBlastThreadsParam ("blastn", min (nDna, dnaLen_total / 2500000)) + " " BLAST_FMT " -out " + tmp + "/blastn > " + logFName + " 2> " + tmp + "/blastn-err", tmp + "/blastn-err");
      		}
    		}
    		else
    		{
  		    OFStream::create (tmp + "/blastx");
  		    OFStream::create (tmp + "/len");
      		if (blastn)
    		    OFStream::create (tmp + "/blastn");
  		  }
   		  amr_report_blastx = "-blastx " + tmp + "/blastx  -dna_len " + tmp + "/len";
  		}
  	}
  	
  	
  	if (hmmChunks)
  	{
  	  concatTextDir (tmp + "/hmmsearch_dir", tmp + "/hmmsearch");
  	  concatTextDir (tmp + "/dom_dir",       tmp + "/dom");
  	}

  	if (tblastnChunks)
  	{
  	  concatTextDir (tmp + "/tblastn_dir",     tmp + "/blastx");
  	//concatTextDir (tmp + "/tblastn_dir.err", tmp + "/tblastn-err");
  	}


  	if (suppress_common)
  	{
			OFStream outF (tmp + "/suppress_prot");
			LineInput f (db + "/AMRProt-suppress.tsv");
			while (f. nextLine ())
			  if (! isLeft (f. line, "#"))
  			{
  			  string org, accver;
  			  istringstream iss (f. line);
  			  iss >> org >> accver;
  			  QC_ASSERT (! accver. empty ());
  			  if (org == organism1)
  			    outF << accver << endl;
  			}
	  }
		
		
	  if (stxTyper)
	  {
  		stderr. section ("Running stxtyper");
 			const Chronometer_OnePass_cerr cop ("stxtyper");
			exec (  fullProg ("stxtyper") 
			      + "  -n " + dna_flat 
			      + prependS (blast_bin, "  --blast_bin ") 
			      + "  -o " + tmp + "/stxtyper"
			      + "  --name " + input_name 
			      + "  --amrfinder"
			      + ifS (print_node, "  --print_node")
			      + "  -q "  // ifS (getVerbosity () == -1, "  -q")
			      + ifS (qc_on, "  --debug")
			      + " > " + logFName
			      , logFName
			      );
	  }
		

    // tmp + "/amr", tmp + "/mutation_all"
		stderr. section ("Making report");
    const string printNode (print_node ? " -print_node" : "");
    const string nameS (" -name " + input_name);
    {
 			const Chronometer_OnePass_cerr cop ("amr_report");
      const string mutation_allS (mutation_all. empty () ? "" : ("-mutation_all " + tmp + "/mutation_all"));      
      const string coreS (add_plus ? "" : " -core");
      const string equidistantS (equidistant ? " -report_equidistant" : "");
  		exec (fullProg ("amr_report") + " -fam " + shellQuote (db + "/fam.tsv") + "  " + amr_report_blastp + "  " + amr_report_blastx
      		  + "  -organism " + strQuote (organism1) 
      		  + "  -mutation "    + shellQuote (db + "/AMRProt-mutation.tsv") 
      		  + "  -susceptible " + shellQuote (db + "/AMRProt-susceptible.tsv") 
      		  + " " + mutation_allS + " "
      		  + force_cds_report + " -pseudo" + coreS + equidistantS + printNode
      		  + (ident == -1 ? noString : "  -ident_min "    + toString (ident)) 
      		  + "  -coverage_min " + toString (cov)
      		  + ifS (suppress_common, " -suppress_prot " + tmp + "/suppress_prot")  
      		  + nameS + qcS + " " + parm + " -log " + logFName + " > " + tmp + "/amr", logFName);
  	}
		if (blastn)
		{
 			const Chronometer_OnePass_cerr cop ("dna_mutation");
      const string mutation_allS (mutation_all. empty () ? "" : ("-mutation_all " + tmp + "/mutation_all.dna")); 
			exec (fullProg ("dna_mutation") + tmp + "/blastn " + shellQuote (db + "/AMR_DNA-" + organism1 + ".tsv") + " " + strQuote (organism1) + " " + mutation_allS 
			      + nameS + printNode + qcS + " -log " + logFName + " > " + tmp + "/amr-snp", logFName);
	    {
  			ofstream f (tmp + "/amr", ios_base::out | ios_base::app);
  			copyText (tmp + "/amr-snp", 1, f);
  	  }
      if (! mutation_all. empty ())
      {
  			ofstream f (tmp + "/mutation_all", ios_base::out | ios_base::app);
  			copyText (tmp + "/mutation_all.dna", 1, f);
  	  }
	  }
	  if (stxTyper)
    {
			ofstream f (tmp + "/amr", ios_base::out | ios_base::app);
			copyText (tmp + "/stxtyper", 1, f);
	  }


    // Column names are from amr_report.cpp

    // Sorting AMR report
    // PD-2244, PD-3230
    {
      StringVector amrSortColumns;
      if (! (emptyArg (dna) && emptyArg (gff)))
        amrSortColumns << contig_colName << start_colName << stop_colName << strand_colName;    
      amrSortColumns << prot_colName << genesymbol_colName;  

      {
        TextTable amrTab (tmp + "/amr");
        amrTab. sort (amrSortColumns);
        amrTab. rows. uniq ();  // PD-4297
        amrTab. qc ();
        Cout out (output);
   		  amrTab. saveText (*out);
      }

      // Sorting mutation_all
      if (! mutation_all. empty ())
      {
        TextTable mutation_allTab (tmp + "/mutation_all");
        mutation_allTab. sort (amrSortColumns);
        mutation_allTab. rows. uniq ();
        mutation_allTab. qc ();
        mutation_allTab. saveFile (mutation_all);
        if (qc_on)
        {
          const TextTable::ColNum i = mutation_allTab. col2num (subtype_colName);
          for (const StringVector& row : mutation_allTab. rows)
            QC_ASSERT (row [i] == "POINT");
        }
      }
    }


    if (! emptyArg (prot_out))
    {
      prepare_fasta_extract (StringVector {prot_colName, genesymbol_colName, elemName_colName}, "prot_out", false);
      exec (fullProg ("fasta_extract") + prot + " " + tmp + "/prot_out -aa" + qcS + " -log " + logFName + " > " + prot_out, logFName);  
    }
    if (! emptyArg (dna_out))
    {
      prepare_fasta_extract (StringVector {contig_colName, start_colName, stop_colName, strand_colName, genesymbol_colName, elemName_colName}, "dna_out", false);
      exec (fullProg ("fasta_extract") + dna_flat + " " + tmp + "/dna_out" + qcS + " -log " + logFName + " > " + dna_out, logFName);  
    }
    if (! emptyArg (dnaFlank5_out))
    {
      prepare_fasta_extract (StringVector {contig_colName, start_colName, stop_colName, strand_colName, genesymbol_colName, elemName_colName}, "dna_out", true);
      //                                   0               1              2             3
      TextTable t (tmp + "/dna_out");
      t. qc ();
      for (StringVector& row : t. rows)
        if (row [3] == "+")
          row [1] = to_string (max (1, stoi (row [1]) - (int) dnaFlank5_size));
        else
          row [2] = to_string (stoi (row [2]) + (int) dnaFlank5_size);
      t. saveHeader = false;
      t. qc ();
      t. saveFile (tmp + "/dnaFlank5_out");
      exec (fullProg ("fasta_extract") + dna_flat + " " + tmp + "/dnaFlank5_out" + qcS + " -log " + logFName + " > " + dnaFlank5_out, logFName);  
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



