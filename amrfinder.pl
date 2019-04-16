#!/usr/bin/env perl 
#
# PUBLIC DOMAIN NOTICE This software/database is "United States Government
# Work" under the terms of the United States Copyright Act. It was written as
# part of the authors' official duties for the United States Government and
# thus cannot be copyrighted. This software/database is freely available to the
# public for use without a copyright notice. Restrictions cannot be placed on
# its present or future use.
# 
# Although all reasonable efforts have been taken to ensure the accuracy and
# reliability of the software and data, the National Center for Biotechnology
# Information (NCBI) and the U.S. Government do not and cannot warrant the
# performance or results that may be obtained by using this software or data.
# NCBI, NLM, and the U.S. Government disclaim all warranties as to performance,
# merchantability or fitness for any particular purpose.
# 
# In any work or product derived from this material, proper attribution of the
# authors as the source of the software or data should be made, using:
# https://ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/ as the
# citation.
# 
use strict;
use warnings;
#use diagnostics;
use Getopt::Long;
use File::Temp qw( tempdir );
use File::Basename;
use Net::FTP;
use Cwd qw( getcwd abs_path );
my $curr_version = '$Revision: 39138 $';
our $DEBUG = 0;

# todo:
# - make --version print database version as well
# - Add interpretation of error messages and quick-abort if
#   a forked process returns with an error
# - Add option to do single-threaded
# - Add test for availability of -num_threads option to blastp and blastx
#   - enable -num_threads 6 for blastp and blastx if it is available

# - rework to add functionality from amrfinder.cpp
# - Add compiling DNA and other blast databases as necessary
# - Add running blastn for point-mutations

my $usage = <<END;

amrfinder [-h|-v]
amrfinder [-options] -p <Protein FASTA>
amrfinder [-options] -n <genomic sequence FASTA>

    -h print help text
    -v print version information 
    -U update AMRFinder database
Options used for either of the operating modes:
    -d <dir> Directory containing the AMR database
    -o <file.tsv> tabfile output to this file instead of STDOUT
    -q quiet mode (don't print status messages to STDERR)
Options relating to protein input (-p):
    -p <protein> Amino-acid sequences to search using BLASTP and HMMER
    -g <gff> GFF file indicating genomic location for proteins in -p <protein>
Options relating to nucleotide sequence input (-n)
    -n <nucleotide> genomic sequence to search using BLASTX
    -i <0.9> Minimum proportion identical translated AA residues
    -c <0.5> Minimum coverage of reference protein sequence
    -t <11> Translation table for blastx

END

my $usage_long = <<END;

amrfinder [-h|-v]
amrfinder [-options] -p <Protein FASTA>
amrfinder [-options] -n <genomic sequence FASTA>


    -h/--help print help text
    -v/--version print version information 
    -U/--update_data update AMRFinder database by downloading latest version
        from the NCBI ftp site
Options used for either of the operating modes:
    -d/--database <dir> Directory containing the AMR database
    -o/--output <file.tsv> tabfile output to this file instead of STDOUT
    -q/--quiet don't print status messages to STDERR
Options relating to protein input (-p):
    -p/--protein <protein> Amino-acid sequences to search using BLASTP and 
        HMMER
    -g <gff> GFF file indicating genomic location for proteins in -p <protein>
Options relating to nucleotide sequence input (-n)
    -n <nucleotide> genomic sequence to search using blastx
    -i <0.9> Minimum proportion identical translated AA residues 
    -c <0.5> Minimum coverage of reference protein sequence for
        a "BLASTX" match vs. a "PARTIALX" match
    -t <11> Translation table for blastx, for options and their meaning see:
        https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
END


if (@ARGV == 0) { print $usage; exit 1; }
my $command_line = "$0 @ARGV";

# Some paths
my $database_dir    = $ENV{AMRFINDER_DB} || '';
my $bin             = ''; # set to path to AMRFinder executables
our $FAMREPORT = "${bin}amr_report";
our $FASTA_CHECK = "${bin}fasta_check";
our $GFF_CHECK = "${bin}gff_check";
our $BLASTP = 'blastp';
our $BLASTX = 'blastx';
our $HMMER  = 'hmmsearch';
our $tempdir;
if ($DEBUG) {
    $tempdir = tempdir( CLEANUP => 0 );
    print STDERR "Temporary files in $tempdir\n";
} else {
    $tempdir = tempdir( CLEANUP => 1 );
}

# Options
my ($outfile, $out_fa, $help, $version); # opts with no defaults
my ($prot_file, $nuc_file, $quiet, $update_data) = ('', '', 0, 0);
#my $min_ident       = 0.9; 
my $min_ident       = undef; # default now set in opts_ok
#my $min_cov         = 0.5;
my $min_cov         = undef; # default now set in opts_ok
#my $trans_table     = 11;
my $trans_table     = undef; # default now set in opts_ok
my $gff             = '';
my $NUM_THREADS     = 4; # parameter to -num_threads or --cpu 
GetOptions(
    'protein=s'             => \$prot_file,
    'nucleotide=s'          => \$nuc_file,
    'database=s'            => \$database_dir,
    'output=s'              => \$outfile,
    'fasta_out=s'           => \$out_fa,
    'gff=s'                 => \$gff,
    'ident_min=f'           => \$min_ident,
    'coverage_min=f'        => \$min_cov,
    'translation_table=i'   => \$trans_table,
    'U|update_data'         => \$update_data,
    'quiet'                 => \$quiet,
    'help|?'                => \$help,
    'version'               => \$version,
    #    'D'                     => \$DEBUG,
);

print $usage_long and exit() if ($help);
print "version $curr_version\n" if ($version);
if ($DEBUG) { $quiet = 0; } # turn quiet mode off if we're in debug mode 
print STDERR "Running $command_line\n" unless ($quiet);
print STDERR "$usage\nERROR: options not understood: @ARGV\n" if (@ARGV);
opts_ok(); # sanity check on options and files before we get started
paths_ok(); # check for exectuables and database files, this may modify globals!
if ($update_data) {
    update_database($database_dir);
    make_databases($database_dir) unless (-e "$database_dir/AMRProt.phr");
    exit 0 unless ($prot_file or $nuc_file); # exit if there's nothing left to do
}
make_databases($database_dir) unless (-e "$database_dir/AMRProt.phr");
my $abs_database_dir = abs_path($database_dir);

# end of config and error checking, now run AMRFinder

my $report;
if ($prot_file) {
    check_prot();
    print STDERR "Running AMRFinder on $prot_file with database $abs_database_dir\n" unless ($quiet);
    $report = run_prot(); 
} else {
    check_nuc();
    $min_ident   ||= 0.9;
    $min_cov     ||= 0.9;
    $trans_table ||= 11;
    print STDERR "Running AMRFinder on $nuc_file with database $abs_database_dir\n" unless ($quiet);
    $report = run_nuc();
}

if ($outfile) { 
    open (OUT, '>', $outfile) or die "Couldn't create $outfile: $!";
    print OUT $report;
} else {
    print $report;
}

### End MAIN ###

##############################################################################
# update_database - download a new AMRFinder database from FTP site if available
# Will update a database in in a subdirectory of AMRFinder run directory or 
# under $AMRFINDER_DB
sub update_database {
    my $database_dir = shift;
    my $amrfinder_dir = dirname(abs_path(__FILE__));
    if ($database_dir !~ m#/latest$# or (-e $database_dir and ! -l "$database_dir/latest")) {
        warn "Updating database directory only works for databases with the default data\n",
            "directory format. Please see https://github.com/ncbi/amr/wiki for details\n",
            "AMRFinder will create new database directories as subdirectories\n",
            "Current database directory is $database_dir\n";
    }
    if ($database_dir !~ m#/latest$#) {
        $database_dir = "$database_dir/latest";
    }
    download_database($database_dir);
    #    if ( get_latest_db_version() gt get_database_version($database_dir)) {
}

##############################################################################
# get_database_version - just figure out what the latest link points to
# $database_dir is the current database directory
# Only returns a version if the database_dir ends in /latest and is a symlink,
# otherwise returns ''
sub get_database_version {
    my $database_dir = shift;
    if ($database_dir =~ /\/latest/ and -l $database_dir) {
        my $dest = readlink($database_dir);
        my $version = basename($dest);
        return $version;
    } else {
        return '';
    }
}


##############################################################################
# download_database
# INPUTS: $database_dir to put the database in, should be of format
#  .../latest where latest is a soft link
sub download_database {
    my $database_dir = shift;
    $database_dir =~ s#((?<!\\)/)+$##; # remove tailng slashes not preceeded by a '\'
    my $base_dir = $database_dir;
    $base_dir =~ s#/latest$##;
    #my ($base_dir) = $database_dir =~ m#^(.*)/latest#;
    my $cwd = getcwd();
    if (! -d "$base_dir") {
        system("mkdir -p $base_dir") and die "Couldn't mkdir $base_dir";
    }
    my $latest_version = get_latest_db_version();
    print STDERR "Dowloading AMRFinder database version $latest_version\n" unless ($quiet);
    my $ftp = Net::FTP->new("ftp.ncbi.nlm.nih.gov", Passive => 1) or die "Cannot connect to hostname $@" ;
    $ftp->login("anonymous", 'AMRFinder@ncbi') or die "Couldn't log in to ftp.ncbi.nlm.nih.gov: ", $ftp->message;
    $ftp->cwd('/pathogen/Antimicrobial_resistance/AMRFinder/data/latest')
        or die "ERROR on FTP site, latest not found, please report to pd-help\@ncbi.nlm.nih.gov ", $ftp->message;
    chdir $base_dir or die "Couldn't chdir to $base_dir: $!";
    print STDERR "Downloading $latest_version into $base_dir/$latest_version\n" if ($DEBUG);
    if (-e $latest_version) {
        warn "$latest_version already exists, overwriting what was there\n";
    } else {
        mkdir $latest_version or die "Couldn't mkdir $latest_version: $!";
    }
    chdir $latest_version or die "Couldn't chdir to $base_dir/$latest_version: $!";
    my @files = $ftp->ls() or die "Couldn't use FTP to ls";
    $ftp->binary() or die "Couldn't switch to binary mode: ", $ftp->message;
    foreach my $file (@files) {
        print STDERR "Fetching $file\n" unless ($quiet);
        $ftp->get($file) or die "Couldn't fetch $file: ", $ftp->message;
    }
    chdir '..' or die "Couldn't chdir to ..: $!";
    if (-e 'latest') {
        unlink 'latest' or die "Couldn't remove link latest from $base_dir: $!" ;
    }
    symlink $latest_version, 'latest' 
        or die "Couldn't create link from $base_dir/latest to $base_dir/$latest_version: $!";
    chdir $cwd or die "Couldn't chdir back to $cwd: $!";
    if ($database_dir !~ m#/latest$#) {
        $database_dir = "$database_dir/latest";
    }
    return 1;
}

    
##############################################################################
# get_latest_db_version - find the latest version of AMRFinder
# OUTPUTS: 
#   connects and sets global $ftp to Net::FTP connection 
#   Returns directory name pointed to by 
# ftp://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinder/data/latest
sub get_latest_db_version {
    my $ftp = Net::FTP->new("ftp.ncbi.nlm.nih.gov", Passive => 1) or die "Cannot connect to hostname $@" ;
    $ftp->login("anonymous", 'AMRFinder@ncbi') or die "Couldn't login to $ftp", $ftp->message();
    $ftp->pasv() or die "Couldn't enter passive mode: ", $ftp->message();
    $ftp->cwd('/pathogen/Antimicrobial_resistance/AMRFinder/data/') 
        or die "Couldn't cd to ftp.ncbi.nlm.nih.gov:/pathogen/Antimicrobial_resistance/AMRFinder/data/: ", $ftp->message();
    my @dirs = $ftp->dir();
    if (@dirs == 0) {
        die "Error listing databases on ftp.ncbi.nlm.nih.gov:/pathogen/Antimicrobial_resistance/AMRFinder/data/\n", $ftp->message();
    }
    my ($latest) = grep index($_, ' latest -> ') != -1, @dirs;
    return '' unless (defined $latest);
    ($latest) = (split(/\s/, $latest))[-1];
    return($latest);
}

##############################################################################
# check_nuc - check nucleotide FASTA file
# reads globals: $nuc_file
# Runs: $FASTA_CHECK
sub check_nuc {
    my @cmd = ($FASTA_CHECK, $nuc_file, "-hyphen");
    system(@cmd) == 0 or
        die "\n"; # fasta_check should print the error message
}

##############################################################################
# check_prot - check protein FASTA file (possibly and GFF)
# reads globals: 
#   $prot_file
#   $gff
# Runs $FASTA_CHECK and possibly $GFF_CHECK
sub check_prot {
    my @cmd = ($FASTA_CHECK, $prot_file, '-aa', '-hyphen');
    system(@cmd) == 0 or
        die "\n"; # fasta_check should print the error message
    if ($gff) {
        # check for locus_tag
        my $locus_tag_opt = '';
        my @cmd = ($GFF_CHECK, $gff, '-fasta', $prot_file);
        open(my $fh, '<', $prot_file)
            or die "Error opening -p $prot_file: $!\n";
        while(<$fh>) {
            if (/^>/) {
                if (/\[locus_tag=/) {
                    push @cmd, '-locus_tag';
                }
                last; # just check the first defline, skip the rest
            }
        }
        close($fh);
        system(@cmd) == 0 or
            die "\n"; # gff_check should print the error message
    }
}

##############################################################################
# run AMRFinder-dna
# reads globals:
# $BLASTX
# $tempdir
# $database_dir
# $trans_table
# $FAMREPORT
# $min_ident
# $min_cov
# $nuc_file
sub run_nuc {
    my $cmd = "$BLASTX -db $database_dir/AMRProt -query $nuc_file -show_gis -word_size 3 -evalue 1e-20 -query_gencode $trans_table "
    . "-seg no -comp_based_stats 0 -max_target_seqs 10000 -outfmt '6 qseqid sseqid length nident qstart qend qlen sstart send slen qseq' "
    . " > $tempdir/blastx 2> $tempdir/blastx.err";
    print STDERR "Running blastx...\n" unless ($quiet);
    system($cmd) == 0 or die "ERROR running '$cmd'";
    if (-s "$tempdir/blastx.err") {
        die "ERROR running $cmd: " . `cat $tempdir/blastx.err`;
    } 
    my $fh;
    $cmd = "$FAMREPORT -fam $database_dir/fam.tab -blastx $tempdir/blastx -ident_min $min_ident -complete_cover_min $min_cov -pseudo";
    print STDERR "Running $cmd\n" if ($DEBUG);
    open ($fh, "-|", $cmd) or die "Error running $cmd: $!";
    my $report = join("", <$fh>);
    return $report;
}


##############################################################################
# run AMRFinder-prot
# reads globals:
# $BLASTP
# $tempdir
# $datbase_dir
# $HMMER
# $gff
# $FAMREPORT
# returns string containing AMR Report
sub run_prot {
    # first blast
    my $cmd = "$BLASTP -task blastp-fast -db $database_dir/AMRProt -query $prot_file -show_gis -word_size 6 -threshold 21 -evalue 1e-20 -comp_based_stats 0 -outfmt '6 qseqid sseqid length nident qstart qend qlen sstart send slen qseq' > $tempdir/blastp 2> $tempdir/blastp.err";
    print STDERR "Running blastp...\n" unless ($quiet);
    my @waitfor = run_forked_system($cmd);
    $cmd = "$HMMER --tblout $tempdir/hmmsearch --noali --domtblout $tempdir/dom --cut_tc -Z 10000 --cpu 6 $database_dir/AMR.LIB $prot_file > $tempdir/hmmer.out 2> $tempdir/hmmer.err";
    print STDERR "Running hmmsearch...\n" unless ($quiet);
    push @waitfor, run_forked_system($cmd);
    while (wait() != -1) {}; # wait for the processes to finish
    if (-s "$tempdir/blastp.err") {
        die "Error returned from BLAST: " . `cat $tempdir/blastp.err`; # ugly
    } elsif (-s "$tempdir/hmmer.err") {
        die "Error returned from hmmsearch: " . `cat $tempdir/hmmer.err`; # ugly
    }
    $gff = "-gff $gff" if ($gff);
    $cmd = "$FAMREPORT -fam $database_dir/fam.tab -blastp $tempdir/blastp $gff -hmmsearch $tempdir/hmmsearch -hmmdom $tempdir/dom -out $tempdir/sseqid ";
    my $fh;
    print STDERR "Running $cmd\n" if ($DEBUG);
    open($fh, '-|', $cmd) or die "Couldn't run $FAMREPORT: $!";
    my $report = join("", <$fh>);
    return($report);
}

##############################################################################
# Used to run HMMER and blastp simultaneously
# Inputs: The shell command to run
# Outputs: Returns the PID of the child process that was forked
sub run_forked_system {
    my $cmd = shift;
    my $child = fork();
    print STDERR "forked and running $cmd\n" if ($DEBUG);
    if (! defined($child)) {
        die "Something went wrong trying to fork and exec $cmd";
    } elsif ($child == 0) { # we're in the child process
        exec($cmd) or die;
    }
    return $child;
}


##############################################################################
# check the options for sanity before running anything
# Reads everything from globals set at beginning of script
sub opts_ok {
    my @errors;
    my @warnings;
    if (! $nuc_file && ! $prot_file && ! $update_data) {
        push @errors, "ERROR: At least one of -p <file> or -n <file> or -U is required\n";
    } elsif ($nuc_file && $prot_file) {
        push @errors, "ERROR: Currently AMRFinder cannot search both protein and nucleotide sequences at the same time\n";
    } elsif ($nuc_file && ! -e $nuc_file) {
        push @errors, "ERROR: $nuc_file not found\n";
    } elsif ($prot_file && ! -e $prot_file) {
        push @errors, "ERROR: $prot_file not found\n";
    }

    if ($gff && $nuc_file) {
        warn "WARNING: -gff $gff option is not used when searching nucleotide sequences\n";
    } 
    if ($out_fa && $nuc_file) {
        push @errors, "ERROR: FASTA output is not implemented for nucleotide searches\n";
    }
    if ($min_ident && $prot_file) {
        push @errors, "ERROR: -i $min_ident option is not implemented for protein searches\n";
    }
    if ($min_cov && $prot_file) {
        push @errors, "ERROR: -c $min_cov option is not implemented for protein searches\n";
    }
    if ($trans_table && $prot_file) {
        warn "WARNING: Translation table not used for protein searches\n";
    }
    # now check if the database exists and is available
    #    if (! -e "$database_dir/AMRProt" or ! -e "$database_dir/fam.tab" or ! -e "$database_dir/AMR.LIB") {
    #    push @errors, "ERROR: Couldn't find a complete AMRFinder database in $database_dir";
    #}
    if (@errors) {
        print $usage;
        print @errors;
        exit 1;
    }
    return 1;
}

##############################################################################
# paths_ok - look for exectuables and paths to things, printing error
# messages if they can't be found
# INPUTS: (Globals)
#   $0
#   $HMMER
#   $BLASTP
#   $BLASTX
#   $database_dir
#   $FAMREPORT
# OUTPUTS: (May alter globals)
#   $database_dir
#   $FAMREPORT
#   $FASTA_CHECK
#   $GFF_CHECK
sub paths_ok {
    my @errors;
    # first check the AMRFinder database
    my $amrfinder_dir = dirname(abs_path(__FILE__));
    if ($database_dir) {
        if (!check_db_dir($database_dir)) {
            my $msg = "Database directory $database_dir does not appear to be in the correct format\n";
            if ($update_data) {
                warn($msg);
            } else {
                die($msg);
            }
        }
    } else {
        if (check_db_dir("$amrfinder_dir/data/latest")) {
            $database_dir = "$amrfinder_dir/data/latest";
        } elsif (check_db_dir("../data")) { # for Slava
            $database_dir = "../data";
        } else {
            if ($database_dir eq '' and $update_data) { 
                # it's ok if we're going to download data first
                warn "No current AMRFinder database.\n";
                $database_dir = "$amrfinder_dir/data/latest";
            } else {
                push @errors, "ERROR: AMRFinder database files not found\n" .
                              "       You can use the -U option to automatically download the latest data\n" .
                              "       and install it in the default location \n" .
                              "       $amrfinder_dir/data";
            }
        }
    }


    # Now check executables
    if ($prot_file) {
        $BLASTP = check_blast($BLASTP, $NUM_THREADS);
        if ($BLASTP =~ /^ERROR/) {
            push @errors, $BLASTP;
        }
        if (system("$HMMER -h > /dev/null")) {
            push @errors, "ERROR executing $HMMER, please make sure hmmsearch is installed and in your path";
        }
    } elsif ($nuc_file) {
        $BLASTX = check_blast($BLASTX, $NUM_THREADS);
        if ($BLASTX =~ /^ERROR/) {
            push @errors, $BLASTX;
        }
    }

    if (-e "$amrfinder_dir/$FAMREPORT" and !system("$amrfinder_dir/$FAMREPORT -help > /dev/null")) {  
        # found famreport
        $bin = "$amrfinder_dir";
        $FAMREPORT      = "${bin}/$FAMREPORT";
        $FASTA_CHECK    = "${bin}/$FASTA_CHECK";
        $GFF_CHECK      = "${bin}/$GFF_CHECK";
    } elsif (system("$FAMREPORT -help > /dev/null")) { # if it's not in our path
        push @errors, "ERROR finding required AMRFinder C++ executables, please make sure they\n" .
            "    are in your path or in the same directory as the amrfinder.pl executable";
    }

    if (@errors) {
        print $usage;
        print join("\n", @errors, '');
        exit 1;
    }
    return 1;
}

# check that blast is executable and return -num_threads option if you can 
# (Killing 2 birds with 1 stone)
# INPUTS: blast_program_name, $num_threads
# OUTPUTS: Returns either:
#   command-line to run (i.e., $blast_program_name -num_threads $num_threads)
#   or
#   an error message beginning with "ERROR executing"...
sub check_blast {
    my $progname = shift;
    my $num_threads = shift;
    open(my $fh, "-|", "$progname -h") 
        or return "ERROR executing $progname, please make sure $progname is installed and in your path";
    while(<$fh>) {
        if (/-num_threads \<Integer, \(.* and \=\<(\d+)/) {
	    $num_threads = $1 if ($1 < $num_threads);
            return "$progname -num_threads $num_threads";
        }
    }
    # if we got here then no num_threads option
    return "$progname";
}
    

##############################################################################
# make_databases - format blast and HMMer databases
# INPUT: The database directory
# OUTPUT: Runs makeblastdb
sub make_databases {
    my $dir = shift;
    my @cmd = ('makeblastdb', '-in', "$dir/AMRProt", '-dbtype', 'prot', '-out', "$dir/AMRProt");
    system(@cmd) == 0 
        or die "ERROR running @cmd: $?";
}


##############################################################################
# just make sure the minimum 3 required database files are present        
# Takes the database directory
# Returns 1 if database directory has required files
# Returns 0 if not
sub check_db_dir {
    my $dir = shift;
    if (-d "$dir" and 
        -f "$dir/AMRProt" and 
        -f "$dir/AMR.LIB" and
        -f "$dir/fam.tab") {
        return 1;
    } else {
        return 0;
    }
}
        
