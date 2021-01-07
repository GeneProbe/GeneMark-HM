#!/usr/bin/perl
# ===============================================
# Authors: Alex Lomsadze
#          Christophe Bonny
#          Francesco Strozzi
#          Mark Borodovsk
# Copyright: GeneProbe Inc.
#            Atlanta, GA, USA
# Copyright: Enterome,
#            94/96 avenue Ledru-Rollin 75011 Paris, France   
#
# Reference:  GeneMark-HM: Improving Gene Prediction in DNA Sequences of Human Microbiome
# under review
#
# Release: January 2021
#
# File: gm_hm.pl
#
# Project: this pipeline predicts protein coding genes 
#          in metagenomes assotiated with human host.
# ===============================================

use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use Cwd qw(abs_path cwd);
use File::Temp qw(tempfile tempdir);
use File::Spec;
use Data::Dumper;
use File::Path qw(remove_tree);
use Parallel::ForkManager;

# ------------------------------------------------
# general configuration
my $VERSION = "2.07";

my $verbose = 0;
my $debug = 0;

my $work_dir = cwd;
my %cfg = ();
# ------------------------------------------------
# command line options

# required
my $seq_file = '';
my $out_file = '';

# output optional
my $format = "gtf";
my $nt_file = '';
my $aa_file = '';

# algorithm optional
my $mgm2_only = '';
my $no_taxa = '';
my $no_gms2 = '';
my $only_gcode11 = '';

# data configuration
my $refdb = "$Bin/refdb/pan";
my $refmod = "$Bin/refmod";

# thresholds taxa
my $identity = 95;
my $evalue = 0.0001;
my $min_hit_size = 90;

# thresholds domain
my $ani_domain = 70;

# threshold gms2
my $genes_for_gms2 = 150;

# threshold gcode
my $min_gcode_contig = 500;

my $GCODE_CONST_A = 1.33;
my $GCODE_CONST_B = 20;

my $CONST_MIN_DIFF_LOGODD = 10;
my $CONST_MIN_LOGODD_FOR_NOT_11 = 20;

# other
my $cores = 1;
my $tmpf = '';
my $clean = '';

# ------------------------------------------------
# scripts and programs called by pipeline
my $parse_blast = "$Bin/parse_blast.pl";
my $probuild = "$Bin/probuild";
my $blastn = "blastn";

# ------------------------------------------------
my $command = '';

Usage() if @ARGV < 1;
ParseCMD();
CheckBeforeRun();

system( "date" ) if $verbose;

# $mgm is used to run program in gcode 11 mode only
# RunMGM() is used by default for gcode 11/4/etc
my $mgm = SetMGM($Bin);          $cfg{'set'}->{'mgm'} = $mgm;

my $gms = SetGMS($Bin);          $cfg{'set'}->{'gms'} = $gms;

my $hmm11 = SetHMM($Bin, "11");  $cfg{'set'}->{'hmm11'} = $hmm11;
my $hmm4  = SetHMM($Bin, "4");   $cfg{'set'}->{'hmm4'}  = $hmm4;
my $hmm25 = SetHMM($Bin, "25");  $cfg{'set'}->{'hmm25'} = $hmm25;
my $hmm15 = SetHMM($Bin, "15");  $cfg{'set'}->{'hmm15'} = $hmm15;

$tmpf = SetTmpFolder($tmpf);     $cfg{'set'}->{'tmpf'} = $tmpf;

if ($debug)
{
	print "cfg-set\n";
	print Dumper($cfg{'set'});
}

# ------------------------------------------------
# data processing starts here

# move to internal FASTA defline
my ($defline_status, $sanitized) = CheckFASTAdefline($seq_file);
my ($sfile, $trace_file) = FastaToInternalFormat( $seq_file, "seq.fasta", "seq.trace" );

$cfg{'workfasta'}->{'sfile'} = $sfile;
$cfg{'workfasta'}->{'trace_file'} = $trace_file;
$cfg{'workfasta'}->{'sanitized'} = $sanitized;
$cfg{'workfasta'}->{'defline_status'} = $defline_status;

if ($debug)
{
	print "cfg-workfasta\n";
	print Dumper($cfg{'workfasta'});
}

# run predictions
if ( $mgm2_only )
{
	print "running pipeline in MetaGeneMark-2 mode only ...\n" if $verbose;
	RunMGM( " --seq $sfile" );
}
else
{
	print "running initial predictions by MetaGeneMark-2 ...\n" if $verbose;

	chdir $tmpf;

	RunMGM( " --seq $sfile  --out mgm_ini.gtf  --format gtf  --NT mgm_ini.fna" );

#	$command  = "$mgm  --seq $sfile  --out mgm_ini.gtf  --format gtf  --NT mgm_ini.fna";
#	RunSystem( $command );

	if ( !$no_taxa )
	{
		print "running BlastN alignments ...\n" if $verbose;
		system( "date" ) if $verbose;

		$command  = "$blastn  -db $refdb  -evalue $evalue  -num_threads $cores  -strand plus ";
		$command .= "  -query mgm_ini.fna  -out blast.out  -outfmt 7  -max_target_seqs 20";
		RunSystem( $command );

		system( "date" ) if $verbose;

		$command  = "$parse_blast  --identity $identity  --in blast.out  --out taxa.set ";
		$command .= " --ver" if $verbose;
		RunSystem( $command );

		$command  = "$parse_blast  --identity $ani_domain  --in blast.out  --out taxa_all.set ";
		RunSystem( $command );
	}
 
	my $run_gms = 0;
	my $run_hmm = 0;

	$run_gms = FindCotigsWithMoreThanNgenes( "mgm_ini.gtf", "for_gms.tbl" ) if !$no_gms2;
	$run_hmm = FindCotigsWithHitsToRef( "taxa.set" ) if !$no_taxa;

	if ($debug)
	{
		print "# run GMS status: $run_gms\n";
		print "# run HMM status: $run_hmm\n";
	}

	if ( $run_gms == 0 and $run_hmm == 0 )
	{
		chdir $work_dir;
		RunMGM( " --seq $sfile" );
	}
	else
	{
		my %taxa_for_hmm;
		my %contigs_for_gms;

		# first run HMM on refhits if any

		if ( $run_hmm )
		{
			print "preparing for taxa GeneMark.hmm-2 runs ... \n" if $verbose;

			if ( ! -d "for_hmm" )
			{
				mkdir "for_hmm" or die "error on create folder for_hmm: $!\n";
			}
			chdir "for_hmm" or die "error on change folder: for_hmm\n";
			ParseTaxa( "../taxa.set" );
			SelectFastaByTaxa( $sfile,  "../taxa.set" );
			chdir $tmpf;

			my $manager = new Parallel::ForkManager( $cores );

			%taxa_for_hmm = GetTaxaNames("taxa.set");
			foreach my $key (keys %taxa_for_hmm)
			{
				$manager->start and next;
				RunHMM( "for_hmm/$key.seq", $key, "for_hmm" );
				$manager->finish;
			}

			$manager->wait_all_children;
		}

		# second run GMS2 on remaining long if any
		
		if ( $run_gms )
		{
			print "preparing for GeneMarkS-2 runs ... \n" if $verbose;
			if ( ! -d "for_gms" )
			{
				mkdir "for_gms" or die "error on create folder for_gms: $!\n";
			}
			chdir "for_gms" or die "error on change folder: for_gms\n";
			ParseContigs( "../for_gms.tbl", $sfile );
			chdir $tmpf;

			my %contig_to_taxa = ();
			my %contig_to_taxa_all = ();

			if ( $run_hmm )
			{
				%contig_to_taxa = LoadList( "taxa.set" );
				%contig_to_taxa_all = LoadList( "taxa_all.set" );
			}

			my %taxa_to_domain = LoadList( $refmod ."/domain.tbl" );

			%contigs_for_gms = GetContigNames("for_gms.tbl");
			foreach my $key (keys %contigs_for_gms)
			{
				if ( exists $contig_to_taxa{$key} )
				{
					print "# skipping GMS2 on contig with taxa hit; $key\n";
					next;
				}

				if ( ! -e $key )
				{
					mkdir $key or die "error on create folder $key: $!\n";
				}
				chdir $key or die "error on change folder: $key\n";

				if ((exists $contig_to_taxa_all{$key}) and (exists $taxa_to_domain{$contig_to_taxa_all{$key}} ))
				{
					if ( $taxa_to_domain{$contig_to_taxa_all{$key}} eq "bac" )
					{
						RunGMS( " --seq ../for_gms/$key  --genome-type bacteria");
					}
					elsif ( $taxa_to_domain{$contig_to_taxa_all{$key}} eq "arc" )
					{
						RunGMS( " --seq ../for_gms/$key  --genome-type archaea" );
					}
					else
						{ die "error, unexpected domain label found: $taxa_to_domain{$contig_to_taxa_all{$key}}\n"; }
				}
				else
				{
					RunGMS( " --seq ../for_gms/$key" );
				}
				
				chdir $tmpf;
			}
		}

		# third run MGM on remaining 

		chdir $work_dir;

		my $ignore_file = "mgm2_ignore";

		if ( $run_gms and $run_hmm)
		{
			$command = "cat $tmpf/for_gms.tbl $tmpf/taxa.set > $tmpf/$ignore_file";
			RunSystem($command); 
		}
		elsif ( $run_gms )
		{
			$ignore_file = "for_gms.tbl";
		}
		else
		{
			$ignore_file = "taxa.set";
		}

		RunMGM( " --seq $sfile --ignore $tmpf/$ignore_file");

		# join predictions from 1-2-3

		if ( $run_hmm )
		{
			foreach my $key (keys %taxa_for_hmm)
			{
				$command = "cat $tmpf/for_hmm/$key.$format >> $out_file";
				RunSystem($command);

				if ($nt_file)
				{
					AppendWithLabelInFastaDefline( "$tmpf/for_hmm/$key.nt", $key, $nt_file );
				}

				if ($aa_file)
				{
					AppendWithLabelInFastaDefline( "$tmpf/for_hmm/$key.aa", $key, $aa_file );
				}
			}
		}

		# add GMS-2 predictions to MGM2

		if ( $run_gms )
		{
			my %ignore_gms;

			if ( $run_hmm )
			{
				%ignore_gms = CreateIgnoreListForGMS( "$tmpf/taxa.set", \%contigs_for_gms );
			}

			foreach my $key (keys %contigs_for_gms)
			{
				next if ( exists $ignore_gms{$key} );

				$command = "cat $tmpf/$key/gms2.$format >> $out_file\n";
				RunSystem($command);

				if ($nt_file)
				{
					$command = "cat $tmpf/$key/nt.fasta >> $nt_file\n";
					RunSystem($command);
				}

				if ($aa_file)
				{
					$command = "cat $tmpf/$key/aa.fasta >> $aa_file\n";
					RunSystem($command);
				}
			}	
		}
	}
}

# return back to original fasta lines
chdir $work_dir;

my %trace;
LoadTrace( $trace_file, \%trace, $defline_status );
GFF_to_original( $out_file, $tmpf, \%trace );
FASTA_to_original( $nt_file, $tmpf, \%trace ) if $nt_file;
FASTA_to_original( $aa_file, $tmpf, \%trace ) if $aa_file;

if ($clean)
{
	 remove_tree($tmpf);
}

system("date") if $verbose;
print "# done\n" if $verbose;
exit 0;

# -----------------------------------------------
sub SelectFastaByTaxa
{
	my $fname = shift;
	my $in_list = shift;

	my %h = LoadList($in_list);

	my $inout_id = '';
	my $OUT;

	open( my $IN, $fname ) or die "error on open file $fname: $!\n";

	while( <$IN> )
	{
		next if /^\s*$/;
		next if /^#/;

		if ( /^>\s*(\S+)/ )
		{
			if ( exists $h{$1} )
			{
				$inout_id = $h{$1};
			}
			else
			{
				$inout_id = '';
			}

			if ( $inout_id )
			{
				open( $OUT, ">>", $inout_id .".seq" ) or die "error on open file $inout_id.seq: $!\n";
				print $OUT $_;
			}
		}
		else
		{
			if ( $inout_id )
			{
				print $OUT $_;
			}
		}
	}
	close $IN;
}
# -----------------------------------------------
sub LoadList
{
	my $fname = shift;

	my %h = ();

	open( my $IN, $fname ) or die "error on open file $fname: $!\n";
	while( <$IN> )
	{
		next if /^\s*$/;
		next if /^#/;

		if ( /^(\S+)\s+(\S+)/ )
		{
			$h{$1} = $2;
		}
	}
	close $IN;

	return %h;
}
# -----------------------------------------------
sub AppendWithLabelInFastaDefline
{
	my $fname = shift;
	my $label = shift;
	my $out_name = shift;

	open( my $IN, $fname ) or die "error on open file $fname: $!\n";
	open( my $OUT, ">>", $out_name ) or die "error on open file $out_name: $!\n";
	while( <$IN> )
	{
		if ( /^>/ )
		{
			chomp $_;
			print $OUT ($_ ." taxa=$label\n");
		}
		else
		{
			print $OUT $_;
		}
	}
	close $IN;
	close $OUT
}
# ------------------------------------------------
sub FastaToInternalFormat
{
	my $fname = shift;
	my $sfile = shift;
	my $trace_file = shift;

	chdir $tmpf;

	my $command = "$probuild --reformat_fasta --in $fname  --out $sfile --trace $trace_file --label _meta";
	RunSystem( $command );

	$sfile = abs_path($sfile);
	$trace_file = abs_path($trace_file);

	chdir $work_dir;

	return ($sfile, $trace_file);
}
# ------------------------------------------------
sub CheckFASTAdefline
{
	my $fname = shift;

	my %h;

	my $is_first_word_good_ID = 1;
	my $sanitized = 0;

	# check if first word in FASTA defline will be good as ID
	
	open( my $IN, $fname ) or die "error on open file $fname $!\n";
	while( my $line = <$IN> )
	{
		if ( $line =~ /^>/ )
		{
			if ( $line =~/^>\s*(\S+)\s*/ )
			{
				my $id = $1;
				my $count_subst = 0;

				$count_subst = $id =~ s/[^\w.\-]/_/g;

				$sanitized = 1 if $count_subst;

				if ( ! exists $h{$id} )
				{
					$h{$id} = 1;
				}
				else
				{
					$is_first_word_good_ID = 0;
					last;
				}
			}
			else { die "error, unexpected format of FASTA defline was found: $line\n"; }
		}
	}
	close $IN;

	if ( $is_first_word_good_ID )
	{
		if ($verbose)
		{
			print "# First word in FASTA definition line is set as sequence identifier\n";;
			print "# Sequence ID/s were modified from original values\n" if $sanitized;
		}

		return ( 1, $sanitized );
	}
	# else ...
	# substitute whitespace symbols in FASTA defline by "_" and check if this creates valid ID

	%h = ();
	$sanitized = 1;

	open( $IN, $fname ) or die "error on open file $fname $!\n";
	while( my $line = <$IN> )
	{
		if ( $line =~ /^>(.+)/ )
		{
			my $id = $1;

			$id =~ s/^\s+//g;
			$id =~ s/\s+$//g;

			$id =~ s/[^\w.\-]/_/g;

			if ( ! exists $h{$id} )
			{
				$h{$id} = 1;
			}
			else { die "error, FASTA sequence ID is not unique: $line\n"; }
		}
	}
	close $IN;

	if ( $verbose)
	{
		print "# Complete FASTA definition line is set as sequence identifier\n";
		print "# Sequence ID/s were modified from original values\n" if $sanitized;
	}

	return ( 2, $sanitized );
}
# ------------------------------------------------
sub LoadTrace
{
	my $fname = shift;
	my $ref = shift;
	my $defline_status = shift;

	open( my $IN, $fname ) or die "error on open file $fname: $!\n";
	while(<$IN>)
	{
		next if /^\s*$/;
		next if /^#/;

		if ( /^>(\S+)\t.+\t>(.+)/ )
		{
			my $key = $1;
			my $value = $2;

			if ( $defline_status == 1 )
			{ 
				if ( $value =~ /^\s*(\S+)\s*/ )
				{
					$ref->{$key} = $1;
					$ref->{$key} =~ s/[^\w.\-]/_/g;
				}
				else { die "error, unexpected FASTA defline format found: $value\n"; }
			}
			else
			{
				if ( $value =~ /^\s*(\S.+\S)\s*/ )
				{
					$ref->{$key} = $1;
					$ref->{$key} =~ s/[^\w.\-]/_/g;
				}
				else { die "error, unexpected FASTA defline format found: $value\n"; }
			}
		}
		else { die "error, unexpected file format found in $fname: $_\n"; }
	}
	close $IN;
}
# ------------------------------------------------
sub FASTA_to_original
{
	my $fname = shift;
	my $scratch = shift;
	my $trace = shift;

	open( my $OUT, ">", "$scratch/tmp_fasta" ) or die "error on open file $scratch/tmp_fasta: $!\n";
	open( my $IN, $fname ) or die "error on open file $fname: $!\n";

	while( my $line = <$IN>)
	{
		if ( $line =~ /^>/ )
		{
			if ( $line =~ /^(>\S+\s+)(\S+)(\s+.+)/ )
			{
				my $pre = $1;
				my $seq_id = $2;
				my $post = $3;

				if ( exists $trace->{$seq_id} )
				{
					$pre =~ s/^>$seq_id/>$trace->{$seq_id}/;

					print $OUT $pre . $trace->{$seq_id} . $post ."\n";
				}
				else { print $OUT $line; }
			}
			else { print $OUT $line; }
		}
		else { print $OUT $line; }
	}

	close $IN;
	close $OUT;

	rename "$scratch/tmp_fasta", $fname;
}
# ------------------------------------------------
sub GFF_to_original
{
	my $fname = shift;
	my $scratch = shift;
	my $trace = shift;

	open( my $OUT, ">", "$scratch/tmp_gff" ) or die "error on open file $scratch/tmp_gff: $!\n";
	open( my $IN, $fname ) or die "error on open file $fname: $!\n";
	while(my $line = <$IN>)
	{
		if ( $line =~ /^##sequence-region\s+(\S+)\s+/ )
		{
			my $seq_id = $1;

			if ( exists $trace->{$seq_id} )
			{
			 	$line =~ s/$seq_id/$trace->{$seq_id}/;
				print $OUT $line;
			}
			else
			{
				print $OUT $line;
			}
		}
		elsif ( $line =~ /#\s+(\d+_meta)\s+/ )
		{
			my $seq_id = $1;

			if ( exists $trace->{$seq_id} )
			{
				$line =~ s/$seq_id/$trace->{$seq_id}/;
				print $OUT $line;
			}
			else
			{
				print $OUT $line;
			}
		}
		elsif( $line =~ /^#/ )
		{
			print $OUT $line;
		}
		elsif ( $line =~ /^(\S+)\t/ )
		{
			my $seq_id = $1;

			if ( exists $trace->{$seq_id} )
			{
				$line =~ s/^$seq_id/$trace->{$seq_id}/;

				if ( $format eq "gtf" )
				{
					$line =~ s/gene_id \"$seq_id/gene_id \"$trace->{$seq_id}/;
					$line =~ s/transcript_id \"$seq_id/transcript_id \"$trace->{$seq_id}/;
				}
				elsif ( $format eq "gff3" )
				{
					$line =~ s/ID=$seq_id/ID=$trace->{$seq_id}/;
					$line =~ s/Parent=$seq_id/Parent=$trace->{$seq_id}/;
				}

				print $OUT $line;
			}
			else
			{
				print $OUT $line;
			}
		}
		else
		{
			print $OUT $line;
		}
	}
	close $IN;
	close $OUT;

	rename "$scratch/tmp_gff", $fname;
}
# ------------------------------------------------
sub CreateIgnoreListForGMS
{
	my $fname = shift;
	my $h_gms_contigs = shift;

	my %ignore;

	open( my $IN, $fname ) or die "error on open file $fname $!\n";
        while( <$IN> )
        {
                next if /^\s*$/;
                next if /^#/;

                if ( /^(\S+)\s+/ )
                {
			my $key = $1;

			if ( exists $h_gms_contigs->{$key} )
			{
				$ignore{$key} += 1;

				print "# ignored GMS for $key\n" if $verbose;
			}
                }
                else { die "error, unexpected line format found: $_"; }
        }
        close $IN;

	return %ignore;
}
# ------------------------------------------------
sub RunHMM
{
	my $seq = shift;
	my $label = shift;
	my $path = shift;

        print "running GeneMark.hmm-2 on $label ...\n" if ($verbose);

        my $command;

	my $g = GetGcodeFromModelFile( "$refmod/$label.mod" );

	if ( $g eq "11" )
	{
		$command = $hmm11;
	}
	elsif ( $g eq "4" )
	{
		$command = $hmm4;
	}
	
	$command .= "  -m $refmod/$label.mod  -o $path/$label.$format";
#al	$command .= " --run_on $path/$label.tbl";
        $command .= "  --seq $seq  --format $format";
        $command .= "  --NT $path/$label.nt" if ( $nt_file );
        $command .= "  --AA $path/$label.aa" if ( $aa_file );

        RunSystem( $command );
}
# ------------------------------------------------
sub GetGcodeFromModelFile
{
	my $fname = shift;

	my $gcode = 0;

	open( my $IN, $fname ) or die "error on open file $fname: $!\n";
	while(<$IN>)
	{
		if ( /GCODE\s+(\d+)\s*/ )
		{
			$gcode = $1;
			last;
		}
	}
	close $IN;

	if ( $gcode != 4 and $gcode != 11 ) { die "error, unsupported genetic code was detected: $gcode\n"; }

	return $gcode;
}
# ------------------------------------------------
sub RunGMS
{
	my $options = shift;

	my $seq = '';
	if ( $options =~ /\s--seq\s+(\S+)/ )
	{
		$seq = $1;
	}

	my $genome_type = "auto";
	if ( $options =~ /\s--genome-type\s+(\S+)/ )
	{
		$genome_type = $1;
	}

	print "running GMS-2 on contig: $seq\n" if $verbose;

	if ( ! -e $seq ) {die"error, file with sequence not found: $seq\n";}

	my $command = "$gms  --seq $seq  --format $format  --output gms2.$format  --genome-type $genome_type";
	$command .= "  --fnn nt.fasta" if ($nt_file);
	$command .= "  --faa aa.fasta" if ($aa_file);

	RunSystem($command);
}
# ------------------------------------------------
sub ParseContigs
{
	my $table = shift;
	my $sequence = shift;

	my $probuild = File::Spec->catfile( $Bin, "probuild" );
	if ( ! -e $probuild ) { die "error, probuild not found\n"; }
	if ( ! -x $probuild ) { die "error, probuild is not executable\n"; }

	my $command = "$probuild  --select_fasta  --seq $sequence  --fasta_list $table  --letters_per_line 60 ";

	RunSystem( $command );
}
# ------------------------------------------------
sub GetContigNames
{
	my $name = shift;

	my %h;

	open( my $IN, $name ) or die "error on open file $name $!\n";
	while( <$IN> )
	{
		next if /^\s*$/;
		next if /^#/;

		if ( /^(\S+)\s+/ )
		{
			if ( exists $h{$1} ) { die "error, duplicated fast ID found in $name\n"; }
			$h{$1} += 1;
		}
		else { die "error, unexpected line format found: $_"; }
	}
	close $IN;

	return %h;
}
# ------------------------------------------------
sub RunMGM_auto
{
	my $option = shift;

	print "running MetaGeneMark-2 auto ...\n" if ($verbose);

	# defoult output names 	
	my $loc_out = $out_file;
	my $loc_format = $format;
	my $loc_nt = $nt_file;
	my $loc_aa = $aa_file;

	# replace defaults, if option is found
	if ( $option =~ /\s--out\s+(\S+)/ )
	{
		$loc_out = $1;
		$option =~ s/\s--out\s+\S+/ /;
	}
	if ( $option =~ /\s--format\s+(\S+)/ )
	{
		$loc_format = $1;
		$option =~ s/\s--format\s+\S+/ /;
	}
	if ( $option =~ /\s--NT\s+(\S+)/ )
	{
		$loc_nt = $1;
		$option =~ s/\s--NT\s+\S+/ /;
	}
	if ( $option =~ /\s--AA\s+(\S+)/ )
	{
		$loc_aa = $1;
		$option =~ s/\s--AA\s+\S+/ /;
	}

	my $command;

	# run with 11
	my $out_11 = $tmpf ."/mgm2_11.gff";
	my $nt_11  = $tmpf ."/mgm2_11.nt";
	my $aa_11  = $tmpf ."/mgm2_11.aa";
	$command = "$Bin/gmhmmp2  -M $Bin/mgm_11.mod  --gid_per_contig";
	$command .= " --out $out_11";
	$command .= " --format $loc_format";
        $command .= " --NT  $nt_11" if $loc_nt;
        $command .= " --AA  $aa_11" if $loc_aa;
        $command .= $option if ( $option );
	RunSystem( $command );

	# run with 4
	my $out_4 = $tmpf ."/mgm2_4.gff";
	my $nt_4  = $tmpf ."/mgm2_4.nt";
	my $aa_4  = $tmpf ."/mgm2_4.aa";
	$command = "$Bin/gmhmmp2  -M $Bin/mgm_4.mod  --gid_per_contig";
        $command .= " --out $out_4";
	$command .= " --format $loc_format";
        $command .= " --NT  $nt_4" if $loc_nt;
        $command .= " --AA  $aa_4" if $loc_aa;
        $command .= $option if ( $option );
	RunSystem( $command );

	# compare 11 vs 4; create combined output
	my %data_11 = ();
	my %data_4 = ();

	LoadMetaData( $out_11, \%data_11 );
	LoadMetaData( $out_4, \%data_4 );
	my %data = ();
	SelectBestGeneticCode( \%data_11, \%data_4, \%data );
	MergeMGM2Predictions( \%data, $out_11, $out_4, $loc_out );
	MergeMGM2_FASTAs( \%data, $nt_11, $nt_4, $loc_nt ) if $loc_nt;
	MergeMGM2_FASTAs( \%data, $aa_11, $aa_4, $loc_aa ) if $loc_aa;
}
# ------------------------------------------------
sub MergeMGM2_FASTAs
{
	my $ref = shift;
        my $fname_11 = shift;
        my $fname_4 = shift;
        my $fname = shift;

	my $found_code = 0;
	my $key = '';

	open( my $OUT, ">", $fname ) or die "error on open file $fname\n";

	open( my $IN_11, $fname_11 ) or die "error on open file $fname_11\n";
	$found_code = 0;
	while(<$IN_11>)
	{
		if ( /^>\s*\S+\s+(\S+)/ )
		{
			$key = $1;
			if ( $ref->{$key} == 11 )
			{
				chomp;
				print $OUT ($_ ." genetic_code=11\n");
				$found_code = 1;
			}
			else
			{
				$found_code = 0;
			}
		}
		elsif ( $found_code )
		{
			print $OUT $_;
		}
	}
	close $IN_11;

	open( my $IN_4, $fname_4 ) or die "error on open file $fname_4\n";
	$found_code = 0;
	while(<$IN_4>)
	{
		if ( /^>\s*\S+\s+(\S+)/ )
		{
			$key = $1;
			if ( $ref->{$key} == 4 )
			{
				chomp;
				print $OUT ($_ ." genetic_code=4\n");
				$found_code = 1;
			}
			else
			{
				$found_code = 0;
			}
		}
		elsif ( $found_code )
		{
			print $OUT $_;
		}
	}
	close $IN_4;

	close $OUT;
}
# ------------------------------------------------
sub MergeMGM2Predictions
{
	my $ref = shift;
	my $fname_11 = shift;
	my $fname_4 = shift;
	my $fname = shift;

	open( my $OUT, ">", $fname ) or die "error on open file $fname\n";

	open( my $IN_11, $fname_11 ) or die "error on open file $fname_11\n";
	while(<$IN_11>)
	{
		if ( /^(\d+_meta)\s+/ )
		{
			print $OUT $_ if ($ref->{$1} == 11);
		}
		elsif ( /#\s+(\d+_meta)\s+/ )
		{
			chomp;
			print $OUT ($_ ."\tgenetic_code 11\n\n") if ($ref->{$1} == 11);
		}
		elsif ( /##sequence-region\s+(\d+_meta)\s+/ )
		{
			if ( defined $ref->{$1} )
			{
				print $OUT ("\n". $_) if ($ref->{$1} == 11);
			}
		}
		elsif ( /\S/ )
		{
			print $OUT $_;
		}
	}
	close $IN_11;

	open( my $IN_4, $fname_4 ) or die "error on open file $fname_4\n";
	while(<$IN_4>)
	{
		if ( /^(\d+_meta)\s+/ )
                {
			print $OUT $_ if ($ref->{$1} == 4);
		}
		elsif ( /#\s+(\d+_meta)\s+/ )
		{
			chomp;
			print $OUT ($_ ."\tgenetic_code 4\n\n") if ($ref->{$1} == 4);
		}
		elsif ( /##sequence-region\s+(\d+_meta)\s+/ )
		{
			if ( defined $ref->{$1} )
			{
				print $OUT ("\n". $_) if ($ref->{$1} == 4);
			}
                }
		elsif( /\S/ )
		{
			print $OUT $_;
		}
	}
	close $IN_4;

	close $OUT;
}
# ------------------------------------------------
sub BestByLogodd
{
	my $score_11 = shift;
	my $score_4 = shift;
	my $contig_length = shift;

	my $best_code = 11;

	if ( $contig_length < $min_gcode_contig )
	{
		return $best_code;
	}

	my $best_score = $score_11;

	if ( $score_4 - $score_11 > $CONST_MIN_DIFF_LOGODD )
	{
		$best_score = $score_4;
		$best_code = 4;
	}

	if( $best_score < $CONST_MIN_LOGODD_FOR_NOT_11 )
	{
		$best_code = 11;
	}

	return $best_code;
}
# ------------------------------------------------
sub BestByAverage
{
	my $score_11 = shift;
	my $score_4 = shift;
	my $contig_length = shift;

	my $best_code = 11;

	if( $contig_length < $min_gcode_contig )
	{
		return $best_code;
	}

	if ( $score_11 > 700 )
	{
		return $best_code;
	}

	if ( $score_4 > $GCODE_CONST_A * $score_11 + $GCODE_CONST_B )
	{
		$best_code = 4;
	}

	return $best_code;
}
# ------------------------------------------------
sub SelectBestGeneticCode
{
	my $ref11 = shift;
	my $ref4 = shift;
	my $ref = shift;

	foreach my $key (keys %{$ref11})
	{
		next if ( ! defined $ref11->{$key}{"logodd"} or  ! defined $ref4->{$key}{"logodd"} ) ;
		
		if ( exists $ref4->{$key} )
		{
#			$ref->{$key} = BestByAverage( $ref11->{$key}{"average"}, $ref4->{$key}{"average"}, $ref11->{$key}{"size"} );
			$ref->{$key} = BestByLogodd( $ref11->{$key}{"logodd"}, $ref4->{$key}{"logodd"}, $ref11->{$key}{"size"} );
		}
		else
			{ die "error, FASTA ID is not found in the second hash: $key\n"; }
	}

	my %h = ();

	$h{'4'} = 0;
	$h{'11'} = 0;

	foreach my $key (keys %{$ref})
	{
		$h{ $ref->{$key} } += 1;
	}

	if ($debug)
	{
		foreach my $key (keys %h)
		{
			print "# gcode assigned $key $h{$key}\n";
		}
	}

	return %h;
}
# ------------------------------------------------
sub LoadMetaData
{
	my $fname = shift;
	my $ref = shift;

	open( my $IN, $fname ) or die "error on open file $fname: $!\n";
	while(<$IN>)
	{
		if ( /^#\s+(\d+_meta)\s+total_logodd\s+(\S+)\s+average_length\s+(\S+)\s+/ )
		{
			my $contig_name = $1;
			my $logodd = $2;
			my $average_length = $3;

			$ref->{$contig_name}{"logodd"} = $logodd;
			$ref->{$contig_name}{"average"} = $average_length;
		}
		elsif ( /##sequence-region\s+(\d+_meta)\s+1\s+(\d+)/ )
		{
			$ref->{$1}{"size"} = $2;
		}
	}
	close $IN;

	if ($debug)
	{
		my $size = scalar (keys %{$ref});
		print "# mgm2 predictions done on file $fname with $size records\n"; 
	}
}
# ------------------------------------------------
sub RunMGM
{
	my $option = shift;

	if ( ! $only_gcode11 )
	{
		RunMGM_auto($option);
		return;
	}

	print "running MetaGeneMark-2 ...\n" if ($verbose);

	my $command; 
	$command  = "$mgm  --out $out_file  --format $format";
	$command .= " --NT $nt_file" if ( $nt_file );
	$command .= " --AA $aa_file" if ( $aa_file );
	$command .= $option if ( $option );

	RunSystem( $command );
}
# ------------------------------------------------
sub GetTaxaNames
{
	my $in_tbl = shift;

	my %h_in;

	open( my $IN, $in_tbl ) or die "error on open file $in_tbl $!\n";
	while( <$IN> )
	{
		next if /^\s*$/;
		next if /^#/;

		if ( /^(\S+)\t(\S+)\s*/ )
		{
			$h_in{$2} += 1
		}
	}
	close $IN;

	return %h_in;
}
# ------------------------------------------------
sub FindCotigsWithHitsToRef
{
	# just check the format 
	
	my $in_tbl = shift;

	my $count = 0;

	open( my $IN, $in_tbl ) or die "error on open file $in_tbl $!\n";
	while( <$IN> )
	{
		next if /^\s*$/;
		next if /^#/;

		if ( /^(\S+)\t(\S+)\s*/ )
		{
			my $sid = $1;
			my $taxa = $2;

			if ( ! -e "$refmod/$taxa.mod" ) { die "error, model file not found: $refmod/$taxa.mod\n"; }

			$count += 1;
		}
		else {die "error, unexpected line format found: $_"; }
	}
	close $IN;

	return $count;
}
# ------------------------------------------------
sub ParseTaxa
{
	my $in_tbl = shift;

	my %h_in;

	open( my $IN, $in_tbl ) or die "error on open file $in_tbl $!\n";
	while( <$IN> )
	{
		next if /^\s*$/;
		next if /^#/;

		if ( /^(\S+)\t(\S+)\s*/ )
		{
			my $sid = $1;
			my $taxa = $2;

			push @{$h_in{$taxa}}, $sid;
		}
		else {die "error, unexpected line format found: $_"; } 
	}
	close $IN;

	my $count = scalar (keys %h_in);

	foreach my $key (keys %h_in)
	{
		my $out_name = $key . ".tbl";

		open( my $OUT, ">", $out_name ) or die "error on open file $out_name $!\n";
		foreach my $id (@{$h_in{$key}})
		{
			print $OUT ($id ."\t". $key ."\n");
		}
		close $OUT
	}

	return $count;
}
# ------------------------------------------------
sub FindCotigsWithMoreThanNgenes
{
	my $name = shift;
	my $out = shift;

	my %h_in;

	open( my $IN, $name ) or die "error on open file $name $!\n";
	while( <$IN> )
	{
		next if /^\s*$/;
		next if /^#/;

		if ( /^(\S+)\t\S+\tCDS\t/ )
		{
			$h_in{$1} += 1;
		}
		else {die "error, unexpected line format found: $_"; }
	}
	close $IN;

	my $count = 0;

	open( my $OUT, ">", $out ) or die "error on open file $out $!\n";
	foreach my $key (keys %h_in)
	{
		if ( $h_in{$key} > $genes_for_gms2 )
		{
			++$count;

			print $OUT ($key ."\t". $h_in{$key} ."\n");
		}
	}
	close $OUT;

	return $count;
}
# ------------------------------------------------
sub CreateTmpWorkspace
{
	my $dir = shift;
	my $label = shift;

	my $tmpdir = tempdir( $label ."_XXXXX", DIR => $dir );
	if ( ! -d  $tmpdir ) { die "error, temporary folder not found: $tmpdir\n"; }
	$tmpdir = ResolvePath($tmpdir);

	return $tmpdir;
}
# ------------------------------------------------
sub RunSystem
{
	my $command = shift;
	print "$command\n" if $debug;
	system( $command ) and die "error on last system call: $command\n";
}
# ------------------------------------------------
sub SetHMM
{
	my $path = shift;
	my $gcode = shift;

	my $gmhmmp2 = File::Spec->catfile( $path, "gmhmmp2" );
	if ( ! -e $gmhmmp2 ) { die "error, gmhmmp2 not found at $path\n"; }
	if ( ! -x $gmhmmp2 ) { die "error, gmhmmp2 is not executable\n"; }

	my $mod;
	if ( $gcode eq "11" )
	{
		$mod = File::Spec->catfile( $path, "mgm_11.mod" );
	}
	elsif ( $gcode eq "4" )
	{
		$mod = File::Spec->catfile( $path, "mgm_4.mod" );
	}
	elsif ( $gcode eq "25" )
	{
		$mod = File::Spec->catfile( $path, "mgm_25.mod" );
	}
	elsif ( $gcode eq "15" )
	{
		$mod = File::Spec->catfile( $path, "mgm_15.mod" );
	}
	else { die "error, unsupported gcode found: $gcode\n"; }

	if ( ! -e $mod ) { die "error, meta parameter file not found: $mod\n"; }

	my $command = "$gmhmmp2 -M $mod --gid_per_contig";

	$cfg{'run_status'}->{"hmm".$gcode} = $command;

	return $command;
}
# ------------------------------------------------
sub SetGMS
{
	my $path = shift;

        my $gms =  File::Spec->catfile( $path, "gms2.pl" );
        if ( ! -e $gms ) { die "error, gms2.pl file is not found\n"; }
        if ( ! -x $gms ) { die "error, gms2.pl file is not executable\n"; }

        my $command = "$gms  --gid ";

	$cfg{'run_status'}->{'gms2'} = $command;

	return $command;
}
# ------------------------------------------------
sub SetMGM
{
	my $path = shift;

	my $gmhmmp2 = File::Spec->catfile( $path, "gmhmmp2" );
	if ( ! -e $gmhmmp2 ) { die "error, gmhmmp2 not found at $path\n"; }
	if ( ! -x $gmhmmp2 ) { die "error, gmhmmp2 is not executable\n"; }

	my $mod = File::Spec->catfile( $path, "mgm_11.mod" );
	if ( ! -e $mod ) { die "error, parameter file mgm_11.mod not found at $path\n"; }

	my $command = "$gmhmmp2 -M $mod --gid_per_contig";

	$cfg{'run_status'}->{'mgm2'} = $command;

	return $command;
}
# ------------------------------------------------
sub CheckBeforeRun
{
	print "CheckBeforeRun\n" if $debug;

	if( !$seq_file ) { die "error, required input file name is missing, check option --seq\n"; }
	if( !$out_file ) { die "error, required output file name is missing, check option --out\n"; }
        
	if (( $format ne "gtf" )and( $format ne "gff3" ))
		{ die "error, unexpected value was specified for option --format: $format\n"; }

	$seq_file = ResolvePath( $seq_file );
	$work_dir = ResolvePath( $work_dir );
	$probuild = ResolvePath( $probuild );
	$parse_blast = ResolvePath( $parse_blast );

	$out_file = File::Spec->rel2abs( $out_file );
	$nt_file  = File::Spec->rel2abs( $nt_file ) if $nt_file;
	$aa_file  = File::Spec->rel2abs( $aa_file ) if $aa_file;

	my %file_names;
	CheckFileNameConflict( $seq_file, \%file_names );
	CheckFileNameConflict( $out_file, \%file_names );
	CheckFileNameConflict( $nt_file,  \%file_names );
	CheckFileNameConflict( $aa_file,  \%file_names );

	if (!$no_taxa)
	{
		if ( ! -e $refmod ) { die "error, path to folder with reference set of models not found: $refmod\n"; }
		if (( ! -e "$refdb.nsq" )and( ! -e "$refdb.00.nsq"  )) { die "error, path to reference sequence database not found: $refdb\n"; }

		$refmod = ResolvePath($refmod);

		my $tmp_name = "$refdb.nsq";
		if ( -e $tmp_name )
		{
			$tmp_name = ResolvePath($tmp_name);
			$tmp_name =~ s/.nsq$//;
		}
		else
		{
			$tmp_name = "$refdb.00.nsq";
			$tmp_name = ResolvePath($tmp_name);
			$tmp_name =~ s/.00.nsq$//;
		}

		$refdb = $tmp_name;
	}

	if ($cores < 1) { die "error, parameter is out of allowed range --cores: $cores\n"; }
	if ($identity < 0 or $identity > 100 ) { die "error, parameter is out of allowed range --identity: $identity\n"; }
	if ($genes_for_gms2 < 30 ) { die "error, parameter is out of allowed range --genes_for_gms2: $genes_for_gms2\n"; }
	if ($evalue < 0) { die "error, value is out of allowed range evalue: $evalue\n"; }
	if ($ani_domain < 0 or $ani_domain > 100 ) { die "error, parameter is out of allowed range --ani_domain: $ani_domain\n"; }

	$cfg{'run_status'}->{'blastn'} = $blastn;
	$cfg{'run_status'}->{'probuild'} = $probuild;
	$cfg{'run_status'}->{'parse_blast'} = $parse_blast;
	$cfg{'run_status'}->{'seq_file'} = $seq_file;
	$cfg{'run_status'}->{'work_dir'} = $work_dir;
	$cfg{'run_status'}->{'bin'} = $Bin;
	$cfg{'run_status'}->{'out_file'} = $out_file;
	$cfg{'run_status'}->{'nt_file'} = $nt_file;
	$cfg{'run_status'}->{'aa_file'} = $aa_file;
	$cfg{'run_status'}->{'refmod'} = $refmod;
	$cfg{'run_status'}->{'refdb'} = $refdb;

	if ($debug)
	{
		print "cfg-run_status\n";
		print Dumper($cfg{'run_status'});
	}
}
# ------------------------------------------------
sub CheckFileNameConflict
{
	my $name = shift;
	my $ref = shift;
	return if (!$name);

	if ( ! exists $ref->{$name} )
	{
		$ref->{$name} = 1;
	}
	else { die "error, file name duplication was detected: $name\n"; }
}
# ------------------------------------------------
sub ResolvePath
{
	my( $name, $path ) = @_;
	return '' if !$name;

	$name = File::Spec->catfile( $path, $name ) if ( defined $path and $path );
	if( ! -e $name ) { die "error, file not found: $name\n"; }
	return abs_path( $name );
}
# ------------------------------------------------
sub SetTmpFolder
{
	my $fname = shift;

	if ( $fname )
        {
                if ( ! -d $fname )
                {
                        mkdir $fname or die "error on create temporary folder $fname\n";
                }
        }
        else
        {
                $fname = CreateTmpWorkspace( $work_dir, "mgm2" );
        }

	$fname = ResolvePath($fname);

	return $fname;
}
# ------------------------------------------------
sub ParseCMD
{
	my $cmd = $0;
	foreach my $str (@ARGV) { $cmd .= ( "\t". $str ); }

	my $opt_results = GetOptions
	(
		'seq=s'     => \$seq_file,
		'out=s'     => \$out_file,
		'nt=s'      => \$nt_file,
		'aa=s'      => \$aa_file,
		'format=s'  => \$format,

		'identity=f'=> \$identity,
		'evalue=f'  => \$evalue,
		'genes_for_gms2=i' => \$genes_for_gms2,
		'ani_domain=f' => \$ani_domain,

		'refdb=s'   => \$refdb,
		'refmod=s'  => \$refmod,
		'tmpf=s'    => \$tmpf,

		'clean'     => \$clean,
		'cores=i'   => \$cores,
		'verbose'   => \$verbose,
		'debug'     => \$debug,

		'mgm2_only' => \$mgm2_only,
		'no_taxa'   => \$no_taxa,
		'no_gms2'   => \$no_gms2,
		'only_gcode11' => \$only_gcode11,
	);

	if( !$opt_results ) { die "error on command line $0\n"; }
	if( @ARGV > 0 ) { die "error, unexpected argument found on command line: @ARGV\n"; }

	$verbose = 1 if $debug;

	print "ParseCMD\n" if $debug;

	# save information for debug
	$cfg{'parameters'}->{'cmd'}      = $cmd;
	$cfg{'parameters'}->{'seq'}      = $seq_file;
	$cfg{'parameters'}->{'out'}      = $out_file;
	$cfg{'parameters'}->{'nt'}       = $nt_file;
	$cfg{'parameters'}->{'aa'}       = $aa_file;
	$cfg{'parameters'}->{'format'}   = $format;
	$cfg{'parameters'}->{'mgm2_only'} = $mgm2_only;
	$cfg{'parameters'}->{'no_taxa'}  = $no_taxa;
	$cfg{'parameters'}->{'no_gms2'}  = $no_gms2;
	$cfg{'parameters'}->{'only_gcode11'} = $only_gcode11,
	$cfg{'parameters'}->{'identity'} = $identity;
	$cfg{'parameters'}->{'evalue'}   = $evalue;
	$cfg{'parameters'}->{'genes_for_gms2'} = $genes_for_gms2;
	$cfg{'parameters'}->{'$ani_domain'} = $ani_domain;
	$cfg{'parameters'}->{'refdb'}    = $refdb;
	$cfg{'parameters'}->{'refmod'}   = $refmod;
	$cfg{'parameters'}->{'tmpf'}     = $tmpf;
	$cfg{'parameters'}->{'clean'}    = $clean;
	$cfg{'parameters'}->{'cores'}    = $cores;
	$cfg{'parameters'}->{'verbose'}  = $verbose;
	$cfg{'parameters'}->{'debug'}    = $debug;
	$cfg{'parameters'}->{'min_gcode_contig'} = $min_gcode_contig;

	if ($debug)
	{
		print "cfg-parameters\n";
		print Dumper($cfg{'parameters'});
	}
}
# ------------------------------------------------
sub Usage
{
        print qq(# -------------------
Usage:  $0  --seq [name]  --out [name]

Required options:

  --seq  [name]
     Name of input file with nucleotide sequence/s in FASTA format.
     First word in FASTA definition line is used as sequence ID.
     Sequence ID must be unique identifier within FASTA file.
  --out  [name]
     Name of output file with coordinates of predicted protein coding genes.

Output options:

  --nt  [name]
     Name of output file with nucleotide sequences of predicted genes
     in FASTA format.
  --aa  [name]
     Name of output file with protein sequences of predicted genes
     in FASTA format.
  --format  [$format]
     Format of output file with gene coordinates.
     GTF and GFF3 are supported: gtf or gff3.

Configuration parameters:

  --identity  [$identity]
     Percent identity threshold in contig to taxa assignment step.
  --evalue  [$evalue]
     E-value for BlastN in contig to taxa assignment step.
  --genes_for_gms2  [$genes_for_gms2]
     Run GeneMarkS-2 on contig if MetaGeneMark-2 
     predicted more than $genes_for_gms2 genes on that contig.
  --refdb  [$refdb]
     Path to database with reference gene sequences.
     Database must be in NCBI blast format.
  --refmod  [$refmod]
     Path to folder with taxa specific gene prediction parameter files
     for GeneMark.hmm-2
  --tmpf  [name]
     Run calculations in specified temporary folder.
     By default temporary folder is created with unique random name
     for each pipeline run.
  --clean
     Remove temporary files and folders.

Algorithm options:

  --mgm2_only
     Predict genes only by MetaGeneMark-2.
  --no_taxa
     Skip GeneMark.hmm-2 predictions with taxa specific
     gene finding parameters.
  --no_gms2
     Skip GeneMarkS-2 predictions.
  --only_gcode11
     Use only genetic code 11

Other parameters:

  --cores  [$cores]
     Number of cores/threads to use by BlastN.
  --verbose
     Output on STDOUT pipeline runtime information.

Developer options:
  --debug
Version  $VERSION
# -------------------
);
	exit 1;
}
# ================== END sub =====================

