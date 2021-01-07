#!/usr/bin/perl
# ---------------------------
# Alex Lomsadze
#
# parse output from BlastN to get the taxa id
# ---------------------------

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

# ------------------------
my $in_blast = '';
my $out_tbl = '';
my $identity_th = 95;
my $length_th = 90;
my $verbose = '';
my $exclude = '';
my $delta_ani = 2;
my $debug = 0;

ParseCMD();
# ------------------------
my %exclude_this;
if ($exclude)
{
	%exclude_this = LoadList($exclude); 
}

my %seq_to_taxa;
ParseBlast( $in_blast, \%seq_to_taxa );

PrintOut( $out_tbl, \%seq_to_taxa );

exit 0;

# ---------------------------
sub PrintOut
{
	my $fname = shift;
	my $h_seq_to_taxa = shift;

	open( my $OUT, ">", $fname ) or die "error on open file $fname: $!\n";
	foreach my $key (keys %{$h_seq_to_taxa})
	{
		print $OUT ( $key ."\t". $h_seq_to_taxa->{$key} ."\n" );
	}
	close $OUT;
}
# ---------------------------
sub ParseBlast
{
	my $fname = shift;
	my $ref = shift;

	my $hit_order;	

	my $count_genes = 0;
	my $count_genes_with_hits = 0;
	my $count_hits_passed = 0;

	# contig id  -> number of genes in contig
	my %contigs;
	# contig id  -> number of genes with hits in contig
	my %contigs_with_hits;
	# contig id  -> number of genes with hits above threshold in contig
	my %contigs_with_above_hits;

	my %chimera;
	my %contigs_genes;

	my %h_contig_gene_hit;

	my $min_identity = 100;
	my $min_length = 10000000000;

	open( my $IN, $fname ) or die "error on open file $fname: $!\n";

	my $id = '';

	# works with %exclude_this 
	# set it tp "0" for exclude block
	my $keep = 1;

	while( <$IN> )
	{
		next if (/^\s*$/);

		if ( /# Query:\s+(\S+)\s*/ )
		{
			$id = $1;
			if ( exists $exclude_this{$id} )
			{
				$keep = 0;
				next;
			}
			else
			{
				$keep = 1;
			}

			$count_genes += 1;
			$hit_order = 0;

			if ( $id =~ /(\S+)_\d+$/ )
			{
				$id = $1;
			}
			else { die "error, unexpected format for gene ID in Query was detected: $id\n"; }

			# count number of genes in each contig

			$contigs{$id} += 1;
		}

		next if ( !$keep );

		next if /^#/;

		my @arr = split;

		if ( scalar @arr != 12 ) { die "error, unxpected line format found: $_\n"; }

                # 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send   evalue bitscore',
                #  0      1      2      3      4        5       6      7    8      9      10     11

		# Somehow BLAST is changing the target ID - uppercases it in some cases
		if ( $arr[1] =~ /^B/ )
		{
			$arr[1] =~ s/^B/b/;
		}
		if ( $arr[1] =~ /^T/ )
		{
			$arr[1] =~ s/^T/t/;
		}

		if ($verbose)
		{
			$min_identity = $arr[2] if ( $min_identity > $arr[2] );
			$min_length = $arr[3] if ( $min_length > $arr[3] );
		}

		my $sid = $arr[0];

		if ( $sid =~ /(\S+)_\d+$/ )
		{
			$sid = $1;
		}
		else { die "error, unexpected format for gene ID was detected: $sid\n"; }

		if ( $id ne $sid ) { die "error, mismatch in id and sid: $id $sid\n"; }
		
		$hit_order += 1;

		next if ( $hit_order > 1 );

		# count number of genes with hits per each contig with at least one hit

		$contigs_with_hits{$sid} += 1;
		$count_genes_with_hits += 1;

		if ( $arr[2] > $identity_th and $arr[3] > $length_th )
		{
			$contigs_with_above_hits{$id} += 1;
			$count_hits_passed += 1;

			if ( $arr[0] =~ /^(\d+_meta)_(\d+)/ )
			{
				$h_contig_gene_hit{$1}{$2} = $arr[1];
			}

			if ( ! exists $ref->{$sid} )
			{
				$ref->{$sid} = $arr[1];
				$contigs_genes{$sid} = $arr[1];
			}
			else
			{
				$contigs_genes{$sid} .= ( "\t". $arr[1] );

				if (( exists $chimera{$sid} )or( $ref->{$sid} ne $arr[1] ))
				{
					$chimera{$sid} += 1;
				}
			}
		}

#		$ref->{$count}{'ref'}   = $arr[1];
#		$ref->{$count}{'mis'}   = $arr[4];
#		$ref->{$count}{'L_q'}   = $arr[6];
#		$ref->{$count}{'R_q'}   = $arr[7];
#		$ref->{$count}{'L_ref'} = $arr[8];
#		$ref->{$count}{'R_ref'} = $arr[9];
#		$ref->{$count}{'blast'} = $_;

#               print Dumper( $ref->{$count} );

	}
	close $IN;


	foreach my $key (keys %chimera)
	{
		my %count_sp_in_contig = ();
		my @sp = split( '\t',  $contigs_genes{$key} );

		foreach my $label (@sp)
		{
			$count_sp_in_contig{$label} += 1;
		}

		my @sp_sorted = sort{ $count_sp_in_contig{$b}<=> $count_sp_in_contig{$a} } keys %count_sp_in_contig;

		$ref->{$key} = $sp_sorted[0];
	}

	# reread the file and look for fix of chimera
	open( $IN, $fname ) or die "error on open file $fname: $!\n";
	my $rescan_this = 0;
	my %score = ();
	while( <$IN> )
	{
		if ( /# Query:\s+(\d+_meta)_/ )
		{
			$id = $1;
			
			if ( exists $chimera{$id} )
			{
				$rescan_this = 1;
			}
			else
			{
				$rescan_this = 0;
			}

			%score = ();
		}

		next if (!$rescan_this);

		next if (/^\s*$/);
		next if /^#/;

		my @arr = split;

		if ( scalar @arr != 12 ) { die "error, unxpected line format found: $_\n"; }

                # 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send   evalue bitscore',
                #  0      1      2      3      4        5       6      7    8      9      10     11

		my $re_seq_id = '';
		my $re_gene_id = 0;

                if ( $arr[0] =~ /^(\S+)_(\d+)$/ )
                {
                        $re_seq_id = $1;
			$re_gene_id = $2;
                }

		if ( exists $h_contig_gene_hit{$re_seq_id}{$re_gene_id} )
		{
			next if ( $h_contig_gene_hit{$re_seq_id}{$re_gene_id} eq $ref->{$re_seq_id} );

			if ( ! exists $score{$re_gene_id} )
			{
				$score{$re_gene_id} = $arr[2];
				print "Working on $re_seq_id  $re_gene_id  $arr[1] $arr[2] $arr[3]\n" if $debug;
			}

			

			if ( ($arr[1] eq $ref->{$re_seq_id}) and $arr[2] > $identity_th and $arr[3] > $length_th )
			{
				if ( $score{$re_gene_id} - $arr[2] < $delta_ani )
				{
					print "# changing from: $h_contig_gene_hit{$re_seq_id}{$re_gene_id}  $score{$re_gene_id}\n" if $debug;
					print "# changing to: $re_seq_id $re_gene_id  $arr[1] $arr[2] $arr[3]\n" if $debug;

					$h_contig_gene_hit{$re_seq_id}{$re_gene_id} = $arr[1];
				}
			}
		}
	}


	# update chimera counts
	my %updated_chimera = ();

	foreach my $skey (keys %chimera)
	{
		foreach my $gkey (keys %{$h_contig_gene_hit{$skey}})
		{
			if ( $seq_to_taxa{$skey} ne $h_contig_gene_hit{$skey}{$gkey} )
			{
				$updated_chimera{$skey} += 1;
			}
		}
	}

	%chimera = %updated_chimera;


	print "\n" if $verbose;

	if ($verbose)
	{
		print "#\n";
		print "# threshold identity: $identity_th\n";
		print "# threshold length: $length_th\n";
		print "# minimum identity detected: $min_identity\n";
		print "# minimum length detected: $min_length\n";
		print "#\n";

		if ( SumHash(\%contigs) !=  $count_genes ) {die "error 1\n";}

		print "# genes in blast input: $count_genes\n";
		print "# contigs with genes total: ". (scalar (keys %contigs)) ."\n";
		print "#\n";

		if ( SumHash(\%contigs_with_hits) != $count_genes_with_hits  ) {die "error 2\n";}

		my $number_of_genes_in_contigs_with_hits = 0;
		foreach my $key (keys %contigs_with_hits)
		{
			$number_of_genes_in_contigs_with_hits += $contigs{$key}
		}

		print "# genes in with hits: $count_genes_with_hits\n";
		print "# contigs with hits: ". (scalar (keys %contigs_with_hits)) ."\n";
		print "# total number of genes in these contigs: $number_of_genes_in_contigs_with_hits\n";
		print "#\n";

		if ( SumHash(\%contigs_with_above_hits) != $count_hits_passed  ) {die "error 3\n";}

		print "# genes with hits passed thresholds: $count_hits_passed\n";
		print "# contigs passed threshold: ". (scalar (keys %{$ref})) ."\n";

		my $genes_in_passed = 0;
		foreach my $key (keys %{$ref})
		{
			$genes_in_passed += $contigs{$key};
		}
		print "# genes in contigs passed threshold: $genes_in_passed\n";

		print "#\n";
		print "# chimera contigs passed threshold: ". (scalar (keys %chimera)) ."\n";

		my $genes_in_chimera = 0;
		my $genes_in_chimera_above_th = 0;
		foreach my $key (keys %chimera)
		{
			$genes_in_chimera += $contigs{$key};
			$genes_in_chimera_above_th += $contigs_with_above_hits{$key};
		}
		print "# total genes in chimera contigs: $genes_in_chimera\n";
		print "# genes in chimera contigs above threshold: $genes_in_chimera_above_th\n";

		foreach my $key (keys %chimera)
		{
			print "CH\t$key\t". $contigs_genes{$key} ."\n";
			print "BCH\t$ref->{$key}\n";
		}
	}
}
# ---------------------------
sub SumHash
{
	my $h = shift;
	my $sum = 0;
	foreach my $key (keys %{$h})
	{
		$sum += $h->{$key};
	}
	return $sum;
}
# ---------------------------
sub LoadList
{
	my $fname = shift;

	my %h;

	open( my $IN, $fname) or die "error on open file $fname: $!\n";
	while(<$IN>)
	{
		next if /^#/;
		next if /^\s*$/;

		if ( /^(\S+)\s*/ )
		{
			$h{$1} += 1;
		}
	}
	close $IN;

	return %h;
}
# ---------------------------
sub ParseCMD
{
	if( @ARGV == 0 ) { Usage(); exit 1; }

	my $opt_results = GetOptions
	(
		'in=s'         => \$in_blast,
		'out=s'        => \$out_tbl,
		'identity=f'   => \$identity_th,
		'delta_ani=f'  => \$delta_ani,
		'length=i'     => \$length_th,
		'verbose'      => \$verbose,
		'exclude=s'    => \$exclude,
	);

	if( !$opt_results ) { die "error on command line\n"; }
	if( @ARGV > 0 ) { die "error, unexpected argument found on command line: @ARGV\n"; }

	if ( !$in_blast ) { die "error, required option not found: --in\n"; }
	if ( !$out_tbl ) { die "error, required option not found: --out\n"; }

	if ( $in_blast eq $out_tbl ) { die "error, ouput and input names match\n"; }
	if ( $identity_th < 0 or $identity_th > 100 ) { die "error, value of option --identity is out out alllowed range 0..100:  $identity_th\n"; }
};
# ------------------------
sub Usage
{
	my $txt =
"# -----------
Usage: $0  --in [file_name]  --out [file_name]

Parse from BlastN output information about taxa id

  --identity  [$identity_th] threshold of identiry levele for taxa call
  --length    [$length_th] threshold of minmum hit length
  --delta_ani [$delta_ani] threshold for second path
  --exclude   [file_name] exclude parsing of genes from this list 

  --verbose
# -----------
";
	print $txt;
	exit 1;
};
# ------------------------

