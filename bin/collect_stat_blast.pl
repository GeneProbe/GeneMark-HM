#!/usr/bin/perl
# ---------------------
# Alex Lomsadze
#
# collect stat from verbose output of mgm2.pl run
# ---------------------

use strict;
use warnings;
use Data::Dumper;

if ( $#ARGV < 0 )
{
	print "Usage: $0 logfile/s\n";
	exit 1;
}

my %d;

foreach my $name (@ARGV)
{
	ParseFile($name, \%d);
}

print "\# 1 genes_predicted_by_mgm2\n";
print "\# 2 contigs_with_predicted_genes\n";
print "\# 3 genes_with_blast_hits\n";
print "\# 4 contigs_with_blast_hits\n";
print "\# 5 genes_with_blast_hits_above_threshold\n";
print "\# 6 contigs_with_blast_hits_above_thresholds\n";
print "\# 7 genes_in_contigs_passed_threshold:\n";
print "\# 8 chimera_contigs_with_blast_hits_above_thresholds\n";
print "\# 9 total_genes_in_chimera_contigs\n";
print "\#10 genes_in_chimera_contigs_above_threshold\n";

print "# 1\t2\t3\t4\t5\t6\t7\t8\t9\t10\n";

foreach my $name (@ARGV)
{
	print $d{$name}{"genes_in_blast_input"} ."\t"; 
	print $d{$name}{"contigs_with_genes_total"} ."\t";
	print $d{$name}{"genes_in_with_hits"} ."\t";
	print $d{$name}{"contigs_with_hits"} ."\t";
	print $d{$name}{"genes_with_hits_passed_thresholds"} ."\t";
	print $d{$name}{"contigs_passed_threshold"} ."\t";
	print $d{$name}{"genes_in_contigs_passed_threshold"} ."\t";
	print $d{$name}{"chimera_contigs_passed_threshold"} ."\t";
	print $d{$name}{"total_genes_in_chimera_contigs"} ."\t";
	print $d{$name}{"genes_in_chimera_contigs_above_threshold"} ."\n";
}

# ---------------------
sub ParseFile
{
	my $fname = shift;
	my $h = shift;

	open( my $IN, $fname) or die "error on open file $fname: $!\n";
	while(my $line = <$IN>)
	{
		if ($line =~ / genes in blast input: (\d+)/)
		{
			$h->{$fname}{"genes_in_blast_input"} = $1;

			$line = <$IN>;
			if ($line =~ / contigs with genes total: (\d+)/)
			{
				$h->{$fname}{"contigs_with_genes_total"} = $1;
			}
			else {die "error in file format: $line\n";}

			$line = <$IN>;

			$line = <$IN>;
			if ($line =~ / genes in with hits: (\d+)/)
			{
				$h->{$fname}{"genes_in_with_hits"} = $1;
			}
			else {die "error in file format: $line\n";}

			$line = <$IN>;
			if ($line =~ / contigs with hits: (\d+)/)
			{
				$h->{$fname}{"contigs_with_hits"} = $1;
			}
			else {die "error in file format: $line\n";}

			$line = <$IN>;
			$line = <$IN>;

			$line = <$IN>;
			if ($line =~ / genes with hits passed thresholds: (\d+)/)
			{
				$h->{$fname}{"genes_with_hits_passed_thresholds"} = $1;
			}
			else {die "error in file format: $line\n";}

			$line = <$IN>;
			if ($line =~ / contigs passed threshold: (\d+)/)
			{
				$h->{$fname}{"contigs_passed_threshold"} = $1;
			}
			else {die "error in file format: $line\n";}

			$line = <$IN>;
			if ($line =~ / genes in contigs passed threshold: (\d+)/)
			{
				$h->{$fname}{"genes_in_contigs_passed_threshold"} = $1;
			}
			else {die "error in file format: $line\n";}

			$line = <$IN>;

			$line = <$IN>;
			if ($line =~ / chimera contigs passed threshold: (\d+)/)
			{
				$h->{$fname}{"chimera_contigs_passed_threshold"} = $1;
			}
			else {die "error in file format: $line\n";}


			$line = <$IN>;
			if ($line =~ / total genes in chimera contigs: (\d+)/)
			{
				$h->{$fname}{"total_genes_in_chimera_contigs"} = $1;
			}
			else {die "error in file format: $line\n";}

			$line = <$IN>;
			if ($line =~ / genes in chimera contigs above threshold: (\d+)/)
			{
				$h->{$fname}{"genes_in_chimera_contigs_above_threshold"} = $1;
			}
			else {die "error in file format: $line\n";}

			last;
		}
	}
	close $IN;

#	print Dumper($h);
}
# ---------------------

