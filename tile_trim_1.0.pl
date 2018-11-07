#!/public/agis/znjy_group/xianglilingzi/software/miniconda/bin/perl

=head1 NAME

tile_trim.pl - remove low quality tile sequences from a FASTQ file

=head1 SYNOPSIS

USAGE: 

  tile_trim.pl <input_file> <tile_number> [ <qc_file> ]

Optional arguments:

    --input_file,-i                     Path to a low quality FASTQ file
    --tile_number,-n                    A list of low quality tile sequences number, comma-separated values (int[,int])
    --qc_file,-q                        Path to a FASTQC quality control result file (fastqc_data.txt)
    --tile_quality_threshold,-c	        Tile quality threshold for saving sequences (default=-7)
    --output_file,-o                    Path to a remaining file
    --discard_file,-d                   Path to a discarded file
    --help,-h                           Help message

=head1 OPTIONS

B<--input_file,-i>

    A low quality FASTQ file.

B<--tile_number,-n>

    A list of low quality tile sequences number, comma-separated values.

B<--qc_file,-q>

    FASTQC quality control result file (fastqc_data.txt) of the low quality FASTQ file to be handled.

B<--tile_quality_threshold,-c>

    Tile quality threshold for saving sequences.

B<--output_file,-o>

    The remaining sequences exclude the low quality tile sequences.

B<--discard_file,-d>

    The discarded sequences (the low quality tile sequences).

B<--help,-h>

    This help message

=head1  CONTACT

    Lingzi Xiangli
    xianglilingzi@icloud.com

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use PerlIO::gzip;

my (%options, @index, %hash);
my ($remain_base, $remain_read, $discard_base, $discard_read)=(0,0,0,0);
GetOptions(\%options,
	'input_file|i=s',
	'tile_number|n=s',	
	'qc_file|q=s',
	'tile_quality_threshold|c=i',	
	'output_file|o=s',
	'discard_file|d=s',
	'help|h'
) || pod2usage();

## display documentation
if ($options{'help'}){
	pod2usage({-exitval => 0, -verbose => 2, -output => \*STDERR});
}

## make sure everything passed was peachy
&check_parameters();

if (defined $options{qc_file}){
	get_index();
}
if (defined $options{tile_number}){
	push @index, split(/,/, $options{tile_number});
}

my $first_line = <IN>;
my $delimiter = (split(/:/, $first_line, 2))[0];
$/="$delimiter";
(my $first_seq = <IN>) =~ s/$delimiter//g;
my $true_seq = (split(/\n/, $first_seq,2))[0];
my $j=0;
foreach my $index (@index){
	if ($first_line =~ /:$index:\d+:\d+\s+/){
		$j=1;
		last;
	}
}
if ($j == 0){
	if (defined $options{output_file}){
		print OUT $first_line.$first_seq;
	}else{
		print STDOUT $first_line.$first_seq;
	}
	$remain_base += length($true_seq);
	$remain_read++;
}else{
	print DISCARD $first_line.$first_seq if (defined $options{discard_file});
	$discard_base += length($true_seq);
	$discard_read++;
}
while (<IN>){
	$_ =~ s/$delimiter//g;
	my $true_seq = (split(/\n/, $_, 3))[1];
	$j=0;
	foreach my $index (@index){
		if (/:$index:\d+:\d+\s+/){
			$j=1;
			last;
		}
	}
	if ($j == 0){
		if (defined $options{output_file}){
			print OUT $delimiter.$_;
		}else{
			print STDOUT $delimiter.$_;
		}
		$remain_base += length($true_seq);
		$remain_read++;
	}else{
		print DISCARD $first_line.$first_seq if (defined $options{discard_file});
		$discard_base += length($true_seq);
		$discard_read++;
	}
}
$/="\n";
close IN;
close OUT;

my $sum_read = $remain_read + $discard_read;
my $sum_base = $remain_base + $discard_base;
my $discard_read_percent = $discard_read / $sum_read *100;
my $discard_base_percent = $discard_base / $sum_base *100;
printf STDERR "%-12s%12d%s%21d%s\n", 'Inputs:', $sum_read, ' reads', $sum_base, ' bases';
printf STDERR "%-12s%12d%s%.3f%s%12d%s%.3f%s\n", 'Removed:', $discard_read, ' reads (', $discard_read_percent, ')', $discard_base, ' bases (', $discard_base_percent, ')';
printf STDERR "%-12s%12d%s%.3f%s%12d%s%.3f%s\n", 'Result:', $remain_read, ' reads (', 100-$discard_read_percent, ')', $remain_base, ' bases (', 100-$discard_base_percent, ')';

sub check_parameters {
	die "you must input a --input_file" unless (defined $options{input_file});
	unless (defined $options{tile_number} || defined $options{qc_file}){
		die "you must specify either --tile_number or --qc_file";
	}
	if ($options{input_file} =~ /\.gz$/){
		open IN, "<:gzip", "$options{input_file}" or die $!;
	}else{
		open IN, "$options{input_file}" or die $!;
	}
	if (defined $options{output_file}){
		if ($options{output_file} =~ /\.gz$/){
			open OUT, ">:gzip", "$options{output_file}";
		}else{
			open OUT, ">$options{output_file}";
		}
	}
	$options{tile_quality_threshold} = -5 unless (defined $options{tile_quality_threshold});
}

sub get_index{
	open QC, "$options{qc_file}" or die "Can't open file $options{qc_file}.\n";
	$/=">>";<QC>;
	while (<QC>){
		if (/^Per tile sequence quality/){
			my @information = split(/\n/, $_);
			foreach (@information){
				if (/^\d/){
					my @array = split(/\s+/, $_);
					$hash{$array[0]}{$array[1]}=$array[2];
				}
			}
			last;
		}
	}
	$/="\n";
	close QC;
	my @tmpindex;
	foreach my $k (keys %hash){
		foreach my $subk (keys %{$hash{$k}}){
			if ($hash{$k}{$subk}<$options{tile_quality_threshold}){
				push @tmpindex, $k;
			}
		}
	}
	my %h;
	@h{@tmpindex}=();
	@index = keys %h;
}
