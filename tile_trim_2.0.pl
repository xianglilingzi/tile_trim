#!/public/agis/znjy_group/xianglilingzi/software/miniconda/bin/perl

=head1 NAME

tile_trim.pl - remove low quality tile sequences from a FASTQ file

=head1 SYNOPSIS

USAGE: 

  tile_trim.pl <input_file> <tile_number> [ <qc_file> ]

Optional arguments:

    --input_file,-i                     Path to a low quality FASTQ file
    --tile_coordinate,-c                A list of low quality tile sequences coordinate, comma-separated values (int:int-int,)
    --qc_file,-q                        Path to a FASTQC quality control result file (fastqc_data.txt)
    --tile_quality_threshold,-t         Tile quality threshold for saving tile sequences (default=-5)
    --whole,-w                          Removed all of the reads of the low quality tile
    --left_position,-pl                 Trim the left bases of this position (1-based)
    --right_position,-pr                Trim the right bases of this position (1-based)
    --min_length,-l                     Length threshold for filtering reads (default=0)
    --output_file,-o                    Path to a remaining file
    --discard_file,-d                   Path to a discarded file
    --log_file,-log                     Path to a log file
    --help,-h                           Help message

=head1 OPTIONS

B<--input_file,-i>

    A low quality FASTQ file.

B<--tile_coordinate,-c>

    A list of low quality tile sequence coordinates. You can only specify tile sequence number. You can also specify bases number,1-based on-base coordinates. The format is just like "1101:10-14,1102:1", "1103:45-49", "1101,1102,1103" or "1104". When you only specify tile sequence number, all of the reads of the specify tile will be removed.

B<--qc_file,-q>

    FASTQC quality control result file (fastqc_data.txt) of the low quality FASTQ file to be handled.

B<--tile_quality_threshold,-t>

    Tile quality threshold for saving sequences. This only works for the FASTQC quality control result file. The default value is -5.

B<--whole,-w>

    When you only specify the "whole" option, all of the reads of the low quality tile will be removed.

B<--left_position,-pl>

    Trim the left bases of this position (saving this position, 1-based).

B<--right_position,-pr>

    Trim the right bases of this position (saving this position, 1-based).

B<--min_length,-l>

    Read shorter than this value will be removed. The default value is 0.

B<--output_file,-o>

    The remaining sequences exclude the low quality tile sequences.

B<--discard_file,-d>

    The discarded sequences (the low quality tile sequences).

B<--log_file,-log>

    The Log file.

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

my (%options, @index, $axis, %hash, %sequence_length, %cut, %leftcut, %rightcut, %leftcutmax, %rightcutmin);
my ($remain_base, $remain_read, $discard_base, $discard_read, $filter_base, $filter_read, $trim_base, $trim_read)=(0,0,0,0,0,0,0,0);
GetOptions(\%options,
	'input_file|i=s',
	'tile_coordinate|c=s',
	'qc_file|q=s',
	'tile_quality_threshold|t=i',
	'whole|w',
	'left_position|pl=i',
	'right_position|pr=i',
	'min_length|l=i',
	'output_file|o=s',
	'discard_file|d=s',
	'log_file|log=s',
	'help|h'
) || pod2usage();

## display documentation;
if ($options{'help'}){
	pod2usage({-exitval => 0, -verbose => 2, -output => \*STDERR});
}

## make sure everything passed was peachy;
&check_parameters();

## get low quality tile sequences coordinate;
if (defined $options{qc_file}){
	&get_coordinate();
}
if (defined $options{tile_coordinate}){
	push my @tmp, split(/,/, $options{tile_coordinate});
	foreach my $t (@tmp){
		push @index, (split(/:/, $t))[0];
		$cut{(split(/:/, $t))[0]} .= (split(/:/, $t))[1]."," if ((split(/:/, $t))[1]);
	}
}

## main function; filter low quality tile sequences;
(my $first_line = <IN>) =~ s/\n//g;
my $delimiter = (split(/:/, $first_line, 2))[0];
$/="$delimiter";
(my $first_seq = <IN>) =~ s/$delimiter//g;
my @true_seq = split(/\n/, $first_seq);
&filter_function($first_line, \@true_seq);
while (<IN>){
	$_ =~ s/$delimiter//g; 
	my @fastq = split /\n/, $_;
	my $first = $delimiter.shift(@fastq);
	&filter_function($first, \@fastq);
}
$/="\n";

## Statistics 
my $sum_read = $remain_read + $discard_read;
my $sum_base = $remain_base + $discard_base;
my $discard_read_percent = $discard_read / $sum_read *100;
my $discard_base_percent = $discard_base / $sum_base *100;
my $filter_read_percent = $filter_read / $sum_read *100;
my $filter_base_percent = $filter_base / $sum_base *100;
my $trim_read_percent = $trim_read / $sum_read *100;
my $trim_base_percent = $trim_base / $sum_base *100;

if (defined $options{log_file}){	
	printf LOG "%-12s%12d%s%21d%s\n", 'Inputs:', $sum_read, ' reads', $sum_base, ' bases';
	printf LOG "%-12s%12d%s%.3f%s%12d%s%.3f%s\n", 'Trimed:', $trim_read, ' reads (', $trim_read_percent, ')', $trim_base, ' bases (', $trim_base_percent, ')';
	printf LOG "%-12s%12d%s%.3f%s%12d%s%.3f%s\n", 'Filtered:', $filter_read, ' reads (', $filter_read_percent, ')', $filter_base, ' bases (', $filter_base_percent, ')';
	printf LOG "%-12s%12d%s%.3f%s%12d%s%.3f%s\n", 'Removed:', $discard_read, ' reads (', $discard_read_percent, ')', $discard_base, ' bases (', $discard_base_percent, ')';
	printf LOG "%-12s%12d%s%.3f%s%12d%s%.3f%s\n", 'Result:', $remain_read, ' reads (', 100-$discard_read_percent, ')', $remain_base, ' bases (', 100-$discard_base_percent, ')';
}
printf STDERR "%-12s%12d%s%21d%s\n", 'Inputs:', $sum_read, ' reads', $sum_base, ' bases';
printf STDERR "%-12s%12d%s%.3f%s%12d%s%.3f%s\n", 'Trimed:', $trim_read, ' reads (', $trim_read_percent, ')', $trim_base, ' bases (', $trim_base_percent, ')';
printf STDERR "%-12s%12d%s%.3f%s%12d%s%.3f%s\n", 'Filtered:', $filter_read, ' reads (', $filter_read_percent, ')', $filter_base, ' bases (', $filter_base_percent, ')';
printf STDERR "%-12s%12d%s%.3f%s%12d%s%.3f%s\n", 'Removed:', $discard_read, ' reads (', $discard_read_percent, ')', $discard_base, ' bases (', $discard_base_percent, ')';
printf STDERR "%-12s%12d%s%.3f%s%12d%s%.3f%s\n", 'Result:', $remain_read, ' reads (', 100-$discard_read_percent, ')', $remain_base, ' bases (', 100-$discard_base_percent, ')';

close IN;
close OUT if (defined $options{output_file});
close DISCARD if (defined $options{discard_file});
close LOG if (defined $options{log_file});


#########################
#######Sub Function######
#########################

sub check_parameters {
	die "you must input a --input_file" unless (defined $options{input_file});
	unless (defined $options{tile_coordinate} || defined $options{qc_file}){
		die "you must specify either --tile_coordinate or --qc_file";
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
	if (defined $options{discard_file}){
		if ($options{discard_file} =~ /\.gz$/){
			open DISCARD, ">:gzip", "$options{discard_file}";
		}else{
			open DISCARD, ">$options{discard_file}";
		}
	}
	if (defined $options{log_file}){
		if ($options{log_file} =~ /\.gz$/){
			open LOG, ">:gzip", "$options{log_file}";
		}else{
			open LOG, ">$options{log_file}";
		}
	}
	$options{tile_quality_threshold} = -5 unless (defined $options{tile_quality_threshold});
	$options{min_length} = 0 unless (defined $options{min_length});
}

sub get_coordinate {
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
	foreach my $k (keys %hash){
		my $i=0;
		my @basecoordinate;
		foreach my $subk (keys %{$hash{$k}}){
			push @basecoordinate, split(/-/, $subk);
			if ($hash{$k}{$subk}<$options{tile_quality_threshold}){
				$cut{$k} .= $subk.",";
				$i=1;
			}
		}
		if ($i == 1){
			$sequence_length{$k} = (sort {$a<=>$b} @basecoordinate)[-1];
			push @index, $k;
		}
	}
}

sub get_base_coordinate {
	my $x = shift;
	my $base = shift;
	my @basearray = split /,/, $base;
	foreach (@basearray){
		my @sortposition = sort {$a<=>$b} split(/-/, $_);
		if ($sortposition[-1] <= $sequence_length{$x}/2){
			$leftcut{$x}.= $sortposition[-1]."+";
		}else{
			$rightcut{$x}.= $sortposition[0]."+";
		}
	}
	$leftcutmax{$x} = (sort {$a<=>$b} (split /\+/, $leftcut{$x}))[-1] if ($leftcut{$x});
	$rightcutmin{$x} = (sort {$a<=>$b} (split /\+/, $rightcut{$x}))[0] if ($rightcut{$x});
}

sub filter_function {
	my ($title, $seq) = @_;
	my $j=0;
	foreach (@index){
		if ($title =~ /:$_:\d+:\d+\s+/){
			$axis = $_;
			$j=1;
			last;
		}
	}
	if ($j == 0){
		if (length(@$seq[0])>=$options{min_length}){
			if (defined $options{output_file}){
				print OUT $title."\n".@$seq[0]."\n".@$seq[1]."\n".@$seq[2]."\n";
			}else{
				print STDOUT $title."\n".@$seq[0]."\n".@$seq[1]."\n".@$seq[2]."\n";
			}
			$remain_base += length(@$seq[0]);
			$remain_read++;
		}else{
			print DISCARD $title."\n".@$seq[0]."\n".@$seq[1]."\n".@$seq[2]."\n" if (defined $options{discard_file});
			$filter_base += length(@$seq[0]);
			$filter_read++;
			$discard_base += length(@$seq[0]);
			$discard_read++;
		}
	}else{
		$sequence_length{$axis} = length(@$seq[0]) unless ($sequence_length{$axis});
		&get_base_coordinate($axis, $cut{$axis}) if ($cut{$axis} && !$rightcutmin{$axis} && !$leftcutmax{$axis});
		if (defined $options{left_position}){
			$leftcutmax{$axis}=$options{left_position}-1 if (!$leftcutmax{$axis} || ($leftcutmax{$axis} && $leftcutmax{$axis}<$options{left_position}-1));
		}
		if (defined $options{right_position}){
			$rightcutmin{$axis}=$options{right_position}+1 if (!$rightcutmin{$axis} || ($rightcutmin{$axis} && $rightcutmin{$axis}>$options{right_position}+1));
		}
		if ($options{whole} || length(@$seq[0])!=$sequence_length{$axis} || (!$rightcutmin{$axis} && !$leftcutmax{$axis})){
			print DISCARD $title."\n".@$seq[0]."\n".@$seq[1]."\n".@$seq[2]."\n" if (defined $options{discard_file});
			$trim_base += length(@$seq[0]);
			$trim_read++;
			$discard_base += length(@$seq[0]);
			$discard_read++;
		}else{
			my ($second_line, $forth_line);
			if ($rightcutmin{$axis} && !$leftcutmax{$axis}){
				$second_line = substr(@$seq[0], 0, $rightcutmin{$axis}-1);
				$forth_line = substr(@$seq[2], 0, $rightcutmin{$axis}-1);
			}elsif (!$rightcutmin{$axis} && $leftcutmax{$axis}){
				$second_line = substr(@$seq[0], $leftcutmax{$axis}, $sequence_length{$axis}-$leftcutmax{$axis});
				$forth_line = substr(@$seq[2], $leftcutmax{$axis}, $sequence_length{$axis}-$leftcutmax{$axis});
			}elsif ($rightcutmin{$axis} && $leftcutmax{$axis}){
				$second_line = substr(@$seq[0], $leftcutmax{$axis}, $rightcutmin{$axis}-1-$leftcutmax{$axis});
				$forth_line = substr(@$seq[2], $leftcutmax{$axis}, $rightcutmin{$axis}-1-$leftcutmax{$axis});
			}
			if ($second_line && length($second_line)>=$options{min_length}){
				if (defined $options{output_file}){
					print OUT $title."\n".$second_line."\n".@$seq[1]."\n".$forth_line."\n";
				}else{
					print STDOUT $title."\n".$second_line."\n".@$seq[1]."\n".$forth_line."\n";
				}
				$remain_base += length($second_line);
				$remain_read++;
				$trim_base += $sequence_length{$axis}-length($second_line);
				$trim_read++;
				$discard_base += $sequence_length{$axis}-length($second_line);
			}elsif(!$second_line){
				print DISCARD $title."\n".@$seq[0]."\n".@$seq[1]."\n".@$seq[2]."\n" if (defined $options{discard_file});
				$trim_base += length(@$seq[0]);
				$trim_read++;
				$discard_base += length(@$seq[0]);
				$discard_read++;
			}elsif($second_line && length($second_line)<$options{min_length}){
				print DISCARD $title."\n".@$seq[0]."\n".@$seq[1]."\n".@$seq[2]."\n" if (defined $options{discard_file});
				$trim_base += $sequence_length{$axis}-length($second_line);
				$trim_read++;
				$filter_base += length($second_line);
				$filter_read++;
				$discard_base += length(@$seq[0]);
				$discard_read++;	
			}
		}
	}
}
