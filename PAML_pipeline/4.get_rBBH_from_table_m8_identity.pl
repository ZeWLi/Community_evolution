#!perl -w

use warnings;
use strict;
use List::Util qw(max);


open (M8,"$ARGV[0]") or die "can't open file $ARGV[0]";
open (RBBH,">$ARGV[1]") or die "can't create file $ARGV[1]";

my %similarity = ();
while(my $line = <M8>){
	$line =~ s/\s+$//ig;
	my @fields = split /\t/,$line;
	my $genome1 = $fields[0];
	my $genome2 = $fields[1];
	$genome1 =~ s/\_scaffold.*$//ig;
	$genome2 =~ s/\_scaffold.*$//ig;
        my @array = ($fields[12],$fields[13]);
        my $max = max @array;
	if(!exists($similarity{$fields[0]}->{$genome2}) && $fields[10]<=1e-5 && $fields[2]>=90 && $fields[3]/$max >=0.5){
		print RBBH "$line\n";
		$similarity{$fields[0]}->{$genome2} = 1;
	}
}
close RBBH;
close M8;
