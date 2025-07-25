#!perl -w

use warnings;
use strict;

open( INPUT,  "<$ARGV[0]" ) or die "cna't open file $ARGV[0]";
open( OUTPUT, ">$ARGV[1]" ) or die "can't create file $ARGV[1]";
my $length = $ARGV[2];

my $contig="";
my $contigLen=0;
my $contigSeq="";
while(my $line=<INPUT>){
	if($line=~/^>.*/){
		if($contigLen>=$length){
			print OUTPUT $contig;
			print OUTPUT $contigSeq;
#			print "$contig	$contigLen\n";
		}
		$contig = $line;
		$contigLen = 0;
		$contigSeq = "";
	}
	else{
		$contigSeq = $contigSeq.$line;
		chomp($line);
		$line=~s/\s+$//ig;
		my @fields = split //,$line,-1;
		$contigLen = $contigLen + length($line);
#		print "1$fields[$#fields]1\n";
#		$contigLen = $contigLen + scalar(@fields);
	}
}
if($contigLen>=$length){
	print OUTPUT $contig;
	print OUTPUT $contigSeq;
#	print "$contig	$contigLen\n";
}
close INPUT;
close OUTPUT;
