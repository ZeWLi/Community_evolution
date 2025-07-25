#!perl -w

use warnings;
use strict;

opendir( DIR, "$ARGV[0]" ) or die "can't open file $ARGV[0]";

foreach my $file ( readdir DIR ) {
        if ( $file ne "." && $file ne ".." ) {
                chdir("$ARGV[0]");
                my $out_dir = $file;
                $out_dir =~ s/\.fa|\.fasta|\.fna|\.fsa_nt//ig;
                mkdir ("$out_dir");                
#gene prediction using prodigal
                my $pro_name = $out_dir . ".pro.fa";
                my $gene_name = $out_dir . ".gene.fa";
                system(
"prodigal -a $pro_name -d $gene_name -f gff -g 11 -i $file -p single -q -m"
                );
                system ("mv $pro_name $gene_name $out_dir");
                chdir("../");
        }
}

chdir("$ARGV[0]");
