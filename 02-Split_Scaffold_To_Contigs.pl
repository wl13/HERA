#!/usr/bin/perl
#Author: huilong du
#Note: splicing the scaffolds into contigs by name count++
use warnings;
use strict;
my $infile=shift;
my $outfile=shift;
my $Enzyme=shift;
if(!defined $infile || !defined $outfile || !defined $Enzyme ){
	print "Usage: perl $0 infile outfile Enzyme\n";
	exit;
}
$Enzyme="N".$Enzyme."N";
open IN,"<$infile" or die $!;
open OUT,">$outfile" or die $!;

my $contigs;
my $count=0;
while(<IN>){
	chomp;
	my $line=$_;
        if($line=~/^>(\S+)/){
               $contigs=$1;
               $count=1;
        }
        else{
	       $line=~ s/$Enzyme/N/g;
               $line =~ s/N*N/N/g;
               $line =~ s/n/N/g;
               my @splitline = split(/N/,$line);
               for(my $i=0;$i<@splitline;$i++){
                        print OUT ">$contigs.$count\n";
                        print OUT "$splitline[$i]\n";
                        $count++;
               }
        }
}
close(OUT);
close(IN);
