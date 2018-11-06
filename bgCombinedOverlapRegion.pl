#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
my($hg19,$otherSpecies,$column);
my $opt=GetOptions(
                        'c|col=i'    => \$column,
                        'hg|human=s' => \$hg19,
                        'o|other=s'  => \$otherSpecies,
                        'h|help'     => sub{&usage;exit(-1);}
                  );

sub usage{
print STDERR <<HELP 
Usage:	perl $0 -hg exon_5juncReads.bed6+ -o rheMac2.tag.final.bed12+,tupBel1.tag.final.bed12+,mm9.tag.final.bed12+ >*.rst 
Output: Output the evolutionary status for the input exons.
        'c|col'    INT      The exon tag column in -hg file which is the same as -o files.[default:4] 
        'hg|human' FILE     The human exons in bed6+ format.
        'o|other'  FILE     *.tag.final.bed12+ file for other species seperated by comma. 
        'help|h'           Print this help message    
HELP
}