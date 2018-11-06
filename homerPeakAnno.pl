#!/bin/usr/perl -w
use strict;
use 5.010;
use Getopt::Long;
my ($inBed,$chrSize,$gtf,$outPrefix);
GetOptions(
            'i|bed=s'		=> \$inBed,
			'g|gtf=s'		=> \$gtf,
			'c|chr=s'		=> \$chrSize,     
			'o|out=s'		=> \$outPrefix,
            'h|help'		=> sub{usage()}
	  ) || usage();
my $outTsv="$outPrefix".".homerAnno.tsv";
my $outLog="$outPrefix".".homerAnno.log";
my $outEnrich="$outPrefix".".enrichment.tsv";
`perl /mnt/share/liym/tools/homer/bin/annotatePeaks.pl $inBed none -gtf $gtf >$outTsv 2>$outLog`;
`grep -B 6 "Counting Tags" $outLog|head -n6 |sed 's/^\t\t//'|cut -f1,2 >R.input.tsv`;
#Peak enrichment test, using binomial test to calculate p value
open IN,"$gtf" or die "Can't open file $gtf:$!";
open CHR,"$chrSize" or die "Can't open file $chrSize:$!";
open EXON,">exon_tmp" or die "Can't open file exon_tmp:$!";
open TSS,">TSS_tmp" or die "Can't open file exon_tmp:$!";
open TTS,">TTS_tmp" or die "Can't open file exon_tmp:$!";
open GENE,">gene_tmp" or die "Can't open file gene_tmp:$!";
`perl /mnt/share/liym/bin/gtfFeatureTobed6.pl -i $gtf|grep "intron"|cut -f1-6 |sort -k1,1 -k2,2n|bedtools merge -i stdin >intron_tmp`;
my $number=0;
while(<IN>){
	chomp;
	my @fields=split /\t/;
	if($fields[2] eq "exon"){
		$number+=1;
		say EXON join "\t",$fields[0],$fields[3]-1,$fields[4],$number,"0",$fields[6];
	}
	if($fields[2] eq "transcript"){
		if($fields[3]-1001<0){
			say GENE join "\t",$fields[0],0,$fields[4]+1000,$number,"0",$fields[6]
		}else{
			say GENE join "\t",$fields[0],$fields[3]-1001,$fields[4]+1000,$number,"0",$fields[6];
		}
		if($fields[6] eq "+"){
			if($fields[3]-1001<0){
				say TSS join "\t",$fields[0],0,$fields[3]+100,$number,"0",$fields[6];
			}else{
				say TSS join "\t",$fields[0],$fields[3]-1001,$fields[3]+100,$number,"0",$fields[6];
			}
			if($fields[4]-101<0){
				say TTS join "\t",$fields[0],0,$fields[4]+1000,$number,"0",$fields[6];
			}else{
				say TTS join "\t",$fields[0],$fields[4]-101,$fields[4]+1000,$number,"0",$fields[6];
			}
		}else{
			if($fields[4]-101<0){
				say TSS join "\t",$fields[0],0,$fields[4]+1000,$number,"0",$fields[6];
			}else{
				say TSS join "\t",$fields[0],$fields[4]-101,$fields[4]+1000,$number,"0",$fields[6];
			}
			if($fields[3]-1001<0){
				say TSS join "\t",$fields[0],0,$fields[3]+100,$number,"0",$fields[6];
			}else{
				say TSS join "\t",$fields[0],$fields[3]-1001,$fields[3]+100,$number,"0",$fields[6];
			}
		}
	}
}
`sort -k1,1 -k2,2n exon_tmp|bedtools merge -i stdin >tmp;mv tmp exon_tmp`;
`sort -k1,1 -k2,2n TSS_tmp|bedtools merge -i stdin >tmp;mv tmp TSS_tmp`;
`sort -k1,1 -k2,2n TTS_tmp|bedtools merge -i stdin >tmp;mv tmp TTS_tmp`;
`cat exon_tmp TSS_tmp TTS_tmp|bedtools subtract -a intron_tmp -b stdin |bedtools merge -i stdin >tmp;mv tmp intron_tmp`;
`cat TSS_tmp TTS_tmp|bedtools subtract -a exon_tmp -b stdin |bedtools merge -i stdin >tmp;mv tmp exon_tmp`;
`bedtools subtract -a TTS_tmp -b TSS_tmp |bedtools merge -i stdin >tmp;mv tmp TTS_tmp`;
my $exonSum=`awk 'BEGIN{sum=0}{sum+=\$3-\$2}END{print sum}' exon_tmp`;
my $intronSum=`awk 'BEGIN{sum=0}{sum+=\$3-\$2}END{print sum}' intron_tmp`;
my $TSSSum=`awk 'BEGIN{sum=0}{sum+=\$3-\$2}END{print sum}' TSS_tmp`;
my $TTSSum=`awk 'BEGIN{sum=0}{sum+=\$3-\$2}END{print sum}' TTS_tmp`;
chomp(my $geneSum=`sort -k1,1 -k2,2n gene_tmp|bedtools merge -i stdin|awk 'BEGIN{sum=0}{sum+=\$3-\$2}END{print sum}'`);
my $total=0;
while(<CHR>){
	chomp;
	my @split=split /\t/;
	$total+=$split[1];
}
my $interSum=$total-$geneSum;
open TMP,">tmp" or die "Can't open file tmp:$!";
print TMP join "","TotalSize\n",$TTSSum,$exonSum,$intronSum,"$interSum\n",$TSSSum;
open SUM,"paste R.input.tsv tmp|" or die "Can't open file: $!";
open RST,">$outEnrich" or die "Can't open file R.input.tsv:$!";
#Calculate the hypothesized probability
chomp(my $totalPeaks=`wc -l $inBed|cut -f1 -d " "`);
my $prob=$totalPeaks/$total;
while(<SUM>){
	chomp;
	my @split=split /\t/;
	if(/TotalSize/){say RST join "\t",$_,"Enrichment","P value";next;}
	my $enrich=($split[1]/$split[2])/$prob;
	my $pvalue;
	if($enrich>1){
		$pvalue=`Rscript /mnt/share/liym/bin/binomTest.R -x=$split[1] -n=$split[2] -p=$prob -a="greater"|cut -f2 -d " "`;
	}else{
		$pvalue=`Rscript /mnt/share/liym/bin/binomTest.R -x=$split[1] -n=$split[2] -p=$prob -a="less"|cut -f2 -d " "`;
	}
	print RST join "\t",$_,$enrich,$pvalue;
}
`Rscript /mnt/share/liym/bin/homerPeakAnnoPie.R R.input.final.tsv ${outPrefix}.pdf`;
`rm R.input.tsv exon_tmp intron_tmp TSS_tmp TTS_tmp tmp gene_tmp`;
close;
sub usage{
print <<HELP;
Usage:perl $0 -i <IN.bed> -f <IN.fa> -g <IN.gtf> -s <mouse> -o <peakAnno>
Author:Yumei Li,2017-5-31
Revision: Yumei Li, 2018-2-27, adding peak enrichment test.
Description:This script will annotate the input bed file by homer and plot pie chart and enrichment barplot.
The enrichment is calculated using the binomial test to compare relative peak enrichment in genomic regions of different categories (Refer to PAVIS, https://doi.org/10.1093/bioinformatics/btt520)
Note: The script will generate tmp files R.input.tsv, exon_tmp, intron_tmp, TSS_tmp, TTS_tmp and tmp. Please ensure there is no files named as these tmp files in your current directory.
Options:
    -i|--bed   FILE    The input bed file.
    -g|--gtf   FILE    The gene structure file in GTF format.
	-c|--chr   FILE    The chromsome size file with two columns.
    -o|--out   STRING  The output prefix for output files.
    -h|--help          Print this help information. 
HELP
    exit(-1);
}