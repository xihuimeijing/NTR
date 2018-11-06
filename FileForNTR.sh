#!/bin/sh

if [ $1 == "-h" ] || [ $# -eq 0 ];then
    echo "This script will output the input file for NTR.R
Usage: sh $0 <IN.bed> week0.bw week1.bw week2.bw week4.bw week6.bw week8.bw output.tsv"
exit 0
fi

cut -f1-3 $1 | bwtool summary stdin $2 week0.tmp -skip-median
cut -f1-3 $1 | bwtool summary stdin $3 week1.tmp -skip-median
cut -f1-3 $1 | bwtool summary stdin $4 week2.tmp -skip-median
cut -f1-3 $1 | bwtool summary stdin $5 week4.tmp -skip-median
cut -f1-3 $1 | bwtool summary stdin $6 week6.tmp -skip-median
cut -f1-3 $1 | bwtool summary stdin $7 week8.tmp -skip-median

cut -f1-3,8 week0.tmp >tmp;mv tmp week0.tmp
cut -f8 week1.tmp>tmp;mv tmp week1.tmp
cut -f8 week2.tmp>tmp;mv tmp week2.tmp
cut -f8 week4.tmp>tmp;mv tmp week4.tmp
cut -f8 week6.tmp>tmp;mv tmp week6.tmp
cut -f8 week8.tmp>tmp;mv tmp week8.tmp
paste week0.tmp week1.tmp week2.tmp week4.tmp week6.tmp week8.tmp >$8
rm *.tmp