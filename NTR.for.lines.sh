#!/bin/sh
grep -v "^@" $1 | cut -f1-6 >${2}.summary.NTR.tsv
for i in `seq 7 506`;do
	j=$(($i+500))
	k=$(($i+1000))
	l=$(($i+1500))
	m=$(($i+2000))
	n=$(($i+2500))
	grep -v "^@" $1|cut -f$i,$j,$k,$l,$m,$n >tmp.tsv
	Rscript /rd1/user/liym/nucleosomeTurnover/scripts/NTR.R -i=tmp.tsv -p=0.0001 -s=1 -e=6 -o=tmp.NTR.tsv
	cut -f7 tmp.NTR.tsv >NTR.tsv
	paste ${2}.summary.NTR.tsv NTR.tsv >tmp
	mv tmp ${2}.summary.NTR.tsv
done
