#!/bin/sh
usage(){
  cat << EOF
Description: Perform IDR anlysis for ChIP-seq replicates and get final peak list.
Options:
    f   FILE    The peak file of first replicate(MACS2 result).
    s   FILE    The peak file of the second replicate(MACS2 result).
    p   FLOAT   The IDR cutoff for the final peak list.
    c   FILE    The chromsome size file.
    o   STRING  The output prefix for all the result files.
    h           Print this help information.
EOF
    exit 0
}
[ $1 ] || usage

while getopts "hf:s:p:c:o:" OPTION
do
    case $OPTION in
        h) usage;;
        f) peak1=$OPTARG;;
        s) peak2=$OPTARG;;
        p) idr=$OPTARG;;
        c) chrSzie=$OPTARG;;
        o) out=$OPTARG;;
        ?) usage;;
    esac
done

Rscript /mnt/share/liym/tools/idrCode/batch-consistency-analysis.r $peak1 $peak2 -1 $out 0 F p.value $chrSzie >${out}.idr.log 2>${out}.idr.err
Rscript /mnt/share/liym/tools/idrCode/batch-consistency-plot.r 1 idr $out

sed 's/"//g' ${out}-overlapped-peaks.txt |sed 's/ /\t/g' >tmp;mv tmp ${out}-overlapped-peaks.txt
grep -v "IDR" ${out}-overlapped-peaks.txt |awk -v OFS="\t" -v ORS="" '{print $2"\t";if($3<=$7){print $3"\t"}else{print $7"\t"}if($4<=$8){print $8"\n"}else{print $4"\n"}}' >${out}.merged.peaks.bed3 
grep -v "IDR" ${out}-overlapped-peaks.txt |awk -v value=$idr '$11<value'|awk -v OFS="\t" -v ORS="" '{print $2"\t";if($3<=$7){print $3"\t"}else{print $7"\t"}if($4<=$8){print $8"\n"}else{print $4"\n"}}' >${out}.IDR${idr}.merged.peaks.bed3
