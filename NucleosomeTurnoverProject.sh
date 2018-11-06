#!/bin/sh
##path: /rd1/user/liym/NTR @venus
#1. Data processing
#1.1 H2BGFP and H3 ChIP-seq data in WT at six time points
mkdir H2BGFP_All && cd H2BGFP_All
ls -d */|while read dir;do
        cd $dir;
        ls *.fq|sed 's/.fq//'|while read file;do
                mkdir $file 
                cd $file
                mv ../${file}.fq .
                sh /mnt/share/liym/nucleosome/scripts/MNaseSeProcessing.sh 8 50 /mnt/share/liym/data/bwa/mm10/mm10.fa /home/share/data/chr.size/mm10.size ${file}.fq >MNaseSeProcessing.log 2>MNaseSeProcessing.err
                cd ../
        done;
        cd ../;
done;
ls -d *|grep -v "figures"|while read dir;do
        cd $dir;
        ls -d ${dir}_*|while read sample;do 
                cd $sample;
                bamtools filter -in out.sorted.bam -out myFiltered.sorted.bam -script /mnt/share/liym/data/bwa/bwaAln.filter.json 2>log/bamtoolsFilter.log;
                bamtools stats -in myFiltered.sorted.bam >myFiltered.sorted.bamStats & 
                samtools index myFiltered.sorted.bam;
                samtools view -bu myFiltered.sorted.bam chr1 chr2 chrX chr3 chr4 chr5 chr6 chr7 chr10 chr8 chr14 chr9 chr11 chr13 chr12 chr15 chr16 chr17 chrY chr18 chr19 chrM |samtools sort -@ 5 - ${sample}.final 2>log/samtoolsSort.log;
                samtools rmdup -s ${sample}.final.bam >${sample}.final.rmDup.bam 2>rmDup.log 
                bamtools stats -in ${sample}.final.bam >${sample}.final.bamStats & 
                cd ../;
        done;
        cd ../;
done;
ls */danposRst/result/pooled/*.wig |grep -v "Input"|while read file;do prefix=$(echo $file|sed 's/bgsub.Fnor.smooth.wig//');wigToBigWig $file /mnt/share/liym/data/chr.size/mm10.filterRandom.size ${prefix}.bw;done;
#1.2 H2BGFP ChIP-seq data in EEDCKO and WT
mkdir 3 4
mkdir 3/fastqc;ls 3/*/*gz|while read file;do fastqc -q $file -o 3/fastqc;done;
mkdir 4/fastqc;ls 4/*/*gz|while read file;do fastqc -q $file -o 4/fastqc;done;
cd ../;
ls -d */|while read dir;do
        cd $dir;
        ls *.fq|sed 's/.fq//'|while read file;do
                mkdir $file 
                cd $file
                mv ../${file}.fq .
                sh /mnt/share/liym/nucleosome/scripts/MNaseSeProcessing.sh 8 50 /mnt/share/liym/data/bwa/mm10/mm10.fa /home/share/data/chr.size/mm10.size ${file}.fq >MNaseSeProcessing.log 2>MNaseSeProcessing.err
                cd ../
        done;
        cd ../;
done;
ls -d *|grep -v "figures"|while read dir;do
        cd $dir;
        ls -d ${dir}_*|while read sample;do 
                cd $sample;
                bamtools filter -in out.sorted.bam -out myFiltered.sorted.bam -script /mnt/share/liym/data/bwa/bwaAln.filter.json 2>log/bamtoolsFilter.log;
                bamtools stats -in myFiltered.sorted.bam >myFiltered.sorted.bamStats & 
                samtools index myFiltered.sorted.bam;
                samtools view -bu myFiltered.sorted.bam chr1 chr2 chrX chr3 chr4 chr5 chr6 chr7 chr10 chr8 chr14 chr9 chr11 chr13 chr12 chr15 chr16 chr17 chrY chr18 chr19 chrM |samtools sort -@ 5 - ${sample}.final 2>log/samtoolsSort.log;
                samtools rmdup -s ${sample}.final.bam >${sample}.final.rmDup.bam 2>rmDup.log 
                bamtools stats -in ${sample}.final.bam >${sample}.final.bamStats &
                cd ../;
        done;
        cd ../;
done;
ls */*H2BGFP*EED*/*final.rmDup.bam|while read file;do
        prefix=$(echo $file|cut -f2 -d '/');
        dir=$(dirname $file);
        java -jar /mnt/share/share/tool/picard-tools-2.2.4/picard.jar CollectInsertSizeMetrics I=$file O=${dir}/${prefix}_insert_size_metrics.txt H=${dir}/${prefix}_insert_size_histogram.pdf >${dir}/${prefix}.collectInsertS.log 2>${dir}/${prefix}.collectInsertS.err
        Rscript /mnt/share/liym/tools/phantompeakqualtools/run_spp.R -c=$file -s=0:5:1500 -savp=${dir}/${prefix}.estFrag.pdf -out=${dir}/${prefix}.spp.estFrag.rst >${dir}/${prefix}.spp.log 2>${dir}/${prefix}.spp.err;
        java -jar /mnt/share/share/tool/picard-tools-2.2.4/picard.jar CollectGcBiasMetrics I=$file O=${dir}/${prefix}.GCbias.metrix.txt CHART=${dir}/${prefix}.GCbias.metrix.pdf S=${dir}/${prefix}.GCbias.summary.metrix.txt R=/mnt/share/liym/data/bwa/mm10/mm10.fa >${dir}/${prefix}.GCbias.log 2> ${dir}/${prefix}.GCbias.err
		Rscript /mnt/share/liym/bin/picard.GCbiasMetrix.plot.R -i=${dir}/${prefix}.GCbias.metrix.txt -y1=0 -y2=2 -o=${dir}/${prefix}.GCbias.metrix.inhouse.pdf
done;
cd 3/danposRst;
ls ../*/*final.rmDup.bam|while read file;do prefix=$(echo $file|cut -f2 -d '/');ln -s $file ${prefix}.final.bam;done;
python /rd1/user/liym/tools/danpos-2.2.2/danpos.py dpos 3_H2BGFP-CM-week0-EEDheto_Rep1.final.bam,3_H2BGFP-CM-week0-EEDheto_Rep2.final.bam,3_H2BGFP-CM-week0-EEDko_Rep1.final.bam -b 3_Input-CM-week0-EED_Rep1.final.bam -c 10000000 --extend 74 -o ./ >danpos.log 2>danpos.err
cd ../../4/danposRst;
ls ../*/*final.rmDup.bam|while read file;do prefix=$(echo $file|cut -f2 -d '/');ln -s $file ${prefix}.final.bam;done;
python /rd1/user/liym/tools/danpos-2.2.2/danpos.py dpos 4_H2BGFP-CM-week4-EEDheto_Rep1.final.bam,4_H2BGFP-CM-week4-EEDheto_Rep2.final.bam,4_H2BGFP-CM-week4-EEDko_Rep1.final.bam,4_H2BGFP-CM-week4-EEDko_Rep2.final.bam -b 4_Input-CM-week4-EED_Rep1.final.bam -c 10000000 --extend 74 -o ./ >danpos.log 2>danpos.err
cd ../;
ls */danposRst/result/pooled/*.wig |grep -v "Input"|while read file;do prefix=$(echo $file|sed 's/bgsub.Fnor.smooth.wig//');wigToBigWig $file /mnt/share/liym/data/chr.size/mm10.filterRandom.size ${prefix}.bw;done;

#1.3 RNA-seq in EEDCKO and WT
mkdir RNA-seq && cd RNA-seq
##EEDheto_2
fqAdapterFilter.pl -a GAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCT,AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATG -i read1_filter.fq -o read1_adaF_filter.fq --a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCC --i2 read2_filter.fq --o2 read2_adaF_filter.fq -r fqAdaFilter.log &
##EEDko_1
fqAdapterFilter.pl -a GAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCT,AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATG -i read1_filter.fq -o read1_adaF_filter.fq --i2 read2_filter.fq --o2 read2_adaF_filter.fq -r fqAdaFilter.log &
##EEDko_2
fqAdapterFilter.pl -a GAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCT,AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATG,GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGC -i read1_filter.fq -o read1_adaF_filter.fq --a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCC,GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCG,GAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCA --i2 read2_filter.fq --o2 read2_adaF_filter.fq -r fqAdaFilter.log &

ls -d *_*/|grep -vE "EEDheto_1|data"|sed 's;/;;'|while read dir;do
        cd $dir;
        /rd1/user/liym/tools/tophat-2.1.1.Linux_x86_64/tophat2 --read-mismatches 4 --read-gap-length 3 --read-edit-dist 4 --min-anchor 4 --splice-mismatches 0 --num-threads 10 --no-coverage-search --segment-length 25 --segment-mismatches 2 --library-type fr-unstranded /mnt/share/liym/data/bowtie2/mm10/mm10 read1_adaF_filter.fq read2_adaF_filter.fq > tophat2.log 2> tophat2.err;
        samtools view -bu -q 50 -@ 5 tophat_out/accepted_hits.bam | samtools sort -@ 5 -m 10G -o ${dir}.uniq.sorted.bam - 2>samtoolsSort.log
        cd ../;
        done;
ls -d */|grep -v "data"|sed 's;/;;'|while read dir;do cd $dir;cufflinks -p 5 -o cufflinks -G /rd1/user/liym/nucleosomeTurnover/data/mm10.refGene.filterRandom.gtf --max-bundle-frags 10000000 --no-update-check ${dir}.uniq.sorted.bam >cufflinks.log 2>cufflinks.err;cd ../;done;

#1.4 Histone modification ChIP-seq data (public)
mkdir ChIP-seq/histoneModif && cd ChIP-seq/histoneModif
cd EEDproject
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGC -i 1_H3K27me3-EEDko-AdultCM_Rep1.fastq -o 1_H3K27me3-EEDko-AdultCM_Rep1.adaF.fastq -r log/1_H3K27me3-EEDko-AdultCM_Rep1.adaF.log &
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGC -i 1_H3K27me3-EEDko-AdultCM_Rep2.fastq -o 1_H3K27me3-EEDko-AdultCM_Rep2.adaF.fastq -r log/1_H3K27me3-EEDko-AdultCM_Rep2.adaF.log &
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGC -i 1_H3K27me3-WT-AdultCM_Rep1.fastq -o 1_H3K27me3-WT-AdultCM_Rep1.adaF.fastq -r log/1_H3K27me3-WT-AdultCM_Rep1.adaF.log &
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGC -i 1_H3K27me3-WT-AdultCM_Rep2.fastq -o 1_H3K27me3-WT-AdultCM_Rep2.adaF.fastq -r log/1_H3K27me3-WT-AdultCM_Rep2.adaF.log &
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGC -i 1_H3K27me3-WT-AdultCM_Rep3.fastq -o 1_H3K27me3-WT-AdultCM_Rep3.adaF.fastq -r log/1_H3K27me3-WT-AdultCM_Rep3.adaF.log &
ls *.fastq|sed 's/.fastq//'|while read file;do 
        bwa aln -t 10 /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}.fastq >${file}.out.sai 2>log/${file}.aln.log;
        bwa samse /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}.out.sai ${file}.fastq 2>log/${file}.samse.log |samtools view -bSu -| /home/liym/.local/bin/samtools sort -@ 5 -o ${file}.out.sorted.bam - ;
        bamtools filter -in ${file}.out.sorted.bam -out ${file}.uniq.sorted.bam -script /mnt/share/liym/data/bwa/bwaAln.filter.json 2>log/${file}.bamtoolsFilter.log
        samtools rmdup -s ${file}.uniq.sorted.bam ${file}.final.rmDup.bam 2>log/${file}.rmDup.log
done;
ls 1_H3K27ac*final.rmDup.bam|while read file;do
        prefix=$(echo $file|sed 's/.final.rmDup.bam//');
        macs2 callpeak -t $file -c 1_Input-H3K27ac-AdultCM_Rep1.final.rmDup.bam -f BAM -g mm --broad --outdir $prefix -n $prefix -B --SPMR >log/${prefix}.log 2>log/${prefix}.err
done;
ls 1_H3K27me3*final.rmDup.bam|while read file;do
        prefix=$(echo $file|sed 's/.final.rmDup.bam//');
        macs2 callpeak -t $file -c 1_Input-AdultCM_Rep1.final.rmDup.bam -f BAM -g mm --broad --outdir ${prefix}_broad -n $prefix -B --SPMR >log/${prefix}_broad.log 2>log/${prefix}_broad.err
done;
ls *.final.rmDup.bam|while read file;do samtools index $file;done;
ls *H3K27me3*final.rmDup.bam|while read file;do
    prefix=$(echo $file|sed 's/.final.rmDup.bam//')
    input=$(grep "Total" 1_Input-AdultCM_Rep1.final.rmDup.bamStats|awk 'BEGIN{sum=0}{sum+=$3}END{print 10000000/sum}')
    factor=$(grep "Total" ${file}Stats|awk 'BEGIN{sum=0}{sum+=$3}END{print 10000000/sum}')
    bamCompare -b1 $file -b2 1_Input-AdultCM_Rep1.final.rmDup.bam --scaleFactors ${factor}:$input --ratio subtract -bs 20 -p 15 --extendReads 140 -o ${prefix}.subtract.bw --outFileFormat bigwig >log/bamCompare.log 2>log/bamCompare.err
done;
ls *H3K27ac*final.rmDup.bam|while read file;do
        prefix=$(echo $file|sed 's/.final.rmDup.bam//')
        input=$(grep "Total" 1_Input-H3K27ac-AdultCM_Rep1.final.rmDup.bamStats|awk 'BEGIN{sum=0}{sum+=$3}END{print 10000000/sum}')
        factor=$(grep "Total" ${file}Stats|awk 'BEGIN{sum=0}{sum+=$3}END{print 10000000/sum}')
        bamCompare -b1 $file -b2 1_Input-H3K27ac-AdultCM_Rep1.final.rmDup.bam --scaleFactors ${factor}:$input --ratio subtract -bs 20 -p 15 --extendReads 140 -o ${prefix}.subtract.bw --outFileFormat bigwig >log/bamCompare.log 2>log/bamCompare.err
done;
mkdir result && cd result 
bedtools intersect -a ../1_H3K27ac-EEDko-AdultCM_Rep1/1_H3K27ac-EEDko-AdultCM_Rep1_peaks.broadPeak -b ../1_H3K27ac-EEDko-AdultCM_Rep2/1_H3K27ac-EEDko-AdultCM_Rep2_peaks.broadPeak -wo |awk -v OFS="\t" -v ORS="" '{print $1"\t";if($2<=$11){print $2"\t";}else{print $11"\t";}if($3<=$12){print $12"\n"}else{print $3"\n"}}'|sort|uniq |sort -k1,1 -k2,2n|bedtools merge >H3K27ac-EEDko.intersect.broadPeaks.bed 
bedtools intersect -a ../1_H3K27ac-WT-AdultCM_Rep1/1_H3K27ac-WT-AdultCM_Rep1_peaks.broadPeak -b ../1_H3K27ac-WT-AdultCM_Rep2/1_H3K27ac-WT-AdultCM_Rep2_peaks.broadPeak -wo |awk -v OFS="\t" -v ORS="" '{print $1"\t";if($2<=$11){print $2"\t";}else{print $11"\t";}if($3<=$12){print $12"\n"}else{print $3"\n"}}'|sort|uniq |sort -k1,1 -k2,2n|bedtools merge >H3K27ac-WT.intersect.broadPeaks.bed
bedtools intersect -a ../1_H3K27me3-EEDko-AdultCM_Rep1.adaF_broad/1_H3K27me3-EEDko-AdultCM_Rep1.adaF_peaks.broadPeak -b ../1_H3K27me3-EEDko-AdultCM_Rep2.adaF_broad/1_H3K27me3-EEDko-AdultCM_Rep2.adaF_peaks.broadPeak -wo |awk -v OFS="\t" -v ORS="" '{print $1"\t";if($2<=$11){print $2"\t";}else{print $11"\t";}if($3<=$12){print $12"\n"}else{print $3"\n"}}'|sort|uniq |sort -k1,1 -k2,2n|bedtools merge >H3K27me3-EEDko.intersect.broadPeaks.bed
bedtools intersect -wo -a ../1_H3K27me3-WT-AdultCM_Rep1.adaF_broad/1_H3K27me3-WT-AdultCM_Rep1.adaF_peaks.broadPeak -b ../1_H3K27me3-WT-AdultCM_Rep2.adaF_broad/1_H3K27me3-WT-AdultCM_Rep2.adaF_peaks.broadPeak|awk -v OFS="\t" -v ORS="" '{print $1"\t";if($2<=$11){print $2"\t";}else{print $11"\t";}if($3<=$12){print $12"\n"}else{print $3"\n"}}'|sort|uniq|sort -k1,1 -k2,2n |bedtools merge|bedtools intersect -wo -a stdin -b ../1_H3K27me3-WT-AdultCM_Rep3.adaF_broad/1_H3K27me3-WT-AdultCM_Rep3.adaF_peaks.broadPeak|awk -v OFS="\t" -v ORS="" '{print $1"\t";if($2<=$5){print $2"\t";}else{print $5"\t";}if($3<=$6){print $6"\n"}else{print $3"\n"}}'|sort -k1,1 -k2,2n|uniq |bedtools merge >H3K27me3-WT.intersect.broadPeaks.bed
ls *broadPeaks.bed |while read file;do
    awk '$1!~"random|GL|Un|chrM"' $file >tmp;
    mv tmp $file;
done;
cd ENCODE
#H3K4me1
wget https://www.encodeproject.org/files/ENCFF481GRM/download/ENCFF481GRM.bed.gz
gunzip -c ENCFF481GRM.bed.gz >H3K4me1.ChIP.8w.mouseHeart.narrowPeak.bed
#H3K4me3
wget https://www.encodeproject.org/files/ENCFF599BFW/download/ENCFF599BFW.bed.gz
gunzip -c ENCFF599BFW.bed.gz >H3K4me3.ChIP.8w.mouseHeart.narrowPeak.bed

#1.5 TF ChIP-seq data (public)
mkdir ChIP-seq/TFs && cd ChIP-seq/TFs
cd GSE29636_JCI2011_SE36 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR223/SRR223492/SRR223492.fastq.gz
wget -q -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR223/SRR223493/SRR223493.fastq.gz
wget -q -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR223/SRR223494/SRR223494.fastq.gz
wget -q -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR223/SRR223495/SRR223495.fastq.gz
fqSeFilter.pl -r 30 1_Input-SRR223494-wholeHeart-2months_Rep2.fastq > 1_Input-SRR223494-wholeHeart-2months_Rep2.filter.fastq 2> fqFilter_input_Rep2.log
mkdir log
ls *.fastq|sed 's/.fastq//'|while read file;do 
    bwa aln -t 3 /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}.fastq >${file}.out.sai 2>log/${file}.aln.log;
    bwa samse /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}.out.sai ${file}.fastq 2>log/${file}.samse.log |samtools view -bSu -| /home/liym/.local/bin/samtools sort -@ 3 -o ${file}.out.sorted.bam - ;
    bamtools filter -in ${file}.out.sorted.bam -out ${file}.uniq.sorted.bam -script /mnt/share/liym/data/bwa/bwaAln.filter.json 2>log/${file}.bamtoolsFilter.log
        samtools rmdup -s ${file}.uniq.sorted.bam ${file}.final.rmDup.bam 2>log/${file}.rmDup.log
done;
macs2 callpeak -t 1_SRR223492-Tbx20-wholeHeart-2months_Rep1.final.rmDup.bam -c 1_Input-SRR223493-wholeHeart-2months_Rep1.final.rmDup.bam -f BAM -g mm --outdir Tbx20_rep1 -n Tbx20 >log/1_SRR223492-Tbx20-wholeHeart-2months_Rep1.macs2.log 2>log/1_SRR223492-Tbx20-wholeHeart-2months_Rep1.macs2.err
macs2 callpeak -t 1_SRR223492-Tbx20-wholeHeart-2months_Rep1.final.rmDup.bam -c 1_Input-SRR223495-wholeHeart-2months_Rep3.final.rmDup.bam -f BAM -g mm --outdir Tbx20_rep3 -n Tbx20 >log/1_SRR223492-Tbx20-wholeHeart-2months_Rep3.macs2.log 2>log/1_SRR223492-Tbx20-wholeHeart-2months_Rep3.macs2.err
macs2 callpeak -t 1_SRR223492-Tbx20-wholeHeart-2months_Rep1.final.rmDup.bam -c 1_Input-SRR223494-wholeHeart-2months_Rep2.filter.final.rmDup.bam -f BAM -g mm --outdir Tbx20_rep2 -n Tbx20 >log/1_SRR223492-Tbx20-wholeHeart-2months_Rep2.macs2.log 2>log/1_SRR223492-Tbx20-wholeHeart-2months_Rep2.macs2.err
cat <(cut -f1-3 Tbx20_rep1/Tbx20_peaks.narrowPeak) <(cut -f1-3 Tbx20_rep2/Tbx20_peaks.narrowPeak) <(cut -f1-3 Tbx20_rep3/Tbx20_peaks.narrowPeak)|sort -k1,1 -k2,2n |bedtools merge -i stdin >Tbx20_combined.peaks.bed3
Rscript /mnt/share/liym/tools/phantompeakqualtools/run_spp.R -c=input.merged.final.bam -s=0:5:1500 -savp=input.estFrag.pdf -out=input.spp.estFrag.rst
Rscript /mnt/share/liym/tools/phantompeakqualtools/run_spp.R -c=1_SRR223492-Tbx20-wholeHeart-2months_Rep1.final.rmDup.bam -s=0:5:1500 -savp=Tbx20.estFrag.pdf -out=Tbx20.spp.estFrag.rst
samtools merge -@ 10 input.merged.final.bam 1_Input-SRR22349*final.rmDup.bam 
input=$(grep "Total" input.merged.final.bamStats|awk 'BEGIN{sum=0}{sum+=$3}END{print 10000000/sum}')
factor=$(grep "Total" 1_SRR223492-Tbx20-wholeHeart-2months_Rep1.final.rmDup.bamStats|awk 'BEGIN{sum=0}{sum+=$3}END{print 10000000/sum}')
samtools index 1_SRR223492-Tbx20-wholeHeart-2months_Rep1.final.rmDup.bam &
samtools index input.merged.final.bam &
/mnt/share/liym/tools/deepTools-2.0.0/bin/bamCompare -b1 1_SRR223492-Tbx20-wholeHeart-2months_Rep1.final.rmDup.bam -b2 input.merged.final.bam --scaleFactors ${factor}:$input --ratio subtract -bs 20 -p 10 --extendReads 70 -o Tbx20.subtract.bw --outFileFormat bigwig >log/bamCompare.log 2>log/bamCompare.err
cd GSE35151_adultHeart
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR400/SRR400044/SRR400044.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR400/SRR400043/SRR400043.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR400/SRR400046/SRR400046.fastq.gz
fqTrimer.pl -r [ATGCN]{14} 1_Input-SRR400044-Heart-control_Rep1.fastq >1_Input-SRR400044-Heart-control_Rep1.trim.fq
fqSeFilter.pl -f 0.1 -b N -l 30 -r 30 -q 0.5 -c 20 -a 20 1_Input-SRR400044-Heart-control_Rep1.trim.fq > 1_Input-SRR400044-Heart-control_Rep1.filter.fq 2> input.fqFilter.log
fqAdapterFilter.pl -a ATCGGAAGAGCTCGTATGCCGTCTTCTGCTTAGAT,ATCGGAAGAGCTCGTATGCCGTCTTCTGCTTGAAA,ATCGGAAGAGCTCGTATGCCGTCTTCTGCTTATAT -i 1_Input-SRR400044-Heart-control_Rep1.filter.fq -o 1_Input-SRR400044-Heart-control_Rep1.filter.adaF.fastq -r 1_Input-SRR400044-Heart-control_Rep1.filter.adaF.log
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTGAACATCTCGTATGCC,GGTGTTGTTGTTGTCTTAGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG,GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTAGATCGGAAGAGCTCGTAT -i 1_SRR400046-Heart-Nkx2-5_Rep1.fastq -o 1_SRR400046-Heart-Nkx2-5_Rep1.filter.fastq -r 1_SRR400046-Heart-Nkx2-5_Rep1.adaFilter.log
fqTrimer.pl -r [ATGCN]{14} 1_SRR400043-Heart-Tbx3_Rep1.fastq >1_SRR400043-Heart-Tbx3_Rep1.trim.fastq 2>1_SRR400043-Heart-Tbx3_Rep1.fqTrim.log
fqSeFilter.pl -f 0.1 -b N -l 30 -r 30 -q 0.5 -c 20 -a 20 -t 10,5 1_SRR400043-Heart-Tbx3_Rep1.trim.fastq > 1_SRR400043-Heart-Tbx3_Rep1.filter.fastq 2>1_SRR400043-Heart-Tbx3_Rep1.fqFilter.log
mkdir log
ls *filter*fastq|sed 's/.fastq//'|while read file;do
        bwa aln -t 5 /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}.fastq >${file}.out.sai 2>log/${file}.aln.log;
        bwa samse /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}.out.sai ${file}.fastq 2>log/${file}.samse.log |samtools view -bSu -| /home/liym/.local/bin/samtools sort -@ 3 -o ${file}.out.sorted.bam - ;
        bamtools filter -in ${file}.out.sorted.bam -out ${file}.uniq.sorted.bam -script /mnt/share/liym/data/bwa/bwaAln.filter.json 2>log/${file}.bamtoolsFilter.log
        samtools rmdup -s ${file}.uniq.sorted.bam ${file}.final.rmDup.bam 2>log/${file}.rmDup.log
        bamtools stats -in ${file}.final.rmDup.bam >${file}.final.rmDup.bamStats
done;
macs2 callpeak -t 1_SRR400043-Heart-Tbx3_Rep1.filter.final.rmDup.bam -c 1_Input-SRR400044-Heart-control_Rep1.filter.adaF.final.rmDup.bam -f BAM -g mm --outdir Tbx3 -n Tbx3 -B --SPMR >log/1_SRR400043-Heart-Tbx3_Rep1.macs2.log 2>log/1_SRR400043-Heart-Tbx3_Rep1.macs2.err
macs2 callpeak -t 1_SRR400046-Heart-Nkx2-5_Rep1.filter.final.rmDup.bam -c 1_Input-SRR400044-Heart-control_Rep1.filter.adaF.final.rmDup.bam -f BAM -g mm --outdir Nkx2-5 -n Nkx2-5 -B --SPMR >log/1_SRR400046-Heart-Nkx2-5_Rep1.macs2.log 2>log/1_SRR400046-Heart-Nkx2-5_Rep1.macs2.err
        Rscript /mnt/share/liym/tools/phantompeakqualtools/run_spp.R -c=1_SRR400043-Heart-Tbx3_Rep1.filter.final.rmDup.bam -s=0:5:1500 -savp=Tbx3.estFrag.pdf -out=Tbx3.spp.estFrag.rst > log/Tbx3.spp.log 2> log/Tbx3.spp.err
        Rscript /mnt/share/liym/tools/phantompeakqualtools/run_spp.R -c=1_SRR400046-Heart-Nkx2-5_Rep1.filter.final.rmDup.bam -s=0:5:1500 -savp=Nkx2-5.estFrag.pdf -out=Nkx2-5.spp.estFrag.rst > log/Nkx2-5.spp.log 2> log/Nkx2-5.spp.err 
input=$(grep "Total" 1_Input-SRR400044-Heart-control_Rep1.filter.adaF.final.rmDup.bamStats|awk '{print 10000000/$3}')
tbx3=$(grep "Total" 1_SRR400043-Heart-Tbx3_Rep1.filter.final.rmDup.bamStats|awk '{print 10000000/$3}')
nkx25=$(grep "Total" 1_SRR400046-Heart-Nkx2-5_Rep1.filter.final.rmDup.bamStats|awk '{print 10000000/$3}')
/mnt/share/liym/tools/deepTools-2.0.0/bin/bamCompare -b1 1_SRR400043-Heart-Tbx3_Rep1.filter.final.rmDup.bam -b2 1_Input-SRR400044-Heart-control_Rep1.filter.adaF.final.rmDup.bam --scaleFactors ${tbx3}:$input --ratio subtract -bs 20 -p 10 --extendReads 109 -o Tbx3.subtract.bw --outFileFormat bigwig >log/bamCompare.log 2>log/bamCompare.err
/mnt/share/liym/tools/deepTools-2.0.0/bin/bamCompare -b1 1_SRR400046-Heart-Nkx2-5_Rep1.filter.final.rmDup.bam -b2 1_Input-SRR400044-Heart-control_Rep1.filter.adaF.final.rmDup.bam --scaleFactors ${nkx25}:$input --ratio subtract -bs 20 -p 10 --extendReads 97 -o Nkx2-5.subtract.bw --outFileFormat bigwig >log/bamCompare.log 2>log/bamCompare.err
cd GSE52123-2014NC
perl /mnt/share/liym/bin/fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGC -i 1_Pol2-WT-VentricularApex-Adult-CM_Rep1.fastq -o 1_Pol2-WT-VentricularApex-Adult-CM_Rep1.filter.fastq -r 1_Pol2-WT-VentricularApex-Adult-CM_Rep1.fqAdaFilter.log
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGC,GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAAAATCTCGTATGC,CGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGT -i 1_Pol2-WT-VentricularApex-Adult-CM_Rep2.fastq -o 1_Pol2-WT-VentricularApex-Adult-CM_Rep2.filter.fastq -r 1_Pol2-WT-VentricularApex-Adult-CM_Rep2.fqAdaFilter.log
ls *.fastq|sed 's/.fastq//'|while read file;do 
    bwa aln -t 5 /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}.fastq >${file}.out.sai 2>log/${file}.aln.log;
    bwa samse /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}.out.sai ${file}.fastq 2>log/${file}.samse.log |samtools view -bSu -| /home/liym/.local/bin/samtools sort -@ 5 -o ${file}.out.sorted.bam - ;
        bamtools filter -in ${file}.out.sorted.bam -out ${file}.uniq.sorted.bam -script /mnt/share/liym/data/bwa/bwaAln.filter.json 2>log/${file}.bamtoolsFilter.log
        samtools rmdup -s ${file}.uniq.sorted.bam ${file}.final.rmDup.bam 2>log/${file}.rmDup.log
done;
macs2 callpeak -t 1_Pol2-WT-VentricularApex-Adult-CM_Rep2.filter.final.rmDup.bam -c 1_Input-Adult-histone-VentricularApex-CM_Rep1.final.rmDup.bam -f BAM -g mm --outdir Pol2_Rep2 -n Pol2 -B --SPMR >log/1_Pol2-WT-VentricularApex-Adult-CM.Rep2.macs2.log 2>log/1_Pol2-WT-VentricularApex-Adult-CM.Rep2.macs2.err &
macs2 callpeak -t 1_Pol2-WT-VentricularApex-Adult-CM_Rep1.filter.final.rmDup.bam -c 1_Input-Adult-histone-VentricularApex-CM_Rep1.final.rmDup.bam -f BAM -g mm --outdir Pol2_Rep1 -n Pol2 -B --SPMR >log/1_Pol2-WT-VentricularApex-Adult-CM.Rep1.macs2.log 2>log/1_Pol2-WT-VentricularApex-Adult-CM.Rep1.macs2.err &
macs2 callpeak -t 2_GATA4-fb-Adult-SRR1025222-SRR1025223_Rep2.final.rmDup.bam -c 2_Input-fb-Ab-Adult-SRR1025224_Rep1.final.rmDup.bam -f BAM -g mm --outdir Gata4_Rep2 -n Gata4 -B --SPMR >log/2_GATA4-fb-Adult-SRR1025222-SRR1025223.Rep2.macs2.log 2>log/2_GATA4-fb-Adult-SRR1025222-SRR1025223.Rep2.macs2.err &
macs2 callpeak -t 2_GATA4-fb-Adult-SRR1025222-SRR1025223_Rep1.final.rmDup.bam -c 2_Input-fb-Ab-Adult-SRR1025224_Rep1.final.rmDup.bam -f BAM -g mm --outdir Gata4_Rep1 -n Gata4 -B --SPMR >log/2_GATA4-fb-Adult-SRR1025222-SRR1025223.Rep1.macs2.log 2>log/2_GATA4-fb-Adult-SRR1025222-SRR1025223.Rep1.macs2.err &
ls *final.rmDup.bam|while read file;do
        prefix=$(echo $file|sed 's/.final.rmDup.bam//');
        Rscript /mnt/share/liym/tools/phantompeakqualtools/run_spp.R -c=$file -s=0:5:1500 -savp=${prefix}.estFrag.pdf -out=${prefix}.spp.estFrag.rst > log/${prefix}.spp.log 2> log/${prefix}.spp.err;
done;
input=$(grep "Total" 2_Input-fb-Ab-Adult-SRR1025224_Rep1.final.rmDup.bamStats|awk 'BEGIN{sum=0}{sum+=$3}END{print 10000000/sum}')
factor=$(grep "Total" 2_GATA4-fb-Adult-SRR1025222-SRR1025223_Rep1.final.rmDup.bamStats|awk 'BEGIN{sum=0}{sum+=$3}END{print 10000000/sum}')
bamCompare -b1 2_GATA4-fb-Adult-SRR1025222-SRR1025223_Rep1.final.rmDup.bam -b2 2_Input-fb-Ab-Adult-SRR1025224_Rep1.final.rmDup.bam --scaleFactors ${factor}:$input --operation subtract -bs 20 -p 10 --extendReads 150 -o GATA4_Rep1.subtract.bw --outFileFormat bigwig >log/bamCompare.Gata4_Rep1.log 2>log/bamCompare.Gata4_Rep1.err
factor=$(grep "Total" 2_GATA4-fb-Adult-SRR1025222-SRR1025223_Rep2.final.rmDup.bamStats|awk 'BEGIN{sum=0}{sum+=$3}END{print 10000000/sum}')
bamCompare -b1 2_GATA4-fb-Adult-SRR1025222-SRR1025223_Rep2.final.rmDup.bam -b2 2_Input-fb-Ab-Adult-SRR1025224_Rep1.final.rmDup.bam --scaleFactors ${factor}:$input --operation subtract -bs 20 -p 10 --extendReads 160 -o GATA4_Rep2.subtract.bw --outFileFormat bigwig >log/bamCompare.Gata4_Rep2.log 2>log/bamCompare.Gata4_Rep2.err
bigwigCompare -b1 GATA4_Rep1.subtract.bw -b2 GATA4_Rep2.subtract.bw --operation mean -bs 20 -p 20 -o GATA4_RepMean.subtract.bw

mkdir result && cd result
bedtools intersect -wo -a ../Pol2_Rep1/Pol2_peaks.narrowPeak -b ../Pol2_Rep2/Pol2_peaks.narrowPeak|awk -v OFS="\t" -v ORS="" '{print $1"\t";if($2<=$12){print $2"\t";}else{print $12"\t";}if($3<=$13){print $13"\n"}else{print $3"\n"}}'|sort|uniq|sort -k1,1 -k2,2n|bedtools merge -i stdin >pol2.intersect.merged.peaks.bed
bedtools intersect -wo -a ../Gata4_Rep1/Gata4_peaks.narrowPeak -b ../Gata4_Rep2/Gata4_peaks.narrowPeak| awk -v OFS="\t" -v ORS="" '{print $1"\t";if($2<=$12){print $2"\t";}else{print $12"\t";}if($3<=$13){print $13"\n"}else{print $3"\n"}}'|sort|uniq |grep -vE "M|Un|random" >Gata4.intersect.merged.peaks.v1.bed 

#1.6 Chromatin modifiers ChIP-seq data (public)
##EED
mkdir ChIP-seq/others/EED
ls *_1.fastq |sed 's/_1.fastq//'|while read file;do
    cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -o ${file}_1.cutAda.fq -p ${file}_2.cutAda.fq -O 10 ${file}_1.fastq ${file}_2.fastq >${file}.cutAda.log 2>${file}.cutAda.err 
done;
ls *_1.cutAda.fq |sed 's/_1.cutAda.fq//'|while read file;do 
        bwa mem -t 10 /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}_1.cutAda.fq ${file}_2.cutAda.fq 2>${file}.bwa.log |samtools view -@ 5 -bSu - |samtools sort -@ 5 -o ${file}.out.sorted.bam - 
done;
ls -d */|grep -v "suz12"|while read dir;do
        cd $dir;
        ls *.out.sorted.bam |while read file;do 
        prefix=$(echo $file|sed 's/.out.sorted.bam//');
        bamtools filter -in $file -out ${prefix}.filtered.sorted.bam -script /mnt/share/liym/data/bwa/bwaMem.filter.json;
        bamtools stats -in ${prefix}.filtered.sorted.bam >${prefix}.filtered.sorted.bamStats &
        samtools rmdup ${prefix}.filtered.sorted.bam ${prefix}.filtered.rmdup.bam 2>${prefix}.rmDup.log
        done;
cd ../;
done;
samtools merge -@ 5 EED_ChIP.merged.sorted.bam EED_ChIP_rep1.SRR4241119.filtered.rmdup.bam EED_ChIP_rep2.SRR4241120.filtered.rmdup.bam
input=$(grep "Total" EED_ChIP_input.SRR4241121.filtered.rmdup.bamStats |awk '{print 10000000/$3}')
factor=$(grep "Total" EED_ChIP.merged.sorted.bamStats|awk '{print 10000000/$3}');
bamCompare -b1 EED_ChIP.merged.sorted.bam -b2 EED_ChIP_input.SRR4241121.filtered.rmdup.bam --scaleFactors ${factor}:$input --ratio subtract -bs 10 -p 10 -o EED.merged.subtract.bw --outFileFormat bigwig >merged.bamCompare.log 2>merged..bamCompare.err
factorIP=$(grep "Total" EED_ChIP.merged.sorted.bamStats|awk '{print 10000000/$3}');
bamCoverage -b EED_ChIP.merged.sorted.bam -o EED_ChIP.merged.signal.bw -of bigwig --scaleFactor $factorIP -bs 20 -p 10 >EED_ChIP.merged.bamCoverage.log 2>EED_ChIP.merged.bamCoverage.err;
factorInput=$(grep "Total" EED_ChIP_input.SRR4241121.filtered.rmdup.bamStats|awk '{print 10000000/$3}');
bamCoverage -b EED_ChIP_input.SRR4241121.filtered.rmdup.bam -o EED_ChIP_input.SRR4241121.signal.bw -of bigwig --scaleFactor $factorInput -bs 20 -p 10 >EED_ChIP_input.bamCoverage.log 2>EED_ChIP_input.bamCoverage.err 

##P300
mkdir ChIP-seq/others/p300_ENCODE
#Input
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX143/SRX143847/SRR489733/SRR489733.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX143/SRX143847/SRR489734/SRR489734.sra
mv SRR489733.fastq Heart-Input-adult-8wks_Rep1.SRR489733.fq
mv SRR489734.fastq Heart-Input-adult-8wks_Rep2.SRR489734.fq 
#P300
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX143/SRX143839/SRR489717/SRR489717.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX143/SRX143839/SRR489718/SRR489718.sra
mv SRR489717.fastq Heart-P300-adult-8wks_Rep1.SRR489717.fq
mv SRR489718.fastq Heart-P300-adult-8wks_Rep2.SRR489718.fq
fqAdapterFilter.pl -a GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTGAAA -i Heart-Input-adult-8wks_Rep1.SRR489733.fastq -o Heart-Input-adult-8wks_Rep1.SRR489733.adaF.fq -r Heart-Input-adult-8wks_Rep1.SRR489733.adaF.log 
fqSeFilter.pl -f 0.1 -l 30 -q 0.5 -c 10 -a 20 Heart-Input-adult-8wks_Rep1.SRR489733.adaF.fq >Heart-Input-adult-8wks_Rep1.SRR489733.Filter.fq 2>Heart-Input-adult-8wks_Rep1.SRR489733.fqFilter.log 
fqSeFilter.pl -f 0.1 -l 30 -q 0.5 -c 10 -a 20 Heart-Input-adult-8wks_Rep2.SRR489734.fastq >Heart-Input-adult-8wks_Rep2.SRR489734.Filter.fq 2>Heart-Input-adult-8wks_Rep2.SRR489734.fqFilter.log
fqSeFilter.pl -f 0.1 -l 30 -q 0.5 -c 10 -a 20 -r 20 Heart-P300-adult-8wks_Rep1.SRR489717.fastq >Heart-P300-adult-8wks_Rep1.SRR489717.Filter.fq 2>Heart-P300-adult-8wks_Rep1.SRR489717.fqFilter.log
fqAdapterFilter.pl -a GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTGAAA -i Heart-P300-adult-8wks_Rep2.SRR489718.fastq -o Heart-P300-adult-8wks_Rep2.SRR489718.adaF.fq -r Heart-P300-adult-8wks_Rep2.SRR489718.adaF.log
fqSeFilter.pl -f 0.1 -l 30 -q 0.5 -c 10 -a 20 Heart-P300-adult-8wks_Rep2.SRR489718.adaF.fq >Heart-P300-adult-8wks_Rep2.SRR489718.Filter.fq 2>Heart-P300-adult-8wks_Rep2.SRR489718.fqFilter.log
ls *.fq|sed 's/.fq//'|while read file;do 
    bwa aln -t 5 /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}.fq >${file}.out.sai 2>${file}.aln.log;
        bwa samse /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}.out.sai ${file}.fq 2>${file}.samse.log |samtools view -bSu -| /home/liym/.local/bin/samtools sort -@ 3 -o ${file}.out.sorted.bam - ;
        bamtools filter -in ${file}.out.sorted.bam -out ${file}.uniq.sorted.bam -script /mnt/share/liym/data/bwa/bwaAln.filter.json 2>${file}.bamtoolsFilter.log
        samtools rmdup -s ${file}.uniq.sorted.bam ${file}.final.rmDup.bam 2>${file}.rmDup.log
done;
macs2 callpeak -t Heart-P300-adult-8wks_Rep1.SRR489717.Filter.final.rmDup.bam -c Heart-Input-adult-8wks_Rep1.SRR489733.Filter.final.rmDup.bam -f BAM -g mm --outdir p300_rep1 -n p300 >Heart-P300-adult-8wks_Rep1.SRR489717.macs2.log 2>Heart-P300-adult-8wks_Rep1.SRR489717.macs2.err
macs2 callpeak -t Heart-P300-adult-8wks_Rep2.SRR489718.Filter.final.rmDup.bam -c Heart-Input-adult-8wks_Rep2.SRR489734.Filter.final.rmDup.bam -f BAM -g mm --outdir p300_rep2 -n p300 >Heart-P300-adult-8wks_Rep2.SRR489718.log 2>Heart-P300-adult-8wks_Rep2.SRR489718.err

##SUZ12
mkdir ChIP-seq/others/suz12
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTCGTATGCCGTCTTCTGCTTG -i SUZ12-ChIP-heart-WT-SRR1297213.fastq -o SUZ12-ChIP-heart-WT-SRR1297213.adaTrim.fq -r SUZ12-ChIP-heart-WT-SRR1297213.adaTrim.log
fqAdapterFilter.pl -a CGGAAGAGCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAA,GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTGAAAAAAA -i Input-ChIP-heart-WT-SRR1297212.fastq -o Input-ChIP-heart-WT-SRR1297212.adaTrim.fq -r Input-ChIP-heart-WT-SRR1297212.adaTrim.log
fqSeFilter.pl -f 0.1 -b N -l 30 -q 0.5 -c 10 -a 20 Input-ChIP-heart-WT-SRR1297212.adaTrim.fq >Input-ChIP-heart-WT-SRR1297212.fqFilter.fq 2>Input-ChIP-heart-WT-SRR1297212.fqFilter.log
bwa aln -t 10 /mnt/share/liym/data/bwa/mm10/mm10.fa SUZ12-ChIP-heart-WT-SRR1297213.adaTrim.fq >SUZ12-ChIP-heart-WT-SRR1297213.sai 2>SUZ12-ChIP-heart-WT-SRR1297213.bwaAln.log
bwa aln -t 10 /mnt/share/liym/data/bwa/mm10/mm10.fa Input-ChIP-heart-WT-SRR1297212.fqFilter.fq >Input-ChIP-heart-WT-SRR1297212.sai 2>Input-ChIP-heart-WT-SRR1297212.bwaAln.log
bwa samse /mnt/share/liym/data/bwa/mm10/mm10.fa SUZ12-ChIP-heart-WT-SRR1297213.sai SUZ12-ChIP-heart-WT-SRR1297213.adaTrim.fq 2>SUZ12-ChIP-heart-WT-SRR1297213.bwaSamse.log|samtools view -bSu - |samtools sort -@ 5 -o SUZ12-ChIP-heart-WT-SRR1297213.out.sorted.bam -
bwa samse /mnt/share/liym/data/bwa/mm10/mm10.fa Input-ChIP-heart-WT-SRR1297212.sai Input-ChIP-heart-WT-SRR1297212.fqFilter.fq 2>Input-ChIP-heart-WT-SRR1297212.bwaSamse.log|samtools view -bSu - |samtools sort -@ 5 -o Input-ChIP-heart-WT-SRR1297212.out.sorted.bam -
ls *.out.sorted.bam |while read file;do 
        prefix=$(echo $file|sed 's/.out.sorted.bam//');
        bamtools stats -in $file >${file}Stats &
        bamtools filter -in $file -out ${prefix}.filtered.bam -script /mnt/share/liym/data/bwa/bwaAln.filter.json 2>${prefix}.bamtoolsFilter.log
        samtools rmdup -s ${prefix}.filtered.bam ${prefix}.final.rmDup.bam 2>${prefix}.rmDup.log
done;
macs2 callpeak -t SUZ12-ChIP-heart-WT-SRR1297213.final.rmDup.bam -c Input-ChIP-heart-WT-SRR1297212.final.rmDup.bam -f BAM -g mm --pvalue 1e-5 --broad-cutoff 0.01 --broad --outdir suz12_broad_4 -n suz12 > SUZ12-ChIP-heart-WT-SRR1297213.macs2_broad_4.log 2> SUZ12-ChIP-heart-WT-SRR1297213.macs2_broad_4.err

##HDAC1
mkdir ChIP-seq/others/HDAC1
ls *_1.fq |sed 's/_1.fq//'|while read file;do cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -o ${file}_1.cutAda.fq -p ${file}_2.cutAda.fq -O 10 ${file}_1.fq ${file}_2.fq >${file}.cutAda.log 2>${file}.cutAda.err;done;
ls *_1.cutAda.fq |sed 's/_1.cutAda.fq//'|while read file;do 
        bwa mem -t 10 /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}_1.cutAda.fq ${file}_2.cutAda.fq 2>${file}.bwa.log |samtools view -@ 5 -bSu - |samtools sort -@ 5 -o ${file}.out.sorted.bam - 
done;
ls *.out.sorted.bam |while read file;do 
        prefix=$(echo $file|sed 's/.out.sorted.bam//');
        bamtools filter -in $file -out ${prefix}.filtered.sorted.bam -script /mnt/share/liym/data/bwa/bwaMem.filter.json;
        bamtools stats -in ${prefix}.filtered.sorted.bam >${prefix}.filtered.sorted.bamStats &
        samtools rmdup ${prefix}.filtered.sorted.bam ${prefix}.filtered.rmdup.bam 2>${prefix}.rmDup.log
done;
macs2 callpeak -t HDAC1_ChIP_EEDhete.lib1643_HVLCNCCXX_L4.filtered.rmdup.bam -c ../HDAC2/HDAC2_input_EEDhete.SRR4241137.filtered.rmdup.bam -f BAMPE -g mm --pvalue 1e-5 --broad-cutoff 1e-3 --broad --outdir HDAC1_rep1_broad_2 -n HDAC1_rep1 >HDAC1_ChIP_EEDhete.lib1643.macs2.log 2>HDAC1_ChIP_EEDhete.lib1643.macs2.err 
macs2 callpeak -t HDAC1_ChIP_EEDhomo.lib1641_HVLCNCCXX_L4.filtered.rmdup.bam -c ../HDAC2/HDAC2_input_EEDhete.SRR4241137.filtered.rmdup.bam -f BAMPE -g mm --pvalue 1e-5 --broad-cutoff 1e-3 --broad --outdir HDAC1_rep2_broad_2 -n HDAC1_rep2 >HDAC1_ChIP_EEDhomo.lib1641.macs2.log 2>HDAC1_ChIP_EEDhomo.lib1641.macs2.err

##HDAC2
mkdir ChIP-seq/others/HDAC2
ls *_1.fastq |sed 's/_1.fastq//'|while read file;do
    cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -o ${file}_1.cutAda.fq -p ${file}_2.cutAda.fq -O 10 ${file}_1.fastq ${file}_2.fastq >${file}.cutAda.log 2>${file}.cutAda.err 
done;
ls *_1.cutAda.fq |sed 's/_1.cutAda.fq//'|while read file;do 
        bwa mem -t 5 /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}_1.cutAda.fq ${file}_2.cutAda.fq 2>${file}.bwa.log |samtools view -@ 5 -bSu - |samtools sort -@ 5 -o ${file}.out.sorted.bam - 
done;
ls *.out.sorted.bam |while read file;do 
        prefix=$(echo $file|sed 's/.out.sorted.bam//');
        bamtools filter -in $file -out ${prefix}.filtered.sorted.bam -script /mnt/share/liym/data/bwa/bwaMem.filter.json;
        bamtools stats -in ${prefix}.filtered.sorted.bam >${prefix}.filtered.sorted.bamStats &
        samtools rmdup ${prefix}.filtered.sorted.bam ${prefix}.filtered.rmdup.bam 2>${prefix}.rmDup.log
done;
macs2 callpeak -t HDAC2_ChIP_EEDhete.SRR4241136.filtered.rmdup.bam -c HDAC2_input_EEDhete.SRR4241137.filtered.rmdup.bam -f BAMPE -g mm --pvalue 1e-5 --broad-cutoff 1e-3 --broad --outdir HDAC2_rep1_broad_2 -n HDAC2_rep1 >HDAC2_ChIP_EEDhete.SRR4241136.filtered.macs2.log 2>HDAC2_ChIP_EEDhete.SRR4241136.filtered.macs2.err
macs2 callpeak -t HDAC2_ChIP_EEDhomo.lib1642_HVLCNCCXX_L4.filtered.rmdup.bam -c HDAC2_input_EEDhete.SRR4241137.filtered.rmdup.bam -f BAMPE -g mm --pvalue 1e-5 --broad-cutoff 1e-3 --broad --outdir HDAC2_rep2_broad_2 -n HDAC2_rep2 >HDAC2_ChIP_EEDhomo.lib1642.macs2.log 2>HDAC2_ChIP_EEDhomo.lib1642.macs2.err 

#1.7 MNase-seq data in EEDCKO and WT
mkdir -p MNase/MNaseSeq-EEDko-EEDhete
cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -o 1_MNase-EEDhetoMF-lib1612_Rep2_1.cutAda.fq -p 1_MNase-EEDhetoMF-lib1612_Rep2_2.cutAda.fq -O 10 -u 1 -U 1 1_MNase-EEDhetoMF-lib1612_Rep2_1.fastq 1_MNase-EEDhetoMF-lib1612_Rep2_2.fastq > 1_MNase-EEDhetoMF-lib1612.cutAda.log 2> 1_MNase-EEDhetoMF-lib1612.cutAda.err
cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -o 1_MNase-EEDkoMF-lib1613_Rep2_1.cutAda.fq -p 1_MNase-EEDkoMF-lib1613_Rep2_2.cutAda.fq -O 10 -u 1 -U 1 1_MNase-EEDkoMF-lib1613_Rep2_1.fastq 1_MNase-EEDkoMF-lib1613_Rep2_2.fastq > 1_MNase-EEDkoMF-lib1613.cutAda.log 2> 1_MNase-EEDkoMF-lib1613.cutAda.err
cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -o 1_MNase-EEDhetoMF-lib1610_Rep1_1.cutAda.fq -p 1_MNase-EEDhetoMF-lib1610_Rep1_2.cutAda.fq -O 10 -u 1 -U 1 1_MNase-EEDhetoMF-lib1610_Rep1_1.fastq 1_MNase-EEDhetoMF-lib1610_Rep1_2.fastq >1_MNase-EEDhetoMF-lib1610.cutAda.log 2>1_MNase-EEDhetoMF-lib1610.cutAda.err
cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -o 1_MNase-EEDkoMF-lib1611_Rep1_1.cutAda.fq -p 1_MNase-EEDkoMF-lib1611_Rep1_2.cutAda.fq -O 10 -u 1 -U 1 1_MNase-EEDkoMF-lib1611_Rep1_1.fastq 1_MNase-EEDkoMF-lib1611_Rep1_2.fastq >1_MNase-EEDkoMF-lib1611.cutAda.log 2>1_MNase-EEDkoMF-lib1611.cutAda.err
bwa mem -t 8 /mnt/share/liym/data/bwa/mm10/mm10.fa 1_MNase-EEDhetoMF-lib1612_Rep2_1.cutAda.fq 1_MNase-EEDhetoMF-lib1612_Rep2_2.cutAda.fq 2> 1_MNase-EEDhetoMF-lib1612.bwa.log | samtools view -bSu - | samtools sort - 1_MNase-EEDhetoMF-lib1612.out.sorted
bwa mem -t 8 /mnt/share/liym/data/bwa/mm10/mm10.fa 1_MNase-EEDkoMF-lib1613_Rep2_1.cutAda.fq 1_MNase-EEDkoMF-lib1613_Rep2_2.cutAda.fq 2> 1_MNase-EEDkoMF-lib1613.log | samtools view -bSu - | samtools sort - 1_MNase-EEDkoMF-lib1613.out.sorted
bwa mem -t 8 /mnt/share/liym/data/bwa/mm10/mm10.fa 1_MNase-EEDhetoMF-lib1610_Rep1_1.cutAda.fq 1_MNase-EEDhetoMF-lib1610_Rep1_2.cutAda.fq 2> 1_MNase-EEDhetoMF-lib1610.bwa.log |samtools view -bSu - | samtools sort - 1_MNase-EEDhetoMF-lib1610.out.sorted 
bwa mem -t 8 /mnt/share/liym/data/bwa/mm10/mm10.fa 1_MNase-EEDkoMF-lib1611_Rep1_1.cutAda.fq 1_MNase-EEDkoMF-lib1611_Rep1_2.cutAda.fq 2> 1_MNase-EEDkoMF-lib1611.bwa.log |samtools view -bSu - | samtools sort - 1_MNase-EEDkoMF-lib1611.out.sorted
ls *.out.sorted.bam|while read file;do
        prefix=$(echo $file|sed 's/.out.sorted.bam//')
        bamtools filter -in $file -out ${prefix}.myFiltered.sorted.bam -script /mnt/share/liym/data/bwa/bwaMem.filter.json 2>log/bamtoolsFilter.log;
        bamtools stats -in ${prefix}.myFiltered.sorted.bam >${prefix}.myFiltered.sorted.bamStats & 
        samtools index ${prefix}.myFiltered.sorted.bam;
        samtools view -@ 5 -bu ${prefix}.myFiltered.sorted.bam chr1 chr2 chrX chr3 chr4 chr5 chr6 chr7 chr10 chr8 chr14 chr9 chr11 chr13 chr12 chr15 chr16 chr17 chrY chr18 chr19 chrM |samtools sort -@ 5 -o ${prefix}.final - 2>log/samtoolsSort.log;
        samtools rmdup ${prefix}.final.bam ${prefix}.final.rmDup.bam 2>${prefix}.rmDup.log 
        bamtools stats -in ${prefix}.final.bam >${prefix}.final.bamStats & 
done;
ls *.out.sorted.bam |sed 's/.out.sorted.bam//'|while read file;do java -jar /mnt/share/share/tool/picard-tools/picard.jar CollectInsertSizeMetrics I=${file}.out.sorted.bam O=${file}_insert_size_metrics.txt H=${file}_insert_size_histogram.pdf >log/${file}.collectInsertS.log 2>log/${file}.collectInsertS.err;done;
ls *.final.rmDup.bam|sed 's/.bam//'|while read file;do bamtools filter -in ${file}.bam -out ${file}.SizeF.bam -script bwaMem.insertSize.json ; bamtools stats -insert -in ${file}.SizeF.bam >${file}.SizeF.bamStats & done;
danpos.py dpos 1_MNase-EEDhetoMF-lib1610.final.rmDup.SizeF.bam,1_MNase-EEDhetoMF-lib1612.final.rmDup.SizeF.bam,1_MNase-EEDkoMF-lib1611.final.rmDup.SizeF.bam,1_MNase-EEDkoMF-lib1613.final.rmDup.SizeF.bam -jd 147 --extend 74 -m 1 -c 10000000 -o danpos_rst > danpos.log 2> danpos.err
ls danpos_rst/pooled/*.wig|while read file;do prefix=$(echo $file|sed 's/.final.rmDup.SizeF.Fnor.smooth.wig//');wigToBigWig -clip $file /rd1/user/liym/nucleosomeTurnover/data/mm10.filterRandom.size ${prefix}.bw 2>wigTobw.log ;done;

#1.8 Brg1 ChIP-seq data (public)
mkdir ChIP-seq/others/public_Brg1_ChIP/Brg1_GSE37151
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX190/SRX190171/SRR577648/SRR577648.sra
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX190/SRX190172/SRR577649/SRR577649.sra
mv SRR577648.fastq mouse_e11.5_heart_smarca4-flag.SRR577648.fq 
mv SRR577649.fastq mouse_e11.5_heart_control.SRR577649.fq
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGC,AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATG -i mouse_e11.5_heart_smarca4-flag.SRR577648.fq -o mouse_e11.5_heart_smarca4-flag.SRR577648.adaF.fq -r mouse_e11.5_heart_smarca4-flag.SRR577648.adaF.log 
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGC,AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATG -i mouse_e11.5_heart_control.SRR577649.fq -o mouse_e11.5_heart_control.SRR577649.adaF.fq -r mouse_e11.5_heart_control.SRR577649.adaF.log
ls *adaF.fq |while read file;do
        prefix=$(echo $file|sed 's/.fq//');
        bwa aln -t 10 /mnt/share/liym/data/bwa/mm10/mm10.fa $file >${prefix}.sai 2>${prefix}.bwaAln.log
        bwa samse /mnt/share/liym/data/bwa/mm10/mm10.fa ${prefix}.sai $file 2>${prefix}.bwaSamse.log|samtools view -bSu - |samtools sort -@ 5 -o ${prefix}.out.sorted.bam -
done;
ls *.out.sorted.bam |while read file;do 
    prefix=$(echo $file|sed 's/.out.sorted.bam//');
    bamtools stats -in $file >${file}Stats &
    bamtools filter -in $file -out ${prefix}.filtered.bam -script /mnt/share/liym/data/bwa/bwaAln.filter.json 2>${prefix}.bamtoolsFilter.log;
    samtools rmdup -s ${prefix}.filtered.bam ${prefix}.final.rmDup.bam 2>${prefix}.rmDup.log
done;
macs2 callpeak -t mouse_e11.5_heart_smarca4-flag.SRR577648.adaF.final.rmDup.bam -c mouse_e11.5_heart_control.SRR577649.adaF.final.rmDup.bam --broad --pvalue 1e-5 --broad-cutoff 1e-3 -f BAM -g mm --outdir macs2_v2 -n Brg1 >macs2_v2.log 2>macs2_v2.err
awk '$8>4' macs2_v2/*broadPeak >mouse_e11.5_heart_smarca4-flag.p1e-4.broadPeak
Rscript /mnt/share/liym/tools/phantompeakqualtools/run_spp.R -c=mouse_e11.5_heart_smarca4-flag.SRR577648.adaF.final.rmDup.bam -s=0:5:1500 -p=10 -savp=mouse_e11.5_heart_smarca4-flag.SRR577648.spp.pdf -out=mouse_e11.5_heart_smarca4-flag.SRR577648.spp.tsv >mouse_e11.5_heart_smarca4-flag.SRR577648.spp.log 2>mouse_e11.5_heart_smarca4-flag.SRR577648.spp.err  
bamCompare -b1 mouse_e11.5_heart_smarca4-flag.SRR577648.adaF.final.rmDup.bam -b2 mouse_e11.5_heart_control.SRR577649.adaF.final.rmDup.bam -o mouse_e11.5_heart_smarca4-flag.SRR577648.subtract.bw --ratio subtract -e 170 -p 10 >bamCompare.log 2>bamCompare.err
computeMatrix reference-point -R mouse_e11.5_heart_smarca4-flag.p1e-4.broadPeak -S *.bw -out PublicPeak.Brg1.ChIPsignals.matrix.gz --referencePoint center --sortRegions no -b 5000 -a 5000 -bs 20 -p 10 
plotHeatmap -m PublicPeak.Brg1.ChIPsignals.matrix.gz -out PublicPeak.Brg1.ChIPsignals.heatmap.pdf --sortUsingSamples 1 --colorList white,red --colorNumber 50 --samplesLabel inhouse.v1.EEDhete.rep1 inhouse.v1.EEDhete.rep2 inhouse.v2.EEDhete.rep1 inhouse.v2.EEDhete.rep2 inhouse.v3.EEDhete.rep1 inhouse.v3.EEDhete.rep2 publicData --regionsLabel "" --refPointLabel "0" --yAxisLabel "ChIP signals" --xAxisLabel "Peak distance(bp)"
computeMatrix reference-point -R ../../data/EEDmerged.peaks.v2.bed3 -S mouse_e11.5_heart_smarca4-flag.SRR577648.subtract.bw -out EEDpeaks.publicBrg1Signals.matrix.gz --referencePoint center --sortRegions no -b 5000 -a 5000 -bs 20 -p 10
plotHeatmap -m EEDpeaks.publicBrg1Signals.matrix.gz -out EEDpeaks.publicBrg1Signals.heatmap.pdf --sortUsingSamples 1 --colorList white,red --colorNumber 50 --samplesLabel Brg1-ChIP --regionsLabel "" --refPointLabel "0" --yAxisLabel "ChIP signals" --xAxisLabel "Peak distance(bp)"
gunzip EEDpeaks.publicBrg1Signals.matrix.gz
Rscript /mnt/share/liym/bin/rescale_linear.R -i=EEDpeaks.publicBrg1Signals.matrix -s=7 -e=506 -n=0,5 -o=EEDpeaks.publicBrg1Signals.rescale.matrix
cat <(head -n1 EEDpeaks.publicBrg1Signals.matrix) EEDpeaks.publicBrg1Signals.rescale.matrix >tmp;mv tmp EEDpeaks.publicBrg1Signals.rescale.matrix
gzip EEDpeaks.publicBrg1Signals.rescale.matrix
plotHeatmap -m EEDpeaks.publicBrg1Signals.rescale.matrix.gz -out EEDpeaks.publicBrg1Signals.rescale.matrix.pdf --colorList white,red --colorNumber 10 --sortUsingSamples 1 --samplesLabel Brg1-ChIP --regionsLabel "" --refPointLabel "0" --yAxisLabel "ChIP signals" --xAxisLabel "Peak distance(bp)"

#1.9 H2BGFP ChIP-seq data in Brg1KD and WT
mkdir ChIP-seq/others/H2BGFP-Brg1KD-week0-20180914
mv lib18167_TKD180900266_1.fq.gz H2BGFP-CHIP-Brg1KD-Adult-lib18167_Rep2_1.fq.gz
mv lib18167_TKD180900266_2.fq.gz H2BGFP-CHIP-Brg1KD-Adult-lib18167_Rep2_2.fq.gz
mv lib18169_TKD180900268_1.fq.gz H2BGFP-CHIP-Control-Adult-lib18169_Rep1_1.fq.gz
mv lib18169_TKD180900268_2.fq.gz H2BGFP-CHIP-Control-Adult-lib18169_Rep1_2.fq.gz
mv lib18170_TKD180900269_1.fq.gz H2BGFP-CHIP-Control-Adult-lib18170_Rep2_1.fq.gz
mv lib18170_TKD180900269_2.fq.gz H2BGFP-CHIP-Control-Adult-lib18170_Rep2_2.fq.gz
mv lib18171_TKD180900270-13_1.fq.gz Input-Brg1KD-Adult-lib18171_1.fq.gz
mv lib18171_TKD180900270-13_2.fq.gz Input-Brg1KD-Adult-lib18171_2.fq.gz
mv lib18172_TKD180900270-12_1.fq.gz Input-Control-Adult-lib18172_1.fq.gz
mv lib18172_TKD180900270-12_2.fq.gz Input-Control-Adult-lib18172_2.fq.gz
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGC -i H2BGFP-CHIP-Control-Adult-lib18170_Rep2_1.fq -o H2BGFP-CHIP-Control-Adult-lib18170_Rep2_adaF_1.fq --a2 GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCG --i2 H2BGFP-CHIP-Control-Adult-lib18170_Rep2_2.fq --o2 H2BGFP-CHIP-Control-Adult-lib18170_Rep2_adaF_2.fq -r H2BGFP-CHIP-Control-Adult-lib18170_Rep2.adaF.log
ls *_1.fq.gz|grep -vE "lib18168|lib18170"|while read file;do
        filePre=$(echo $file|sed 's/_1.fq.gz//');
        trim_galore -q 20 --fastqc --illumina --stringency 6 -e 0 --gzip --length 50 --trim-n --paired ${filePre}_1.fq.gz ${filePre}_2.fq.gz >log/${filePre}_trim.log 2>log/${filePre}_trim.err
        done;
ls *_1_val_1.fq.gz|while read file;do
     prefix=$(echo $file|sed 's/_1_val_1.fq.gz//'|sed 's/_adaF//')
     mkdir $prefix && cd $prefix;
        bwa mem -t 15 /mnt/share/liym/data/bwa/mm10/mm10.filterRandom.fa ../${prefix}*_1_val_1.fq.gz ../${prefix}*_2_val_2.fq.gz 2>${prefix}.bwa.log |samtools view -@ 10 -bSu - |/home/liym/.local/bin/samtools sort -@ 10 -o ${prefix}.out.sorted.bam -;
        bamtools filter -in ${prefix}.out.sorted.bam -out ${prefix}.filtered.sorted.bam -script /mnt/share/liym/data/bwa/bwaMem.filter.json 2>${prefix}.bamtoolsFilter.log;
        samtools view ${prefix}.filtered.sorted.bam |grep -E "SA|XA" |cut -f1 >multAligName.txt;
        bamtools stats -insert -in ${prefix}.filtered.sorted.bam >${prefix}.filtered.sorted.bamStats & 
        perl /mnt/share/liym/bin/bam_filterByReadname.pl -i ${prefix}.filtered.sorted.bam -r multAligName.txt |samtools view -@ 20 -bSu - |/home/liym/.local/bin/samtools sort -@ 10 -o ${prefix}.filtered.uniq.bam -
        bamtools stats -insert -in ${prefix}.filtered.uniq.bam >${prefix}.filtered.uniq.bamStats &
        java -jar /mnt/share/share/tool/picard-tools-2.2.4/picard.jar MarkDuplicates I=${prefix}.filtered.uniq.bam O=${prefix}.markDup.final.bam M=${prefix}.markDup.metrix.txt REMOVE_DUPLICATES=true >${prefix}.markDup.log 2>${prefix}.markDup.err
        java -jar /mnt/share/share/tool/picard-tools-2.2.4/picard.jar CollectGcBiasMetrics I=${prefix}.markDup.final.bam O=${prefix}.GCbias.metrix.txt CHART=${prefix}.GCbias.metrix.pdf S=${prefix}.GCbias.summary.metrix.txt R=/mnt/share/liym/data/bwa/mm10/mm10.filterRandom.fa > ${prefix}.GCbias.log 2> ${prefix}.GCbias.err
        cd ../;
done;
ls */*GCbias.metrix.txt|sed 's/.txt//'|while read file;do Rscript /mnt/share/liym/bin/picard.GCbiasMetrix.plot.R -i=${file}.txt -y1=0 -y2=7 -o=${file}.inhouse.pdf;done;
Rscript /mnt/share/liym/bin/picard.GCbiasMetrix.plot.R -i=H2BGFP-CHIP-Brg1KD-Adult-lib18166_Rep1/H2BGFP-CHIP-Brg1KD-Adult-lib18166_Rep1.GCbias.metrix.txt -y1=0 -y2=35 -o=H2BGFP-CHIP-Brg1KD-Adult-lib18166_Rep1/H2BGFP-CHIP-Brg1KD-Adult-lib18166_Rep1.GCbias.metrix.inhouse.pdf
Rscript /mnt/share/liym/bin/picard.GCbiasMetrix.plot.R -i=H2BGFP-CHIP-Brg1KD-Adult-lib18168_Rep3/H2BGFP-CHIP-Brg1KD-Adult-lib18168_Rep3.GCbias.metrix.txt -y1=0 -y2=35 -o=H2BGFP-CHIP-Brg1KD-Adult-lib18168_Rep3/H2BGFP-CHIP-Brg1KD-Adult-lib18168_Rep3.GCbias.metrix.inhouse.pdf
#Calculate nucleosome occupancy 
danpos.py dpos H2BGFP-CHIP-Brg1KD-Adult-lib18166_Rep1/H2BGFP-CHIP-Brg1KD-Adult-lib18166_Rep1.markDup.final.bam,H2BGFP-CHIP-Brg1KD-Adult-lib18167_Rep2/H2BGFP-CHIP-Brg1KD-Adult-lib18167_Rep2.markDup.final.bam,H2BGFP-CHIP-Brg1KD-Adult-lib18168_Rep3/H2BGFP-CHIP-Brg1KD-Adult-lib18168_Rep3.markDup.final.bam -b Input-Brg1KD-Adult-lib18171/Input-Brg1KD-Adult-lib18171.markDup.final.bam -c 10000000 --extend 74 -o danposRst_Brg1kd_week0 >danposRst_Brg1kd_week0.log 2>danposRst_Brg1kd_week0.err 
danpos.py dpos H2BGFP-CHIP-Control-Adult-lib18169_Rep1/H2BGFP-CHIP-Control-Adult-lib18169_Rep1.markDup.final.bam,H2BGFP-CHIP-Control-Adult-lib18170_Rep2/H2BGFP-CHIP-Control-Adult-lib18170_Rep2.markDup.final.bam -b Input-Control-Adult-lib18172/Input-Control-Adult-lib18172.markDup.final.bam -c 10000000 --extend 74 -o danposRst_control_week0 >danposRst_control_week0.log 2>danposRst_control_week0.err
ls danposRst_*/pooled/*.wig|sed 's/.wig//'|while read file;do wigToBigWig -clip ${file}.wig /mnt/share/liym/data/chr.size/mm10.filterRandom.size ${file}.bw 2>wigTobw.err; done;
#Calculate fragment size distribution 
ls H2BGFP-CHIP-*/*final.bam|while read file;do
        dir=$(dirname $file);
        java -jar /mnt/share/share/tool/picard-tools-2.2.4/picard.jar CollectInsertSizeMetrics I=$file O=${dir}/${dir}_insert_size_metrics.txt H=${dir}/${dir}_insert_size_histogram.pdf >${dir}/${dir}.collectInsertS.log 2>${dir}/${dir}.collectInsertS.err
done;
#Generate signals by using deeptools 
ls */*markDup.final.bam|grep -v "Input"|while read file;do 
        prefix=$(echo $file|sed 's/.markDup.final.bam//');
        condition=$(dirname $file|cut -f3 -d '-');
        factor=$(grep "Total" ${file}Stats|awk '{print 10000000/$3}')
        input=$(grep "Total" Input-${condition}-*/*markDup.final.bamStats|awk '{print 10000000/$3}')
        bamCompare -b1 $file -b2 Input-${condition}-*/*.markDup.final.bam --scaleFactors ${factor}:$input --ratio subtract -bs 20 -p 15 --extendReads -o ${prefix}.subtract.bw --outFileFormat bigwig
done;

mkdir ChIP-seq/others/H2BGFP-Brg1KD-week3-lib1826-Final
ls lib1826-6_HCV5TDMXX_L1_*|sed 's/lib1826-6_HCV5TDMXX_L1_//'|while read file;do mv lib1826-6_HCV5TDMXX_L1_$file H2BGFP_Input_scramble_li1826-6_$file;done;
ls lib1826-5_HCV5TDMXX_L1_*|sed 's/lib1826-5_HCV5TDMXX_L1_//'|while read file;do mv lib1826-5_HCV5TDMXX_L1_$file H2BGFP_Input_Brg1kd_li1826-5_$file;done;
ls lib1826-4_HCV5TDMXX_L1_*|sed 's/lib1826-4_HCV5TDMXX_L1_//'|while read file;do mv lib1826-4_HCV5TDMXX_L1_$file H2BGFP_ChIP_scramble_li1826-4_rep1_$file;done;
ls lib1826-3_HCV5TDMXX_L1_*|sed 's/lib1826-3_HCV5TDMXX_L1_//'|while read file;do mv lib1826-3_HCV5TDMXX_L1_$file H2BGFP_ChIP_scramble_li1826-3_rep2_$file;done;
ls lib1826-2_HCV5TDMXX_L1_*|sed 's/lib1826-2_HCV5TDMXX_L1_//'|while read file;do mv lib1826-2_HCV5TDMXX_L1_$file H2BGFP_ChIP_Brg1kd_li1826-2_rep1_$file;done;
ls lib1826-1_HCV5TDMXX_L1_*|sed 's/lib1826-1_HCV5TDMXX_L1_//'|while read file;do mv lib1826-1_HCV5TDMXX_L1_$file H2BGFP_ChIP_Brg1kd_li1826-1_rep2_$file;done;
ls *_1.fq.gz|while read file;do
        filePre=$(echo $file|sed 's/_1.fq.gz//');
        trim_galore -q 20 --fastqc --illumina --stringency 6 -e 0 --gzip --length 50 --trim-n --paired ${filePre}_1.fq.gz ${filePre}_2.fq.gz >log/${filePre}_trim.log 2>log/${filePre}_trim.err
done;

ls *_1_val_1.fq.gz|while read file;do
     prefix=$(echo $file|sed 's/_1_val_1.fq.gz//')
     mkdir $prefix && cd $prefix
        bwa mem -t 15 /mnt/share/liym/data/bwa/mm10/mm10.filterRandom.fa ../${prefix}_1_val_1.fq.gz ../${prefix}_2_val_2.fq.gz 2>${prefix}.bwa.log |samtools view -@ 10 -bSu - |samtools sort -@ 10 - ${prefix}.out.sorted;
        bamtools filter -in ${prefix}.out.sorted.bam -out ${prefix}.filtered.sorted.bam -script /mnt/share/liym/data/bwa/bwaMem.filter.json 2>${prefix}.bamtoolsFilter.log;
        samtools view ${prefix}.filtered.sorted.bam |grep -E "SA|XA" |cut -f1 >multAligName.txt;
        bamtools stats -insert -in ${prefix}.filtered.sorted.bam >${prefix}.filtered.sorted.bamStats & 
        perl /mnt/share/liym/bin/bam_filterByReadname.pl -i ${prefix}.filtered.sorted.bam -r multAligName.txt |samtools view -@ 20 -bSu - |/home/liym/miniconda3/bin/samtools sort -@ 20 -o ${prefix}.filtered.uniq.bam -
        bamtools stats -insert -in ${prefix}.filtered.uniq.bam >${prefix}.filtered.uniq.bamStats &
        java -jar /mnt/share/share/tool/picard-tools-2.2.4/picard.jar MarkDuplicates I=${prefix}.filtered.uniq.bam O=${prefix}.markDup.final.bam M=${prefix}.markDup.metrix.txt REMOVE_DUPLICATES=true >${prefix}.markDup.log 2>${prefix}.markDup.err
        java -jar /mnt/share/share/tool/picard-tools-2.2.4/picard.jar CollectGcBiasMetrics I=${prefix}.markDup.final.bam O=${prefix}.GCbias.metrix.txt CHART=${prefix}.GCbias.metrix.pdf S=${prefix}.GCbias.summary.metrix.txt R=/mnt/share/liym/data/bwa/mm10/mm10.filterRandom.fa > ${prefix}.GCbias.log 2> ${prefix}.GCbias.err
        cd ../;
done;
ls H2BGFP_ChIP*/*final.bam|while read file;do
        dir=$(dirname $file);
        java -jar /mnt/share/share/tool/picard-tools-2.2.4/picard.jar CollectInsertSizeMetrics I=$file O=${dir}/${dir}_insert_size_metrics.txt H=${dir}/${dir}_insert_size_histogram.pdf >${dir}/${dir}.collectInsertS.log 2>${dir}/${dir}.collectInsertS.err
done;
ls */*GCbias.metrix.txt|sed 's/.txt//'|while read file;do  Rscript /mnt/share/liym/bin/picard.GCbiasMetrix.plot.R -i=${file}.txt -y1=0 -y2=4 -o=${file}.inhouse.pdf;done;

#Calculate nucleosome occupancy 
danpos.py dpos H2BGFP_ChIP_Brg1kd_li1826-1_rep2/H2BGFP_ChIP_Brg1kd_li1826-1_rep2.markDup.final.bam,H2BGFP_ChIP_Brg1kd_li1826-2_rep1/H2BGFP_ChIP_Brg1kd_li1826-2_rep1.markDup.final.bam -b H2BGFP_Input_Brg1kd_li1826-5/H2BGFP_Input_Brg1kd_li1826-5.markDup.final.bam -c 10000000 --extend 74 -o danposRst_Brg1kd_week3 >danposRst_Brg1kd_week3.log 2>danposRst_Brg1kd_week3.err 
danpos.py dpos H2BGFP_ChIP_scramble_li1826-3_rep2/H2BGFP_ChIP_scramble_li1826-3_rep2.markDup.final.bam,H2BGFP_ChIP_scramble_li1826-4_rep1/H2BGFP_ChIP_scramble_li1826-4_rep1.markDup.final.bam -b H2BGFP_Input_scramble_li1826-6/H2BGFP_Input_scramble_li1826-6.markDup.final.bam -c 10000000 --extend 74 -o danposRst_scramble_week3 >danposRst_scramble_week3.log 2>danposRst_scramble_week3.err
ls danposRst_*/pooled/*.wig|sed 's/.wig//'|while read file;do wigToBigWig -clip ${file}.wig /mnt/share/liym/data/chr.size/mm10.filterRandom.size ${file}.bw 2>wigTobw.err & done;
#Merged bam for two replicates 
mkdir mergedBam && cd mergedBam
samtools merge H2BGFP_ChIP_Brg1kd.merged.bam ../H2BGFP_ChIP_Brg1kd_li1826-1_rep2/H2BGFP_ChIP_Brg1kd_li1826-1_rep2.markDup.final.bam ../H2BGFP_ChIP_Brg1kd_li1826-2_rep1/H2BGFP_ChIP_Brg1kd_li1826-2_rep1.markDup.final.bam
samtools merge H2BGFP_ChIP_scramble.merged.bam ../H2BGFP_ChIP_scramble_li1826-3_rep2/H2BGFP_ChIP_scramble_li1826-3_rep2.markDup.final.bam ../H2BGFP_ChIP_scramble_li1826-4_rep1/H2BGFP_ChIP_scramble_li1826-4_rep1.markDup.final.bam
danpos.py dpos H2BGFP_ChIP_Brg1kd.merged.bam -b ../H2BGFP_Input_Brg1kd_li1826-5/H2BGFP_Input_Brg1kd_li1826-5.markDup.final.bam -c 10000000 --extend 74 -m 1 -o danposRst_Brg1kd_week3 >danposRst_Brg1kd_week3.log 2>danposRst_Brg1kd_week3.err
danpos.py dpos H2BGFP_ChIP_scramble.merged.bam -b ../H2BGFP_Input_scramble_li1826-6/H2BGFP_Input_scramble_li1826-6.markDup.final.bam -c 10000000 --extend 74 -m 1 -o danposRst_scramble_week3 >danposRst_scramble_week3.log 2>danposRst_scramble_week3.err
ls */pooled/*.wig |grep -v "Input"|sed 's/.wig//'|while read file;do wigToBigWig -clip ${file}.wig /mnt/share/liym/data/chr.size/mm10.filterRandom.size ${file}.bw;done;

#1.10 MNase-seq data in Brg1KD and WT
mkdir MNase/MNaseSeq-Brg1KD-week0-20180914
ls lib18173*|while read file;do prefix=$(echo $file|cut -f1 -d '_');suffix=$(echo $file|cut -f4- -d '_');mv $file MNase-Brg1KD-week0-${prefix}_Rep1_${suffix};done;
ls lib18174*|while read file;do prefix=$(echo $file|cut -f1 -d '_');suffix=$(echo $file|cut -f4- -d '_');mv $file MNase-Brg1KD-week0-${prefix}_Rep2_${suffix};done;
ls lib18175*|while read file;do prefix=$(echo $file|cut -f1 -d '_');suffix=$(echo $file|cut -f4- -d '_');mv $file MNase-Control-week0-${prefix}_Rep1_${suffix};done;
ls lib18176*|while read file;do prefix=$(echo $file|cut -f1 -d '_');suffix=$(echo $file|cut -f4- -d '_');mv $file MNase-Control-week0-${prefix}_Rep2_${suffix};done;
ls *_1.fq.gz|while read file;do
        prefix=$(echo $file|sed 's/_1.fq.gz//');
        mkdir $prefix;
        cd $prefix;
        mv ../${prefix}_* .;
        mkdir log;
        trim_galore -q 20 --fastqc --illumina --stringency 6 -e 0 --gzip --length 50 --trim-n --paired ${prefix}_1.fq.gz ${prefix}_2.fq.gz >log/${prefix}_trim.log 2>log/${prefix}_trim.err 
        cd ../;
done;

ls -d */|sed 's;/;;'|while read prefix;do
        cd $prefix;
        bwa mem -t 15 /mnt/share/liym/data/bwa/mm10/mm10.filterRandom.fa ${prefix}_1_val_1.fq.gz ${prefix}_2_val_2.fq.gz 2>log/${prefix}.bwa.log |samtools view -@ 10 -bSu - |samtools sort -@ 10 -o ${prefix}.out.sorted.bam -;
        java -jar /mnt/share/share/tool/picard-tools-2.2.4/picard.jar CollectInsertSizeMetrics I=${prefix}.out.sorted.bam O=${prefix}_insert_size_metrics.txt H=${prefix}_insert_size_histogram.pdf >log/${prefix}.collectInsertS.log 2>log/${prefix}.collectInsertS.err
        bamtools filter -in ${prefix}.out.sorted.bam -out ${prefix}.filtered.sorted.bam -script /mnt/share/liym/data/bwa/bwaMem.filter.json 2>log/${prefix}.bamtoolsFilter.log;
        samtools view ${prefix}.filtered.sorted.bam |grep -E "SA|XA" |cut -f1 >multAligName.txt;
        bamtools stats -insert -in ${prefix}.filtered.sorted.bam >${prefix}.filtered.sorted.bamStats & 
        perl /mnt/share/liym/bin/bam_filterByReadname.pl -i ${prefix}.filtered.sorted.bam -r multAligName.txt |samtools view -@ 20 -bSu - |samtools sort -@ 20 -o ${prefix}.filtered.uniq.bam -
        bamtools stats -insert -in ${prefix}.filtered.uniq.bam >${prefix}.filtered.uniq.bamStats &
        java -jar /mnt/share/share/tool/picard-tools-2.2.4/picard.jar MarkDuplicates I=${prefix}.filtered.uniq.bam O=${prefix}.markDup.final.bam M=${prefix}.markDup.metrix.txt REMOVE_DUPLICATES=true >log/${prefix}.markDup.log 2>log/${prefix}.markDup.err
        bamtools stats -insert -in ${prefix}.markDup.final.bam >${prefix}.markDup.final.bamStats
        perl /mnt/share/liym/nucleosome/scripts/peReadsFragmentMid.pl -i ${prefix}.markDup.final.bam -f 147,147 >n147_mid.txt 
        perl /mnt/share/liym/nucleosome/scripts/region_dinucleotideFreq.pl -b n147_mid.txt -n AA,AT,TA,TT -u -150 -d 150 -f /mnt/share/liym/data/genome/mm10/mm10.fa -c >AT_freq.tsv 2>AT_freq.log &
        perl /mnt/share/liym/nucleosome/scripts/region_dinucleotideFreq.pl -b n147_mid.txt -n GG,GC,CG,CC -u -150 -d 150 -f /mnt/share/liym/data/genome/mm10/mm10.fa -c >GC_freq.tsv 2>GC_freq.log
        if [ -s AT_freq.tsv && -s GC_freq.tsv ];then
            /mnt/share/liym/nucleosome/scripts/dinucleotide_plot.R AT_freq.tsv GC_freq.tsv ${prefix}.dinucleotide.pdf 2>log/dinucleotide_R.log
        fi
        cd ../;
done;
ls */*.markDup.final.bam|sed 's/.bam//'|while read file;do 
        dir=$(dirname $file);
        #bamtools filter -in ${file}.bam -out ${file}.SizeF.bam -script /mnt/share/liym/data/bwa/bwaMem.insertSize.json & 
        #bamtools stats -insert -in ${file}.SizeF.bam >${file}.SizeF.bamStats  
        cd $dir;
        danpos.py dpos ../${file}.bam -c 10000000 --extend 74 -m 1 -o danposRst >danpos2.log 2>danpos2.err
        cd ../;
done;
ls */*/*/*.wig |sed 's/.wig//'|while read file;do wigToBigWig -clip ${file}.wig /mnt/share/liym/data/chr.size/mm10.filterRandom.size ${file}.bw;done;

#2. HTR calcualtion (WT samples; 6 time points; corresponding to Figure1, Figure2 and part of Figure 3)
mkdir H2BGFP_All/figures && cd H2BGFP_All/figures
mkdir fig1
toolRunner.sh wigmath.Average -f --step 10 -p 2 -o week0-new.wig /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/1/danposRst/day1_input2/pooled/1_H2bGFP-day1_Rep1.final.bgsub.Fnor.smooth.wig /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/1/danposRst/day1_input2/pooled/1_H2bGFP-day1_Rep2.final.bgsub.Fnor.smooth.wig > log/week0-new.wigmath.log
wigToBigWig -clip week0-new.wig /rd1/user/liym/nucleosomeTurnover/data/mm10.filterRandom.size week0-new.bw
toolRunner.sh wigmath.Average -f --step 10 -p 2 -o week1.wig /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/1/danposRst/result/pooled/1_H2bGFP-week1_Rep1.final.bgsub.Fnor.smooth.wig /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/1/danposRst/result/pooled/1_H2bGFP-week1_Rep2.final.bgsub.Fnor.smooth.wig > week1.wigmath.log
toolRunner.sh wigmath.Average -f --step 10 -p 2 -o week2.wig /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/1/danposRst/result/pooled/1_H2bGFP-week2_Rep1.final.bgsub.Fnor.smooth.wig /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/1/danposRst/result/pooled/1_H2bGFP-week2_Rep2.final.bgsub.Fnor.smooth.wig > week2.wigmath.log
toolRunner.sh wigmath.Average -f --step 10 -p 2 -o week4.wig /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/1/danposRst/result/pooled/1_H2bGFP-week4_Rep1.final.bgsub.Fnor.smooth.wig /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/2/danposRst/result/pooled/2_H2BGFP-CM-week4-2015_Rep1.final.bgsub.Fnor.smooth.wig /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/2/danposRst/result/pooled/2_H2BGFP-CM-week4-2015_Rep2.final.bgsub.Fnor.smooth.wig > week4.wigmath.log
ls *.wig |sed 's/.wig//'|while read file;do wigToBigWig -clip ${file}.wig /mnt/share/liym/data/chr.size/mm10.filterRandom.size ${file}.bw 2>log/${file}.wigTobw.log;done;
ls *.bw |sed 's/.bw//'|while read file;do bwtool summary /rd1/user/liym/nucleosomeTurnover/data/mm10.1kbIntervals.filterGap.bed ${file}.bw ${file}.1kb.tsv;done;
toolRunner.sh wigmath.Average -f --step 10 -p 2 -o week0-new.wig /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/1/danposRst/day1_input2/pooled/1_H2bGFP-day1_Rep1.final.bgsub.Fnor.smooth.wig /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/1/danposRst/day1_input2/pooled/1_H2bGFP-day1_Rep2.final.bgsub.Fnor.smooth.wig > log/week0-new.wigmath.log
wigToBigWig -clip week0-new.wig /rd1/user/liym/nucleosomeTurnover/data/mm10.filterRandom.size week0-new.bw
ls *.bw |sed 's/.bw//'|while read file;do bwtool summary /rd1/user/liym/nucleosomeTurnover/data/mm10.1kbIntervals.filterGap.bed ${file}.bw ${file}.1kb.tsv;done;
mkdir geneReion && cd geneReion
/mnt/share/liym/tools/deepTools-2.0.0/bin/computeMatrix scale-regions -R ../../fig2/EEDhete.Meanfpkm.bed6 -S ../H3.bw ../week0-new.bw ../week1.bw ../week2.bw ../week4.bw ../week6.bw ../week8.bw -out geneRegions.matrix.gz --outFileNameMatrix geneRegions.matrix.tsv --regionBodyLength 10000 --startLabel TSS --endLabel TTS -b 5000 -a 5000 -bs 40 --skipZeros -p 15 >log/computeMatrix.log 2>log/computeMatrix.err
###Start R
data<-read.delim(file="geneRegions.matrix",header=F,comment.char="@")
test<-as.matrix(data[,c(7:3506)])
tmp<-test[,c(1:500)]           
tmp[tmp>2.0702414]=2.0702414
test[,c(1:500)]<-tmp
tmp<-test[,c(501:1000)]
tmp[tmp>1.704329]=1.704329  
test[,c(501:1000)]<-tmp
tmp<-test[,c(1001:1500)]
tmp[tmp>1.718702]=1.718702
test[,c(1001:1500)]<-tmp
tmp<-test[,c(1501:2000)]
tmp[tmp>1.364334]=1.364334
test[,c(1501:2000)]<-tmp
tmp<-test[,c(2001:2500)]
tmp[tmp>1.005376]=1.005376
test[,c(2001:2500)]<-tmp
tmp<-test[,c(2501:3000)]
tmp[tmp>1.005376]=1.005376
test[,c(2501:3000)]<-tmp
tmp<-test[,c(3001:3500)]
tmp[tmp>1.005376]=1.005376
test[,c(3001:3500)]<-tmp
test[test<0.01]=0.01
write.table(test,file="final.matrix.v2",quote=F,sep="\t",na = "nan",row.names=F,col.names=F)
###End R
paste (grep -v "@" geneRegions.matrix|cut -f1-6) final.matrix.v2 >tmp;
cat header tmp >final.matrix.v2
gzip -k final.matrix.v2
plotHeatmap -m final.matrix.v2.gz -out final.matrix.rowMean.sort.heatmap.v3.pdf --sortRegions descend --sortUsing mean --colorList white,red --startLabel TSS --endLabel TTS -x "" -y "" --samplesLabel H3 0w 1w 2w 4w 6w 8w --heatmapHeight 16 -max 1 -min 0 -yMin 0 -yMax 1 
plotHeatmap -m final.matrix.v2.pol2Sort.gz -out final.matrix.rowMean.sort.heatmap.v9.pdf --sortRegions descend --sortUsing mean --colorMap YlGnBu --startLabel TSS --endLabel TTS -x "" -y "" --samplesLabel H3 0w 1w 2w 4w 6w 8w --heatmapHeight 16 --zMax 0.8 --zMin 0 --missingDataColor yellow
mkdir H2BGFP_All/figures/fig2
###HTR and gene expression
awk -v OFS="\t" '{if($1!="tracking_id"){print $1":"$7,$10}}' /rd1/user/liym/nucleosomeTurnover/RNA-seq/EEDheto_1/cufflinks/genes.fpkm_tracking >EEDheto_1.fpkm.tsv
awk -v OFS="\t" '{if($1!="tracking_id"){print $1":"$7,$10}}' /rd1/user/liym/nucleosomeTurnover/RNA-seq/EEDheto_2/cufflinks/genes.fpkm_tracking >EEDheto_2.fpkm.tsv
cut -f2,3,4,12 /rd1/user/liym/nucleosomeTurnover/data/mm10.refGene.filterRandom.gpe|sort|uniq|awk -v OFS="\t" '{print $4":"$1":"$3,$2}' >mm10.refGene.start.tsv 
join.pl -i1 EEDheto_1.fpkm.tsv -i2 EEDheto_2.fpkm.tsv |awk -v OFS="\t" '{split($1,a,":");split(a[3],b,"-");print a[2],b[1],b[2],a[1]":"a[2]":"b[1],($2+$4)/2}'|join.pl -f1 4 -i2 mm10.refGene.start.tsv -f2 1|awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$7}'|sort|uniq >EEDhete.Meanfpkm.bed6 
awk -v OFS="\t" '{if($6=="+"){print $1,$2-500,$2+500,$4,$5,$6}else{print $1,$3-500,$3+500,$4,$5,$6}}' EEDhete.Meanfpkm.bed6 >EEDhete.Meanfpkm.TSSfl500.bed6 
ls ../fig1/week*.bw |grep -v "week0.bw"|while read file;do prefix=$(basename $file|sed 's/.bw//'|sed 's/-new//');cut -f1-4 EEDhete.Meanfpkm.TSSfl500.bed6 |bigWigAverageOverBed $file stdin ${prefix}.TSSfl500.H2BGFP.tsv;done;
paste <(cut -f1,5 week0.TSSfl500.H2BGFP.tsv) <(cut -f5 week1.TSSfl500.H2BGFP.tsv) <(cut -f5 week2.TSSfl500.H2BGFP.tsv) <(cut -f5 week4.TSSfl500.H2BGFP.tsv) <(cut -f5 week6.TSSfl500.H2BGFP.tsv) <(cut -f5 week8.TSSfl500.H2BGFP.tsv)|join.pl -f1 1 -i2 EEDhete.Meanfpkm.TSSfl500.bed6 -f2 4|cut -f1-7,12|sort -k8 -gr >TSS.fl500.H2BGFP.fpkmSort.tsv
Rscript /rd1/user/liym/nucleosomeTurnover/scripts/NTR.R -i=TSS.fl500.H2BGFP.fpkmSort.tsv -p=1e-4 -s=2 -e=7 -o=TSS.fl500.H2BGFP.fpkmSort.NTR.tsv 
###HTR and histone modifications
ls /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/*.bed |grep -vE "H3K27ac|H3K27me3"|while read file;do prefix=$(basename $file|sed 's/.ChIP.8w.mouseHeart.narrowPeak.bed//');computeMatrix reference-point -R $file -S ../fig1/week0-new.bw ../fig1/week1.bw ../fig1/week2.bw ../fig1/week4.bw ../fig1/week6.bw ../fig1/week8.bw -out ${prefix}.H2BGFP.gz --referencePoint center -b 5000 -a 5000 -bs 20 --sortRegions no -p 20 >computeMatrix.log 2>computeMatrix.err;done;
computeMatrix reference-point -R /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/result/H3K27ac-WT.intersect.broadPeaks.bed -S ../fig1/week0-new.bw ../fig1/week1.bw ../fig1/week2.bw ../fig1/week4.bw ../fig1/week6.bw ../fig1/week8.bw -out H3K27ac.H2BGFP.gz --referencePoint center -b 5000 -a 5000 -bs 20 --sortRegions no -p 20 >computeMatrix.log 2>computeMatrix.err
computeMatrix reference-point -R /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/result/H3K27me3-WT.intersect.broadPeaks.bed -S ../fig1/week0-new.bw ../fig1/week1.bw ../fig1/week2.bw ../fig1/week4.bw ../fig1/week6.bw ../fig1/week8.bw -out H3K27me3.H2BGFP.gz --referencePoint center -b 5000 -a 5000 -bs 20 --sortRegions no -p 20 >computeMatrix.log 2>computeMatrix.err
ls *.H2BGFP.gz|sed 's/.gz//'|while read file;do 
        mkdir $file;
        mv ${file}.gz $file;
        cd $file;
        gunzip ${file}.gz;
        ../NTR.for.lines.sh $file 
        Rscript /mnt/share/liym/bin/columnMean.R -i=*summary.NTR.tsv -s=7 -e=506 -o=colMean.NTR.tsv
done;
paste H3*/colMean.NTR.tsv >HM.lines.NTR.tsv
ls /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/*.bed|grep -vE "H3K27ac|H3K27me3"|while read file;do 
        prefix=$(basename $file|sed 's/.ChIP.8w.mouseHeart.narrowPeak.bed//');
        #/rd1/user/liym/nucleosomeTurnover/scripts/FileForNTR.sh $file ../fig1/week0-new.bw ../fig1/week1.bw ../fig1/week2.bw ../fig1/week4.bw ../fig1/week6.bw ../fig1/week8.bw ${prefix}.H2BGFP.tsv;
        Rscript /rd1/user/liym/nucleosomeTurnover/scripts/NTR.R -i=${prefix}.H2BGFP.tsv -p=0.0001 -s=4 -e=9 -o=${prefix}.NTR.tsv 
        done;
/rd1/user/liym/nucleosomeTurnover/scripts/FileForNTR.sh /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/result/H3K27me3-WT.intersect.broadPeaks.bed ../fig1/week0-new.bw ../fig1/week1.bw ../fig1/week2.bw ../fig1/week4.bw ../fig1/week6.bw ../fig1/week8.bw H3K27me3.peaks.H2BGFP.tsv 
Rscript /rd1/user/liym/nucleosomeTurnover/scripts/NTR.R -i=H3K27me3.peaks.H2BGFP.tsv -p=0.0001 -s=4 -e=9 -o=H3K27me3.peaks.NTR.tsv 
##Intergenic background
sh /rd1/user/liym/NTR/scripts/FileForNTR.sh /rd1/user/liym/NTR/data/mm10.refGene.intergenic.1kb.bed3 ../../fig3/week0.bw ../../fig3/week1.bw ../../fig3/week2.bw ../../fig3/week4.bw ../../fig3/week6.bw ../../fig3/week8.bw mm10.refGene.intergenic.1k.H2BGFP.tsv
grep -v "NA" mm10.refGene.intergenic.1k.H2BGFP.tsv >tmp;mv tmp mm10.refGene.intergenic.1k.H2BGFP.tsv
Rscript /rd1/user/liym/NTR/scripts/NTR.R -i=mm10.refGene.intergenic.1k.H2BGFP.tsv -p=0.00001 -s=4 -e=9 -o=mm10.refGene.intergenic.1k.NTR.tsv
###HTR and TF
bedtools intersect -wao -a ../../fig2/enhancers.pol2.status.NTR.tsv -b ../TFs.sum.bed3 >enhancers.TFs.intersect.bed
awk '$7>0' enhancers.TFs.intersect.bed |awk -v OFS="" -v ORS="" '{print $1"\t";if($2<$7){print $7"\t";}else{print $2"\t"}if($3<$8){print $3"\n"}else{print $8"\n"}}'|awk -v OFS="\t" '{print $1,$2,$3,"wt"}' |sort -k1,1 -k2,2n |bedtools merge >enhancers.TFs.status.bed3+
/rd1/user/liym/NTR/scripts/FileForNTR.sh enhancers.TFs.status.bed3+ ../week0.bw ../week1.bw ../week2.bw ../week4.bw ../week6.bw ../week8.bw enhancers.TFs.wt.H2BGFP.tsv 
Rscript /rd1/user/liym/NTR/scripts/NTR.R -i=enhancers.TFs.wt.H2BGFP.tsv -s=4 -e=9 -p=0.0001 -o=enhancers.TFs.wt.NTR.tsv 
awk '$7<0' enhancers.TFs.intersect.bed |cut -f1-3,5|sort|uniq|awk -v OFS="\t" '{print $1,$2,$3,"wo",$4}' >enhancers.TFs.status.NTR.tsv
awk -v OFS="\t" '{print $1,$2,$3,"wt",$10}' enhancers.TFs.wt.NTR.tsv >>enhancers.TFs.status.NTR.tsv

bedtools multiinter -i <(sort -k1,1 -k2,2n ../../fig2/enhancers.pol2.status.NTR.tsv) <(sort -k1,1 -k2,2n ../Gata4.narrowPeak.bed) <(sort -k1,1 -k2,2n ../Nkx2-5.narrowPeak.bed) <(sort -k1,1 -k2,2n ../Tbx20.narrowPeak.bed) <(sort -k1,1 -k2,2n ../Tbx3.narrowPeak.bed)|awk -v OFS="\t" '{if($5~"1" && $5~","){print $1,$2,$3,$4}}' >enhancers.TFnum.bed3+
/rd1/user/liym/NTR/scripts/FileForNTR.sh enhancers.TFnum.bed3+ ../week0.bw ../week1.bw ../week2.bw ../week4.bw ../week6.bw ../week8.bw enhancers.TFnum.H2BGFP.tsv
Rscript /rd1/user/liym/NTR/scripts/NTR.R -i=enhancers.TFnum.H2BGFP.tsv -s=4 -e=9 -p=0.0001 -o=enhancers.TFnum.NTR.tsv
paste enhancers.TFnum.bed3+ <(cut -f10 enhancers.TFnum.NTR.tsv) >enhancers.TFnum.NTR.final.tsv
####TF signals
ls ../*.signal.bw |sed 's/.signal.bw//'|sed 's;../;;'|while read file;do bwtool summary enhancers.TFnum.bed3+ ../${file}.signal.bw enhancers.TFnum.${file}.signal.tsv;done;
paste *TFnum*signal.tsv |cut -f1-3,8,17,16,25|bedtools intersect -wo -a stdin -b enhancers.TFnum.NTR.final.tsv|awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$11,$12}' >enhancers.TFnum.allTFsignal.NTR.tsv
####Intergenic regions 
/rd1/user/liym/NTR/scripts/FileForNTR.sh /rd1/user/liym/NTR/data/mm10.refGene.intergenic.1kb.bed3 ../week0.bw ../week1.bw ../week2.bw ../week4.bw ../week6.bw ../week8.bw mm10.refGene.intergenic.1kb.H2BGFP.tsv
grep -v "NA" mm10.refGene.intergenic.1kb.H2BGFP.tsv >tmp;mv tmp mm10.refGene.intergenic.1kb.H2BGFP.tsv
Rscript /rd1/user/liym/NTR/scripts/NTR.R -i=mm10.refGene.intergenic.1kb.H2BGFP.tsv -s=4 -e=9 -p=0.0001 -o=mm10.refGene.intergenic.1kb.NTR.tsv 
cut -f1-3,10,11 mm10.refGene.intergenic.1kb.NTR.tsv >tmp;mv tmp mm10.refGene.intergenic.1kb.NTR.tsv
####NTR lines for TF peaks
ls *Peak.bed |while read file;do prefix=$(echo $file|sed 's/.narrowPeak.bed//');computeMatrix reference-point -R $file -S week0.bw week1.bw week2.bw week4.bw week6.bw week8.bw -out ${prefix}.H2BGFP.gz --referencePoint center -b 5000 -a 5000 -bs 20 --sortRegions no -p 20 >${prefix}.computeMatrix.log 2>${prefix}.computeMatrix.err;done;
ls *.H2BGFP.gz|sed 's/.H2BGFP.gz//'|while read file;do 
        mkdir $file;
        mv ${file}.H2BGFP.gz $file;
        cd $file;
        gunzip ${file}.H2BGFP.gz;
        /rd1/user/liym/nucleosomeTurnover/scripts/NTR.for.lines.sh ${file}.H2BGFP $file 
        Rscript /mnt/share/liym/bin/columnMean.R -i=${file}.summary.NTR.tsv -s=7 -e=506 -o=${file}.colMean.NTR.tsv
        cd ../;
done;
paste */*colMean.NTR.tsv >Gata4.Nkx2-5.Tbx20.Tbx3.agg.txt 
##HTR for enhancers bound or not bound by chromatin modifiers
mkdir H2BGFP_All/figures/fig4
ls *.peaks.bed|sed 's/.peaks.bed//'|while read file;do cut -f1-3 ${file}.peaks.bed|bedtools intersect -wao -a enhancers.NTR.bed3+ -b stdin >${file}.enhancers.intersect.bed;done;
cd suz12
cut -f1-3 /rd1/user/liym/nucleosomeTurnover/ChIP-seq/others/suz12/suz12_broad_4/suz12_peaks.broadPeak |bedtools intersect -wao -a enhancers.NTR.bed3+ -b stdin >suz12.enhancers.intersect.V2.bed
bedtools intersect -v -a /rd1/user/liym/nucleosomeTurnover/ChIP-seq/others/suz12/suz12_broad_4/suz12_peaks.broadPeak -b ../enhancers.NTR.bed3+ |cut -f1-3 >suz12.unique.V2.bed3 
awk '$6>0' suz12.enhancers.intersect.V2.bed |awk -v OFS="" -v ORS="" '{print $1"\t";if($2<=$6){print $6"\t";}else{print $2"\t";}if($3<=$7){print $3"\n"}else{print $7"\n"}}' >suz12.enhancers.shareRegion.V2.bed3 
awk '$6>0' suz12.enhancers.intersect.V2.bed |awk -v OFS="" -v ORS="" '{print $1"\t";if($2<=$6){print $2"\t";}else{print $6"\t";}if($3<=$7){print $7"\n"}else{print $3"\n"}}' >suz12.enhancers.mergeRegion.V2.bed3 
#P300
/rd1/user/liym/nucleosomeTurnover/scripts/macs2Peaks.union.sh narrow /rd1/user/liym/nucleosomeTurnover/ChIP-seq/TF/ENCODE/p300_rep1/p300_peaks.narrowPeak /rd1/user/liym/nucleosomeTurnover/ChIP-seq/TF/ENCODE/p300_rep2/p300_peaks.narrowPeak p300.merged.peaks.V0.bed
prefix=p300
version=V0
bedtools intersect -wao -a ../enhancers.NTR.bed3+ -b ${prefix}.merged.peaks.${version}.bed >${prefix}.enhancers.intersect.${version}.bed; 
awk '$6>0' ${prefix}.enhancers.intersect.${version}.bed |awk -v OFS="" -v ORS="" '{print $1"\t";if($2<=$6){print $6"\t";}else{print $2"\t";}if($3<=$7){print $3"\n"}else{print $7"\n"}}' |awk '$1!~"random|GL|Un"' >${prefix}.enhancers.shareRegion.${version}.bed3;
awk '$6>0' ${prefix}.enhancers.intersect.${version}.bed |awk -v OFS="" -v ORS="" '{print $1"\t";if($2<=$6){print $2"\t";}else{print $6"\t";}if($3<=$7){print $7"\n"}else{print $3"\n"}}' |awk '$1!~"random|GL|Un"' >${prefix}.enhancers.mergeRegion.${version}.bed3;
bedtools intersect -v -a ${prefix}.merged.peaks.${version}.bed -b ../enhancers.NTR.bed3+ |cut -f1-3| awk '$1!~"random|GL|Un"'>${prefix}.unique.${version}.bed3;
ls *.bed3|while read file;do prefix=$(echo $file|sed 's/.bed3//');/rd1/user/liym/nucleosomeTurnover/scripts/FileForNTR.sh $file ../week0.bw ../week1.bw ../week2.bw ../week4.bw ../week6.bw ../week8.bw ${prefix}.H2BGFP.tsv;Rscript /rd1/user/liym/nucleosomeTurnover/scripts/NTR.R -i=${prefix}.H2BGFP.tsv -p=0.0001 -s=4 -e=9 -o=${prefix}.NTR.tsv 2>NTR.err;done;
#EED
/rd1/user/liym/nucleosomeTurnover/scripts/macs2Peaks.union.sh broad /rd1/user/liym/nucleosomeTurnover/ChIP-seq/others/EED/EED_rep1_broad_2/EED_rep1_peaks.broadPeak /rd1/user/liym/nucleosomeTurnover/ChIP-seq/others/EED/EED_rep2_broad_2/EED_rep2_peaks.broadPeak EED.merged.peaks.V1.bed 
prefix=EED
version=V1
#HDAC1
cut -f1-3 /rd1/user/liym/nucleosomeTurnover/ChIP-seq/others/HDAC1/EEDhete_broad_2/HDAC1_peaks.broadPeak >HDAC1.merged.peaks.V1.bed
prefix=HDAC1
version=V1
#HDAC2
cut -f1-3 /rd1/user/liym/nucleosomeTurnover/ChIP-seq/others/HDAC2/EEDhete_broad_2/HDAC2_peaks.broadPeak >HDAC2.merged.peaks.V1.bed
prefix=HDAC2
version=V1
#EED peaks and H3K27me3
bedtools intersect -wao -a /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/result/H3K27me3-WT.intersect.broadPeaks.bed -b eed/EED.merged.peaks.V1.bed |awk '$1!~"random|GL|Un"' >H3K27me3.EED.intersect.bed; 
bedtools intersect -v -a /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/result/H3K27me3-WT.intersect.broadPeaks.bed -b eed/EED.merged.peaks.V1.bed |awk '$1!~"random|GL|Un"' >H3K27me3.unique.bed3;
bedtools intersect -v -a eed/EED.merged.peaks.V1.bed -b /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/result/H3K27me3-WT.intersect.broadPeaks.bed |awk '$1!~"random|GL|Un"' >EED.unique.bed3;
awk '$5>0' H3K27me3.EED.intersect.bed |awk -v OFS="" -v ORS="" '{print $1"\t";if($2<=$5){print $5"\t";}else{print $2"\t";}if($3<=$6){print $3"\n"}else{print $6"\n"}}' >H3K27me3.EED.shareRegion.bed3; 
ls *.bed3|while read file;do prefix=$(echo $file|sed 's/.bed3//');/rd1/user/liym/nucleosomeTurnover/scripts/FileForNTR.sh $file week0.bw week1.bw week2.bw week4.bw week6.bw week8.bw ${prefix}.H2BGFP.tsv;Rscript /rd1/user/liym/nucleosomeTurnover/scripts/NTR.R -i=${prefix}.H2BGFP.tsv -p=0.0001 -s=4 -e=9 -o=${prefix}.NTR.tsv 2>NTR.err;done;

#3 HTR change in EED CKO
mkdir H2BGFP_All/figures/fig5
##Up-regulated H3K27ac peaks
ls ../../3/danposRst/result/pooled/*.wig |grep -v "Input"|while read file;do prefix=$(basename $file|sed 's/.final.bgsub.Fnor.smooth.wig//');wigToBigWig -clip $file /mnt/share/liym/data/chr.size/mm10.chrom.sizes ${prefix}.bw 2>wigTobw.log & done;
ls ../../4/danposRst/result/pooled/*.wig|grep -v "Input"|while read file;do prefix=$(basename $file|sed 's/.final.bgsub.Fnor.smooth.wig//');wigToBigWig -clip $file /mnt/share/liym/data/chr.size/mm10.chrom.sizes ${prefix}.bw 2>wigTobw.log; done;

Rscript /mnt/share/liym/bin/vennFor2BedRegion.R -b1=H3K27ac-WT.intersect.broadPeaks.bed -b2=H3K27ac-EEDko.intersect.broadPeaks.bed -n1=WT -n2=KO -c=3
bedtools intersect -wao -a H3K27ac-EEDko.intersect.broadPeaks.bed -b H3K27ac-WT.intersect.broadPeaks.bed |awk '$5<0'|cut -f1-3 >EEDko.specific.peaks.bed3 
bedtools intersect -wao -b H3K27ac-EEDko.intersect.broadPeaks.bed -a H3K27ac-WT.intersect.broadPeaks.bed |awk '$5<0'|cut -f1-3 >WT.specific.peaks.bed3
bedtools intersect -wao -b H3K27ac-EEDko.intersect.broadPeaks.bed -a H3K27ac-WT.intersect.broadPeaks.bed |awk '$5>0'|awk -v OFS="\t" -v ORS="" '{print $1"\t";if($2<=$5){print $5"\t"}else{print $2"\t"}if($3<=$6){print $3"\n"}else{print $6"\n"}}' >WT.EEDko.common.peaks.bed3 
ls *.peaks.bed3|sed 's/.peaks.bed3//'|while read file;do
        bwtool agg 5000:5000 ${file}.peaks.bed3 3_H2BGFP-CM-week0-EEDheto_Rep1.bw,3_H2BGFP-CM-week0-EEDheto_Rep2.bw,4_H2BGFP-CM-week4-EEDheto_Rep1.bw,4_H2BGFP-CM-week4-EEDheto_Rep2.bw,3_H2BGFP-CM-week0-EEDko_Rep1.bw,4_H2BGFP-CM-week4-EEDko_Rep1.bw,4_H2BGFP-CM-week4-EEDko_Rep2.bw ${file}.H2BGFP.profile.txt
        bwtool agg 5000:5000 ${file}.peaks.bed3 MNase.WT.norm.bw,MNase.EEDKO.norm.bw ${file}.NC.profile.txt
        Rscript /mnt/share/liym/bin/lines.R -i=${file}.H2BGFP.profile.txt -x="Peak relative position" -y="Normalized H2BGFP signal" -c="#FF0000FF,#FFDB00FF,#49FF00FF,#00FF92FF,#0092FFFF,#4900FFFF,#FF00DBFF" -l="WT-week0-rep1,WT-week0-rep2,WT-week4-rep1,WT-week4-rep2,EEDko-week0,EEDko-week4-rep1,EEDko-week4-rep2" -o=${file}.H2BGFP.profile.pdf
        Rscript /mnt/share/liym/bin/lines.R -i=${file}.NC.profile.txt -x="Peak relative position" -y="Normalized nucleosome occupancy" -c="red,blue" -l="WT,EEDko" -o=${file}.NC.profile.pdf
        Rscript run.R $file
done;

ls /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/*H3K27ac*.bw |while read file;do prefix=$(basename $file|sed 's/.subtract.bw//');bwtool summary WT.EEDko.common.peaks.bed3 $file WT.EEDko.common.${prefix}.signal.tsv -skip-median & done;
 paste <(cut -f1-3,8 WT.EEDko.common.1_H3K27ac-EEDko-AdultCM_Rep1.signal.tsv) <(cut -f8 WT.EEDko.common.1_H3K27ac-EEDko-AdultCM_Rep2.signal.tsv) <(cut -f8 WT.EEDko.common.1_H3K27ac-WT-AdultCM_Rep1.signal.tsv) <(cut -f8 WT.EEDko.common.1_H3K27ac-WT-AdultCM_Rep2.signal.tsv)|awk -v OFS="\t" '{print $1,$2,$3,($4+$5)/2,($6+$7)/2}' >WT.EEDko.common.EEDko-WT.signals.tsv             
 awk '$4/$5>1.5' WT.EEDko.common.EEDko-WT.signals.tsv |cut -f1-3 >EEDko.increase.bed3
 awk '$4/$5<0.5' WT.EEDko.common.EEDko-WT.signals.tsv |cut -f1-3 >EEDko.decrease.bed3
 awk '$4/$5>=0.5 && $4/$5<=1.5' WT.EEDko.common.EEDko-WT.signals.tsv |cut -f1-3 >WT.EEDko.unchanged.bed3

ls /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/*H3K27ac*.bw |while read file;do prefix=$(basename $f
ile|sed 's/.subtract.bw//');bwtool summary /rd1/user/liym/nucleosomeTurnover/data/mm10.1kbIntervals.filterGap.bed4 $file 1kbIntervals.${prefix}.signal.tsv -skip-median & done;
paste <(cut -f1-3,8 1kbIntervals.1_H3K27ac-EEDko-AdultCM_Rep1.signal.tsv) <(cut -f8 1kbIntervals.1_H3K27ac-EEDko-AdultCM_Rep2.signal.tsv) <(cut -f8 1kbIntervals.1_H3K27ac-WT-AdultCM_Rep1.signal.tsv) <(cut -f8 1kbIntervals.1_H3K27ac-WT-AdultCM_Rep2.signal.tsv) |awk -v OFS="\t" '{print $1,$2,$3,($4+$5)/2,($6+$7)/2}' >1kbIntervals.H3K27ac.EEDko.WT.signals.tsv 

ls /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/danposRst/H3K27ac/pooled/*.bw |while read file;do prefix=$(basename $file|sed 's/.bw//');bwtool summary WT.EEDko.common.peaks.bed3 $file ${prefix}.signals.commonPeaks.tsv & done;
paste <(cut -f1-3,8 1_H3K27ac-EEDko-AdultCM_Rep1.signals.commonPeaks.tsv) <(cut -f8 1_H3K27ac-EEDko-AdultCM_Rep2.signals.commonPeaks.tsv) <(cut -f8 1_H3K27ac-WT-AdultCM_Rep1.signals.commonPeaks.tsv) <(cut -f8 1_H3K27ac-WT-AdultCM_Rep2.signals.commonPeaks.tsv)|awk -v OFS="\t" '{print $1,$2,$3,($4+$5)/2,($6+$7)/2}' >commonPeaks.signals.danpos.tsv 
ls /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/danposRst/H3K27ac/pooled/*.bw |while read file;do pre
fix=$(basename $file|sed 's/.bw//');bwtool summary /rd1/user/liym/nucleosomeTurnover/data/mm10.1kbIntervals.filterGap.bed4 $file ${prefix}.signals.1kb.tsv & done;
#scatter plot
cat EEDko.specific.peaks.bed3 WT.EEDko.common.peaks.bed3 WT.specific.peaks.bed3 >H3K27ac/H3K27ac.allPeaks.bed3

ls /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/*H3K27ac*.bw |while read file;do prefix=$(basename $file|sed 's/.subtract.bw//');bwtool summary H3K27ac.allPeaks.bed3 $file peaks.${prefix}.signal.tsv -skip-median; done;
ls ../*H2BGFP*bw |while read file;do prefix=$(basename $file|sed 's/.bw//');bwtool summary H3K27ac.allPeaks.bed3 $file peaks.${prefix}.H2BGFP.tsv -skip-median;done;
paste *.signal.tsv|cut -f1-3,8,16,24,32|awk -v OFS="\t" '{print $1,$2,$3,($4+$5)/2-($6+$7)/2}' >peaks.ko.wt.signals.diff.tsv 
paste *H2BGFP.tsv|cut -f1-3,8,16,24,32,40,48,56|awk -v OFS="\t" '{print $1,$2,$3,($4+$5)/2-($7+$8)/2,$6-($9+$10)/2}' |awk -v OFS="\t" '{print $1,$2,$3,$5-$4}' >peaks.ko.wt.NTR.diff.tsv 
bedtools intersect -wo -a peaks.ko.wt.signals.diff.tsv -b peaks.ko.wt.NTR.diff.tsv|cut -f1-4,8 >peaks.ko.wt.signal.NTR.diff.tsv
paste *.signal.tsv|cut -f1-3,8,16,24,32|awk -v OFS="\t" '{if($6+$7!=0){print $1,$2,$3,($4+$5)/($6+$7)}else{print $1,$2,$3,"inf"}}'|bedtools intersect -wo -a stdin -b peaks.ko.wt.NTR.diff.tsv|cut -f1-4,8 >peaks.ko.wt.signalFC.NTR.diff.tsv
paste *.signal.tsv|cut -f1-3,8,16,24,32|bedtools intersect -wo -a stdin -b peaks.ko.wt.NTR.diff.tsv|cut -f1-7,11 >peaks.ko.wt.repSignal.NTR.diff.tsv
ls /rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/*.bw |while read file;do prefix=$(basename $file|sed 's/.bw//');cut -f1-3 peaks.ko.wt.signalFC.NTR.diff.tsv|bwtool summary -skip-median stdin $file peaks.${prefix}.NC.tsv;done;
awk '$4>1.5' peaks.ko.wt.signalFC.NTR.diff.tsv|cut -f1-3 >peaks.ko.wt.signalFC.gt1.5.bed3
bwtool agg 5000:5000 peaks.ko.wt.signalFC.gt1.5.bed3 /rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDhetoMF-lib1610.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDhetoMF-lib1612.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDkoMF-lib1611.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDkoMF-lib1613.bw peaks.ko.wt.signalFC.gt1.5.NCagg.txt

##Down-regulated H3K27me3 peaks
bedtools intersect -wao -a H3K27me3-WT.peak.bed -b H3K27me3-EEDko.peak.bed |awk '$5<0' |cut -f1-3 >EEDko.loss.peak.bed3
bedtools intersect -wao -a H3K27me3-WT.peak.bed -b H3K27me3-EEDko.peak.bed |awk '$5>0' |awk -v OFS="" -v ORS="" '{print $1"\t";if($2<=$5){print $5"\t"}else{print $2"\t"}if($3<=$6){print $3"\n"}else{print $6"\n"}}' >EEDko.common.peak.bed3
 ls /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/1_H3K27me3*.bw|while read file;do prefix=$(basename $file|sed 's/.adaF.subtract.bw//');bwtool summary -skip-median EEDko.common.peak.bed3 $file EEDko.common.${prefix}.signal.tsv;done;
paste EEDko.common.*.tsv|cut -f1-3,8,16,24,32,40|awk -v OFS="\t" '{print $1,$2,$3,($4+$5)/2,($6+$7+$8)/3}'|awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$4/$5}' >EEDko.common.signals.ko.wt.FC.tsv 
cat EEDko.loss.peak.bed3 <(awk '$6<0.5' EEDko.common.signals.ko.wt.FC.tsv|cut -f1-3) >EEDko.loss.final.bed3
awk '$6>=0.5 && $6<=2' EEDko.common.signals.ko.wt.FC.tsv |cut -f1-3 >EEDko.unchange.final.bed3
#NTR
ls ../*.bw |while read file;do prefix=$(basename $file|sed 's/.bw//');bwtool summary -skip-median EEDko.unchange.final.bed3 $file EEDko.unchange.${prefix}.tsv;bwtool summary -skip-median EEDko.loss.final.bed3 $file EEDko.loss.final.${prefix}.tsv;done;
paste EEDko.unchange.*H2BGFP*.tsv |cut -f1-3,8,16,24,32,40,48,56|awk -v OFS="\t" '{print $1,$2,$3,($4+$5)/2-($7+$8)/2,$6-($9+$10)/2}' >EEDko.unchange.wt.ko.NTR.bed3+
paste EEDko.loss.*H2BGFP*.tsv |cut -f1-3,8,16,24,32,40,48,56|awk -v OFS="\t" '{print $1,$2,$3,($4+$5)/2-($7+$8)/2,$6-($9+$10)/2}' >EEDko.loss.wt.ko.NTR.bed3+
bwtool agg 5000:5000 EEDko.unchange.final.bed3 ../3_H2BGFP-CM-week0-EEDheto_Rep1.bw,../3_H2BGFP-CM-week0-EEDheto_Rep2.bw,../3_H2BGFP-CM-week0-EEDko_Rep1.bw,../4_H2BGFP-CM-week4-EEDheto_Rep1.bw,../4_H2BGFP-CM-week4-EEDheto_Rep2.bw,../4_H2BGFP-CM-week4-EEDko_Rep1.bw,../4_H2BGFP-CM-week4-EEDko_Rep2.bw EEDko.unchange.agg.txt 
bwtool agg 5000:5000 EEDko.loss.final.bed3 ../3_H2BGFP-CM-week0-EEDheto_Rep1.bw,../3_H2BGFP-CM-week0-EEDheto_Rep2.bw,../3_H2BGFP-CM-week0-EEDko_Rep1.bw,../4_H2BGFP-CM-week4-EEDheto_Rep1.bw,../4_H2BGFP-CM-week4-EEDheto_Rep2.bw,../4_H2BGFP-CM-week4-EEDko_Rep1.bw,../4_H2BGFP-CM-week4-EEDko_Rep2.bw EEDko.loss.agg.txt
awk -v OFS="\t" '{print $1,($2+$3)/2-($5+$6)/2,$4-($7+$8)/2}' EEDko.loss.agg.txt >EEDko.loss.wt.ko.NTR.agg.txt 
awk -v OFS="\t" '{print $1,($2+$3)/2-($5+$6)/2,$4-($7+$8)/2}' EEDko.unchange.agg.txt >EEDko.unchange.wt.ko.NTR.agg.txt 
 ls /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/1_H3K27me3*.bw|while read file;do prefix=$(basename $file|sed 's/.adaF.subtract.bw//');cut -f1-3 EEDko.unchange.final.bed3|bwtool summary -skip-median stdin $file EEDko.unchange.${prefix}.signal.tsv;done;
  ls /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/1_H3K27me3*.bw|while read file;do prefix=$(basename $file|sed 's/.adaF.subtract.bw//');cut -f1-3 EEDko.loss.final.bed3|bwtool summary -skip-median stdin $file EEDko.loss.${prefix}.signal.tsv;done;
paste *loss*signal.tsv|cut -f1-3,8,16,24,32,40|awk -v OFS="\t" '{print $1,$2,$3,($4+$5)/2,($6+$7+$8)/3,($4+$5)/2-($6+$7+$8)/3}' >EEDko.loss.ko.wt.diff.tsv 
paste *unchange*signal.tsv|cut -f1-3,8,16,24,32,40|awk -v OFS="\t" '{print $1,$2,$3,($4+$5)/2,($6+$7+$8)/3,($4+$5)/2-($6+$7+$8)/3}' >EEDko.unchange.ko.wt.diff.tsv 
awk -v OFS="\t" '{print $1,$2,$3,$5-$4}' EEDko.loss.wt.ko.NTR.bed3+ |bedtools intersect -wo -a stdin -b EEDko.loss.ko.wt.diff.tsv |cut -f1-4,10 >EEDko.loss.ko.wt.NTR.signalDiff.tsv 
awk -v OFS="\t" '{print $1,$2,$3,$5-$4}' EEDko.unchange.wt.ko.NTR.bed3+ |bedtools intersect -wo -a stdin -b EEDko.unchange.ko.wt.diff.tsv |cut -f1-4,10 >EEDko.enchange.ko.wt.NTR.signalDiff.tsv 
cat <(awk '$6>2' EEDko.common.signals.ko.wt.FC.tsv|cut -f1-3) <(bedtools intersect -wo -v -a H3K27me3-EEDko.p     eak.bed -b H3K27me3-WT.peak.bed) >EEDko.gain.final.bed3 
ls /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/1_H3K27me3*.bw|while read file;do prefix=$(basename $file|sed 's/.adaF.subtract.bw//');bwtool summary -skip-median EEDko.gain.final.bed3 $file EEDko.gain.${prefix}.signal.tsv;done;
ls ../*.bw |while read file;do prefix=$(basename $file|sed 's/.bw//');bwtool summary -skip-median EEDko.gain.final.bed3 $file EEDko.gain.${prefix}.tsv;done;
paste EEDko.gain.*H2BGFP*.tsv |cut -f1-3,8,16,24,32,40,48,56|awk -v OFS="\t" '{print $1,$2,$3,($4+$5)/2-($7+$8)/2,$6-($9+$10)/2}' |awk -v OFS="\t" '{print $1,$2,$3,$5-$4}'>EEDko.gain.wt.ko.NTRdiff.bed3+
paste *gain*signal.tsv|cut -f1-3,8,16,24,32,40|awk -v OFS="\t" '{if($6+$7+$8!=0){print $1,$2,$3,(($4+$5)/2)/(($6+$7+$8)/3)}else{print $1,$2,$3,"inf"}}' >EEDko.gain.ko.wt.FC.tsv
paste *loss*signal.tsv|cut -f1-3,8,16,24,32,40|awk -v OFS="\t" '{if($6+$7+$8!=0){print $1,$2,$3,(($4+$5)/2)/(($6+$7+$8)/3)}else{print $1,$2,$3,"inf"}}' >EEDko.loss.ko.wt.FC.tsv
paste *unchange*signal.tsv|cut -f1-3,8,16,24,32,40|awk -v OFS="\t" '{if($6+$7+$8!=0){print $1,$2,$3,(($4+$5)/2)/(($6+$7+$8)/3)}else{print $1,$2,$3,"inf"}}' >EEDko.unchange.ko.wt.FC.tsv

#All peaks
bedtools intersect -wo -a <(cat EEDko.gain.ko.wt.FC.tsv EEDko.loss.ko.wt.FC.tsv EEDko.unchange.ko.wt.FC.tsv) -b <(cat EEDko.gain.wt.ko.NTRdiff.bed3+ <(awk -v OFS="\t" '{print $1,$2,$3,$5-$4}' EEDko.loss.wt.ko.NTR.bed3+) <(awk -v OFS="\t" '{print $1,$2,$3,$5-$4}' EEDko.unchange.wt.ko.NTR.bed3+))|cut -f1-4,8 >peaks.ko.wt.signalFC.NTR.diff.tsv                                    
cat <(paste *unchange*signal.tsv|cut -f1-3,8,16,24,32,40) <(paste *gain*signal.tsv|cut -f1-3,8,16,24,32,40) <(paste *loss*signal.tsv|cut -f1-3,8,16,24,32,40) >all.peaks.ko.wt.repSignal.tsv 
bedtools intersect -wo -f 1.0 -r -a all.peaks.ko.wt.repSignal.tsv -b <(cat EEDko.gain.wt.ko.NTRdiff.bed3+ <(awk -v OFS="\t" '{print $1,$2,$3,$5-$4}' EEDko.loss.wt.ko.NTR.bed3+) <(awk -v OFS="\t" '{print $1,$2,$3,$5-$4}' EEDko.unchange.wt.ko.NTR.bed3+))|cut -f1-8,12 >peaks.ko.wt.repSignal.NTR.diff.tsv
###In R
subset2<-read.delim(file="peaks.ko.wt.repSignal.NTR.diff.tsv",header=F)
subset2[subset2$V4<0,4]=0
subset2[subset2$V5<0,5]=0
subset2[subset2$V6<0,6]=0
subset2[subset2$V7<0,7]=0
subset2[subset2$V8<0,8]=0
subset2<-as.data.frame(cbind(subset2[,c(1:3)],rowMeans(subset2[,c(4,5)])/rowMeans(subset2[,c(6,7,8)]),subset2$V9))
write.table(subset2,file="peaks.ko.wt.SignalsFCminusToZero.NTR.diff.tsv",quote=F,row.names=F,col.names=F,sep="\t")
###End R
bedtools intersect -v -a <(awk '$4<0.5' peaks.ko.wt.SignalsFCminusToZero.NTR.diff.tsv) -b <(awk '$4>5 && $5>5 && $6>5 && $7>5' peaks.wt.ko.NC.tsv) |cut -f1-3 >peaks.ko.wt.signalsFC.lt0.5.bed
bwtool agg 5000:5000 peaks.ko.wt.signalsFC.lt0.5.bed /rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDhetoMF-lib1610.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDhetoMF-lib1612.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDkoMF-lib1611.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDkoMF-lib1613.bw peaks.ko.wt.signalsFC.lt0.5.NCagg.txt

#4 Interaction between EED and BRG1; HTR change in Brg1 KD
##Co-localization of EED and BRG1 
cd revision_201809
ln -s ../ChIP-seq/others/data/EEDmerged.peaks.v2.bed3 EEDmerged.peaks.v2.bed3
ln -s ../ChIP-seq/others/public_Brg1_ChIP/Brg1_GSE37151/mouse_e11.5_heart_smarca4-flag.p1e-4.broadPeak mouse_e11.5_heart_smarca4-flag.p1e-4.broadPeak bedtools intersect -a EEDmerged.peaks.v2.bed3 -b mouse_e11.5_heart_smarca4-flag.p1e-4.broadPeak -wo |awk -v OFS="" -v ORS="" '{print $1;if($2<$5){print "\t"$5}else{print "\t"$2}if($3<$6){print "\t"$3}else{print "\t"$6}print "\n"}' >EED.Brg1.commonReion.bed3
sh /rd1/user/liym/NTR/scripts/FileForNTR.sh EEDmerged.peaks.v2.bed3 ../H2BGFP_All/figures/fig3/week0.bw ../H2BGFP_All/figures/fig3/week1.bw ../H2BGFP_All/figures/fig3/week2.bw ../H2BGFP_All/figures/fig3/week4.bw ../H2BGFP_All/figures/fig3/week6.bw ../H2BGFP_All/figures/fig3/week8.bw EEDmerged.peaks.v2.H2BGFP.tsv
Rscript /rd1/user/liym/NTR/scripts/NTR.R -i=EEDmerged.peaks.v2.H2BGFP.tsv -p=0.0001 -s=4 -e=9 -o=EEDmerged.peaks.v2.NTR.tsv
sh /rd1/user/liym/NTR/scripts/FileForNTR.sh EED.Brg1.commonReion.bed3 ../H2BGFP_All/figures/fig3/week0.bw ../H2BGFP_All/figures/fig3/week1.bw ../H2BGFP_All/figures/fig3/week2.bw ../H2BGFP_All/figures/fig3/week4.bw ../H2BGFP_All/figures/fig3/week6.bw ../H2BGFP_All/figures/fig3/week8.bw EED.Brg1.commonReion.H2BGFP.tsv
Rscript /rd1/user/liym/NTR/scripts/NTR.R -i=EED.Brg1.commonReion.H2BGFP.tsv -p=0.0001 -s=4 -e=9 -o=EED.Brg1.commonReion.NTR.tsv
bedtools intersect -a EEDmerged.peaks.v2.NTR.tsv -b mouse_e11.5_heart_smarca4-flag.p1e-4.broadPeak -v|cut -f1-3,10 >EED.unique.NTR.tsv
computeMatrix reference-point -R EEDmerged.peaks.v2.bed3 -S /rd1/user/liym/NTR/ChIP-seq/others/EED/EED.merged.subtract.bw /rd1/user/liym/NTR/ChIP-seq/others/public_Brg1_ChIP/Brg1_GSE37151/mouse_e11.5_heart_smarca4-flag.SRR577648.subtract.bw -out EEDpeaks.EED.Brg1Signals.matrix.gz --referencePoint center --sortRegions no -b 5000 -a 5000 -bs 20 -p 10
plotHeatmap -m EEDpeaks.EED.Brg1Signals.matrix.gz -out EEDpeaks.EED.Brg1Signals.heatmap.pdf --sortUsingSamples 1 --colorList white,red --colorNumber 50 --samplesLabel EED-ChIP BRG1-ChIP --regionsLabel "" --refPointLabel "0" --yAxisLabel "ChIP signals" --xAxisLabel "Peak distance(bp)"
gunzip EEDpeaks.EED.Brg1Signals.matrix.gz
Rscript /mnt/share/liym/bin/rescale_linear.R -i=EEDpeaks.EED.Brg1Signals.matrix -s=7 -e=1006 -n=0,5 -o=EEDpeaks.EED.Brg1Signals.rescale.matrix
cat <(head -n1 EEDpeaks.EED.Brg1Signals.matrix) EEDpeaks.EED.Brg1Signals.rescale.matrix >tmp;mv tmp EEDpeaks.EED.Brg1Signals.rescale.matrix
gzip EEDpeaks.EED.Brg1Signals.rescale.matrix
plotHeatmap -m EEDpeaks.EED.Brg1Signals.rescale.matrix.gz -out EEDpeaks.EED.Brg1Signals.rescale.matrix.pdf --colorList white,red --colorNumber 10 --sortUsingSamples 1 --samplesLabel EED-ChIP Brg1-ChIP --regionsLabel "" --refPointLabel "0" --yAxisLabel "ChIP signals" --xAxisLabel "Peak distance(bp)"
grep -v "^@" EEDpeaks.EED.Brg1Signals.rescale.matrix|cut -f7-506|Rscript /mnt/share/liym/bin/columnMean.R -o=EEDpeak.EEDsignals.agg.txt
grep -v "^@" EEDpeaks.EED.Brg1Signals.rescale.matrix|cut -f507-1006|Rscript /mnt/share/liym/bin/columnMean.R -o=EEDpeak.Brg1signals.agg.txt
paste EEDpeak.EEDsignals.agg.txt EEDpeak.Brg1signals.agg.txt |awk -v OFS="\t" 'BEGIN{i=-5000}{print i,$1,$2;i+=20}' >EEDpeaks.EED.Brg1Signals.rescale.agg.txt
Rscript /mnt/share/liym/bin/lines.R -i=EEDpeaks.EED.Brg1Signals.rescale.agg.txt -y1=1.5 -y2=3 -x="Distance to EED peak center (bp)" -y="ChIP signals" -c="red,blue" -l="EED,BRG1" -o=EEDpeaks.EED.Brg1Signals.rescale.agg.pdf

##HTR and nucleosome occupancy change in Brg1 KD
cd ChIP-seq/others/H2BGFP-Brg1KD-week0-20180914/result
ln -s /rd1/user/liym/NTR/ChIP-seq/others/H2BGFP-Brg1KD-week3-lib1826-Final/mergedBam/danposRst_Brg1kd_week3/pooled/H2BGFP_ChIP_Brg1kd.merged.bgsub.Fnor.smooth.bw NTR/H2BGFP_ChIP_Brg1kd.week3.bw
ln -s /rd1/user/liym/NTR/ChIP-seq/others/H2BGFP-Brg1KD-week3-lib1826-Final/danposRst_Brg1kd_week3/pooled/H2BGFP_ChIP_Brg1kd_li1826-2_rep1_H2BGFP_ChIP_Brg1kd_li1826-2_rep1.markDup.final.bgsub.Fnor.smooth.bw NTR/H2BGFP_ChIP_Brg1kd_Rep1.week3.bw
ln -s /rd1/user/liym/NTR/ChIP-seq/others/H2BGFP-Brg1KD-week3-lib1826-Final/danposRst_Brg1kd_week3/pooled/H2BGFP_ChIP_Brg1kd_li1826-1_rep2_H2BGFP_ChIP_Brg1kd_li1826-1_rep2.markDup.final.bgsub.Fnor.smooth.bw NTR/H2BGFP_ChIP_Brg1kd_Rep2.week3.bw
ln -s /rd1/user/liym/NTR/ChIP-seq/others/H2BGFP-Brg1KD-week3-lib1826-Final/mergedBam/danposRst_scramble_week3/pooled/H2BGFP_ChIP_scramble.merged.bgsub.Fnor.smooth.bw NTR/H2BGFP_ChIP_Control.week3.bw
ln -s /rd1/user/liym/NTR/ChIP-seq/others/H2BGFP-Brg1KD-week3-lib1826-Final/danposRst_scramble_week3/pooled/H2BGFP_ChIP_scramble_li1826-4_rep1_H2BGFP_ChIP_scramble_li1826-4_rep1.markDup.final.bgsub.Fnor.smooth.bw NTR/H2BGFP_ChIP_Control_Rep1.week3.bw
ln -s /rd1/user/liym/NTR/ChIP-seq/others/H2BGFP-Brg1KD-week3-lib1826-Final/danposRst_scramble_week3/pooled/H2BGFP_ChIP_scramble_li1826-3_rep2_H2BGFP_ChIP_scramble_li1826-3_rep2.markDup.final.bgsub.Fnor.smooth.bw NTR/H2BGFP_ChIP_Control_Rep2.week3.bw
ls /rd1/user/liym/NTR/ChIP-seq/others/data/*.bed3|grep -v "BlackList"|while read file;do
        prefix=$(basename $file|sed 's/.bed3//');
        #bwtool summary $file danposRst_Brg1kd_week0/pooled/H2BGFP-CHIP-Brg1KD-Adult.merged.bgsub.Fnor.smooth.bw NTR/${prefix}.Brg1KD.H2BGFP.tsv -skip-median &
        #bwtool summary $file danposRst_control_week0/pooled/H2BGFP-CHIP-Control-Adult.merged.bgsub.Fnor.smooth.bw NTR/${prefix}.Control.H2BGFP.tsv -skip-median
        bwtool summary $file NTR/H2BGFP_ChIP_Brg1kd.week3.bw NTR/${prefix}.Brg1KD.week3.H2BGFP.tsv -skip-median &
        bwtool summary $file NTR/H2BGFP_ChIP_Control.week3.bw NTR/${prefix}.Control.week3.H2BGFP.tsv -skip-median 
done;
cd NTR;
ls *.week3.H2BGFP.tsv|sed 's/.week3.H2BGFP.tsv//'|while read file;do
        paste ${file}.*H2BGFP.tsv|cut -f1-3,8,16|awk -v OFS="\t" '{print $0,$4-$5}' >${file}.NTR.tsv;
done;
cd NTR;
ls ls /rd1/user/liym/NTR/ChIP-seq/others/data/*.bed3|grep -v "BlackList"|while read file;do
        prefix=$(basename $file|sed 's/.bed3//');
        sh /mnt/share/liym/bin/bwMeanForMultiFiles.sh -r $file -b ../H2BGFP-CHIP-Brg1KD-Adult-lib18167_Rep2.bw,../H2BGFP-CHIP-Control-Adult-lib18169_Rep1.bw,../H2BGFP-CHIP-Control-Adult-lib18170_Rep2.bw,H2BGFP_ChIP_Brg1kd_Rep1.week3.bw,H2BGFP_ChIP_Brg1kd_Rep2.week3.bw,H2BGFP_ChIP_Control_Rep1.week3.bw,H2BGFP_ChIP_Control_Rep2.week3.bw -o ${prefix}.H2BGFP.sum.tsv;
        awk -v OFS="\t" '{print $0,$4-($7+$8)/2,($5+$6)/2-($9+$10)/2}' ${prefix}.H2BGFP.sum.tsv >${prefix}.kd.wt.NTR.tsv;
done;


