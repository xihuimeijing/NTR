#!/bin/sh

#This is the whole pipeline for nucleosome turnover project processing.(All data are in venus)

#1. Data processing
##1.1 ChIP-seq data
##1.1.1 H2BGFP/H3 ChIP-seq data for NTR calculation
cd /rd1/user/liym/nucleosomeTurnover/H2BGFP_All
###Pre-mapping process
cd 1/
ls *.bz2|sed 's/.bz2//'|while read file;do 
    bzip2 -d ${file}.bz2;
    fastqc -q -o fastqc/ $file
    done;
perl /mnt/share/liym/bin/fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGC -i 1_H2bGFP-week8_Rep1.fq -o 1_H2bGFP-week8_Rep1.adaFilter.fq -r 1_H2bGFP-week8_Rep1.adaFilter.log
ls *.fq|sed 's/.fq//'|while read file;do fqSeFilter.pl -f 0.1 -b N -l 50 -q 0.5 -c 20 -a 20 ${file}.fq >${file}.filter.fq 2>log/${file}.fqFilter.log ;done;
cd ../2
perl /mnt/share/liym/bin/fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCC -i 2_H2BGFP-CM-week0-2015_Rep1.fq -o 2_H2BGFP-CM-week0-2015_Rep1.filter.fq -r 2_H2BGFP-CM-week0-2015_Rep1.filter.log 
perl /mnt/share/liym/bin/fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCC -i 2_H2BGFP-CM-week0-2015_Rep2.fq -o 2_H2BGFP-CM-week0-2015_Rep2.filter.fq -r 2_H2BGFP-CM-week0-2015_Rep2.filter.log 
perl /mnt/share/liym/bin/fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCC -i 2_H2BGFP-CM-week0-2015_Rep3.fq -o 2_H2BGFP-CM-week0-2015_Rep3.filter.fq -r 2_H2BGFP-CM-week0-2015_Rep3.filter.log
perl /mnt/share/liym/bin/fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCC -i 2_H2BGFP-CM-week4-2015_Rep1.fq -o 2_H2BGFP-CM-week4-2015_Rep1.filter.fq -r 2_H2BGFP-CM-week4-2015_Rep1.filter.log
perl /mnt/share/liym/bin/fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCC -i 2_H2BGFP-CM-week4-2015_Rep2.fq -o 2_H2BGFP-CM-week4-2015_Rep2.filter.fq -r 2_H2BGFP-CM-week4-2015_Rep2.filter.log
perl /mnt/share/liym/bin/fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCC -i 2_Input-CM-2015_Rep1.fq -o 2_Input-CM-2015_Rep1.filter.fq -r 2_Input-CM-2015_Rep1.log
###Mapping and post mapping filter
cd 
ls -d */|grep -v "figures"||while read dir;do
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
###Calculate signals by DANPOS2
cd 1/
mkdir danposRst && cd danposRst
danpos.py dpos 1_H2bGFP-day0_Rep1.final.bam,1_H2bGFP-day0_Rep2.final.bam,1_H2bGFP-day1_Rep1.final.bam,1_H2bGFP-day1_Rep2.final.bam,1_H2bGFP-week12_Rep1.final.bam,1_H2bGFP-week1_Rep1.final.bam,1_H2bGFP-week1_Rep2.final.bam,1_H2bGFP-week2_Rep1.final.bam,1_H2bGFP-week2_Rep2.final.bam,1_H2bGFP-week4_Rep1.final.bam,1_H2bGFP-week6_Rep1.final.bam,1_H2bGFP-week6_Rep2.final.bam,1_H2bGFP-week8_Rep1.final.bam,1_H2bGFP-week8_Rep2.final.bam,1_H3_Rep1.final.bam -b 1_Input_Rep1.final.bam -c 10000000 --extend 74 -o ./ >danpos.log 2>danpos.err
danpos.py dpos 1_H2bGFP-day1_Rep1.final.bam,1_H2bGFP-day1_Rep2.final.bam -b ../../2/danposRst/2_Input-CM-2015_Rep1.final.bam -c 10000000 --extend 74 -o day1_input2 >day1_input2.log 2>day1_input2.err
danpos.py dpos 1_H2bGFP-day0_Rep1.final.bam,1_H2bGFP-day0_Rep2.final.bam,1_H2bGFP-week1_Rep1.final.bam,1_H2bGFP-week1_Rep2.final.bam,1_H2bGFP-week2_Rep1.final.bam,1_H2bGFP-week2_Rep2.final.bam,1_H2bGFP-week4_Rep1.final.bam,1_H2bGFP-week6_Rep2.final.bam,1_H2bGFP-week8_Rep1.final.bam,1_H3_Rep1.final.bam -b ../../2/danposRst/2_Input-CM-2015_Rep1.final.bam -c 10000000 --extend 74 -o result_input2 >danpos_input2.log 2>danpos_input2.err
cd ../2
danpos.py dpos 2_H2BGFP-CM-week0-2015_Rep1.final.bam,2_H2BGFP-CM-week0-2015_Rep2.final.bam,2_H2BGFP-CM-week0-2015_Rep3.final.bam,2_H2BGFP-CM-week4-2015_Rep1.final.bam,2_H2BGFP-CM-week4-2015_Rep2.final.bam -b 2_Input-CM-2015_Rep1.final.bam -c 10000000 --extend 74 -o ./ >danpos.log 2>danpos.err
cd ../3/
danpos.py dpos 3_H2BGFP-CM-week0-EEDheto_Rep1.final.bam,3_H2BGFP-CM-week0-EEDheto_Rep2.final.bam,3_H2BGFP-CM-week0-EEDko_Rep1.final.bam -b 3_Input-CM-week0-EED_Rep1.final.bam -c 10000000 --extend 74 -o ./ >danpos.log 2>danpos.err
cd ../4/
danpos.py dpos 4_H2BGFP-CM-week4-EEDheto_Rep1.final.bam,4_H2BGFP-CM-week4-EEDheto_Rep2.final.bam,4_H2BGFP-CM-week4-EEDko_Rep1.final.bam,4_H2BGFP-CM-week4-EEDko_Rep2.final.bam -b 4_Input-CM-week4-EED_Rep1.final.bam -c 10000000 --extend 74 -o ./ >danpos.log 2>danpos.err
cd ../5/
danpos.py dpos 5_H2BGFP-CM-week2.5-banding_Rep1.final.bam,5_H2BGFP-CM-week2.5-banding_Rep2.final.bam,5_H2BGFP-CM-week2.5-sham_Rep1.final.bam,5_H2BGFP-CM-week2.5-sham_Rep2.final.bam -b 5_Input-CM-week2.5-TAC_Rep1.final.bam -c 10000000 --extend 74 -o ./ >danpos.log 2>danpos.err

##1.1.2 Histone modification/TF/chromatin reulators ChIP-seq
cd /rd1/user/liym/nucleosomeTurnover/ChIP-seq
###1.1.2.1 TFs
cd TF
cd GSE29636_JCI2011_SE36 #Tbx20
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR223/SRR223492/SRR223492.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR223/SRR223493/SRR223493.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR223/SRR223494/SRR223494.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR223/SRR223495/SRR223495.fastq.gz
fqSeFilter.pl -r 30 1_Input-SRR223494-wholeHeart-2months_Rep2.fastq > 1_Input-SRR223494-wholeHeart-2months_Rep2.filter.fastq 2> fqFilter_input_Rep2.log
ls *.fastq|sed 's/.fastq//'|while read file;do 
    bwa aln -t 3 /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}.fastq >${file}.out.sai 2>log/${file}.aln.log;
    bwa samse /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}.out.sai ${file}.fastq 2>log/${file}.samse.log |samtools view -bSu -| /home/share/local/bin/samtools sort -@ 3 -o ${file}.out.sorted.bam - ;
    bamtools filter -in ${file}.out.sorted.bam -out ${file}.uniq.sorted.bam -script /mnt/share/liym/data/bwa/bwaAln.filter.json 2>log/${file}.bamtoolsFilter.log
    samtools rmdup -s ${file}.uniq.sorted.bam ${file}.final.rmDup.bam 2>log/${file}.rmDup.log
    done;
####Call peaks
macs2 callpeak -t 1_SRR223492-Tbx20-wholeHeart-2months_Rep1.final.rmDup.bam -c 1_Input-SRR223493-wholeHeart-2months_Rep1.final.rmDup.bam 1_Input-SRR223494-wholeHeart-2months_Rep2.filter.final.rmDup.bam 1_Input-SRR223495-wholeHeart-2months_Rep3.final.rmDup.bam -f BAM -g mm --outdir Tbx20 -n Tbx20 -B --SPMR >log/1_SRR223492-Tbx20-wholeHeart-2months_Rep1.macs2.log 2>log/1_SRR223492-Tbx20-wholeHeart-2months_Rep1.macs2.err
macs2 callpeak -t 1_SRR223492-Tbx20-wholeHeart-2months_Rep1.final.rmDup.bam -c 1_Input-SRR223493-wholeHeart-2months_Rep1.final.rmDup.bam -f BAM -g mm --outdir Tbx20_rep1 -n Tbx20 >log/1_SRR223492-Tbx20-wholeHeart-2months_Rep1.macs2.log 2>log/1_SRR223492-Tbx20-wholeHeart-2months_Rep1.macs2.err
macs2 callpeak -t 1_SRR223492-Tbx20-wholeHeart-2months_Rep1.final.rmDup.bam -c 1_Input-SRR223495-wholeHeart-2months_Rep3.final.rmDup.bam -f BAM -g mm --outdir Tbx20_rep3 -n Tbx20 >log/1_SRR223492-Tbx20-wholeHeart-2months_Rep3.macs2.log 2>log/1_SRR223492-Tbx20-wholeHeart-2months_Rep3.macs2.err
macs2 callpeak -t 1_SRR223492-Tbx20-wholeHeart-2months_Rep1.final.rmDup.bam -c 1_Input-SRR223494-wholeHeart-2months_Rep2.filter.final.rmDup.bam -f BAM -g mm --outdir Tbx20_rep2 -n Tbx20 >log/1_SRR223492-Tbx20-wholeHeart-2months_Rep2.macs2.log 2>log/1_SRR223492-Tbx20-wholeHeart-2months_Rep2.macs2.err
cat <(cut -f1-3 Tbx20_rep1/Tbx20_peaks.narrowPeak) <(cut -f1-3 Tbx20_rep2/Tbx20_peaks.narrowPeak) <(cut -f1-3 Tbx20_rep3/Tbx20_peaks.narrowPeak)|sort -k1,1 -k2,2n |bedtools merge -i stdin >Tbx20_combined.peaks.bed3
####Generate signals
Rscript /mnt/share/liym/tools/phantompeakqualtools/run_spp.R -c=input.merged.final.bam -s=0:5:1500 -savp=input.estFrag.pdf -out=input.spp.estFrag.rst
Rscript /mnt/share/liym/tools/phantompeakqualtools/run_spp.R -c=1_SRR223492-Tbx20-wholeHeart-2months_Rep1.final.rmDup.bam -s=0:5:1500 -savp=Tbx20.estFrag.pdf -out=Tbx20.spp.estFrag.rst
samtools merge -@ 10 input.merged.final.bam 1_Input-SRR22349*final.rmDup.bam 
input=$(grep "Total" input.merged.final.bamStats|awk 'BEGIN{sum=0}{sum+=$3}END{print 10000000/sum}')
factor=$(grep "Total" 1_SRR223492-Tbx20-wholeHeart-2months_Rep1.final.rmDup.bamStats|awk 'BEGIN{sum=0}{sum+=$3}END{print 10000000/sum}')
/mnt/share/liym/tools/deepTools-2.0.0/bin/bamCompare -b1 1_SRR223492-Tbx20-wholeHeart-2months_Rep1.final.rmDup.bam -b2 input.merged.final.bam --scaleFactors ${factor}:$input --ratio subtract -bs 20 -p 10 --extendReads 70 -o Tbx20.subtract.bw --outFileFormat bigwig >log/bamCompare.log 2>log/bamCompare.err
cd GSE35151_adultHeart #Tbx3 and Nkx2.5
#Input
fqTrimer.pl -r [ATGCN]{14} 1_Input-SRR400044-Heart-control_Rep1.fastq >1_Input-SRR400044-Heart-control_Rep1.trim.fq
fqSeFilter.pl -f 0.1 -b N -l 30 -r 30 -q 0.5 -c 20 -a 20 1_Input-SRR400044-Heart-control_Rep1.trim.fq > 1_Input-SRR400044-Heart-control_Rep1.filter.fq 2> input.fqFilter.log
fqAdapterFilter.pl -a ATCGGAAGAGCTCGTATGCCGTCTTCTGCTTAGAT,ATCGGAAGAGCTCGTATGCCGTCTTCTGCTTGAAA,ATCGGAAGAGCTCGTATGCCGTCTTCTGCTTATAT -i 1_Input-SRR400044-Heart-control_Rep1.filter.fq -o 1_Input-SRR400044-Heart-control_Rep1.filter.adaF.fastq -r 1_Input-SRR400044-Heart-control_Rep1.filter.adaF.log
#Gata4
fqAdapterFilter.pl -a GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTAGA,GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTGGA,GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTGAA -i 1_SRR400045-Heart-Gata4_Rep1.fastq -o 1_SRR400045-Heart-Gata4_Rep1.filter.fastq -r 1_SRR400045-Heart-Gata4_Rep1.fqFilter.log
#Nkx2.5
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTGAACATCTCGTATGCC,GGTGTTGTTGTTGTCTTAGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG,GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTAGATCGGAAGAGCTCGTAT -i 1_SRR400046-Heart-Nkx2-5_Rep1.fastq -o 1_SRR400046-Heart-Nkx2-5_Rep1.filter.fastq -r 1_SRR400046-Heart-Nkx2-5_Rep1.adaFilter.log
#Tbx3
fqTrimer.pl -r [ATGCN]{14} 1_SRR400043-Heart-Tbx3_Rep1.fastq >1_SRR400043-Heart-Tbx3_Rep1.trim.fastq 2>1_SRR400043-Heart-Tbx3_Rep1.fqTrim.log
fqSeFilter.pl -f 0.1 -b N -l 30 -r 30 -q 0.5 -c 20 -a 20 -t 10,5 1_SRR400043-Heart-Tbx3_Rep1.trim.fastq > 1_SRR400043-Heart-Tbx3_Rep1.filter.fastq 2>1_SRR400043-Heart-Tbx3_Rep1.fqFilter.log
ls *filter*fastq|sed 's/.fastq//'|while read file;do
    bwa aln -t 5 /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}.fastq >${file}.out.sai 2>log/${file}.aln.log;
    bwa samse /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}.out.sai ${file}.fastq 2>log/${file}.samse.log |samtools view -bSu -| /home/share/local/bin/samtools sort -@ 3 -o ${file}.out.sorted.bam - ;
    bamtools filter -in ${file}.out.sorted.bam -out ${file}.uniq.sorted.bam -script /mnt/share/liym/data/bwa/bwaAln.filter.json 2>log/${file}.bamtoolsFilter.log
    samtools rmdup -s ${file}.uniq.sorted.bam ${file}.final.rmDup.bam 2>log/${file}.rmDup.log
done;
macs2 callpeak -t 1_SRR400043-Heart-Tbx3_Rep1.filter.final.rmDup.bam -c 1_Input-SRR400044-Heart-control_Rep1.filter.adaF.final.rmDup.bam -f BAM -g mm --outdir Tbx3 -n Tbx3 -B --SPMR >log/1_SRR400043-Heart-Tbx3_Rep1.macs2.log 2>log/1_SRR400043-Heart-Tbx3_Rep1.macs2.err
macs2 callpeak -t 1_SRR400045-Heart-Gata4_Rep1.filter.final.rmDup.bam -c 1_Input-SRR400044-Heart-control_Rep1.filter.adaF.final.rmDup.bam -f BAM -g mm --outdir Gata4 -n Gata4 -B --SPMR >log/1_SRR400045-Heart-Gata4_Rep1.macs2.log 2>log/1_SRR400045-Heart-Gata4_Rep1.macs2.err &
macs2 callpeak -t 1_SRR400046-Heart-Nkx2-5_Rep1.filter.final.rmDup.bam -c 1_Input-SRR400044-Heart-control_Rep1.filter.adaF.final.rmDup.bam -f BAM -g mm --outdir Nkx2-5 -n Nkx2-5 -B --SPMR >log/1_SRR400046-Heart-Nkx2-5_Rep1.macs2.log 2>log/1_SRR400046-Heart-Nkx2-5_Rep1.macs2.err
Rscript /mnt/share/liym/tools/phantompeakqualtools/run_spp.R -c=1_SRR400043-Heart-Tbx3_Rep1.filter.final.rmDup.bam -s=0:5:1500 -savp=Tbx3.estFrag.pdf -out=Tbx3.spp.estFrag.rst > log/Tbx3.spp.log 2> log/Tbx3.spp.err
Rscript /mnt/share/liym/tools/phantompeakqualtools/run_spp.R -c=1_SRR400046-Heart-Nkx2-5_Rep1.filter.final.rmDup.bam -s=0:5:1500 -savp=Nkx2-5.estFrag.pdf -out=Nkx2-5.spp.estFrag.rst > log/Nkx2-5.spp.log 2> log/Nkx2-5.spp.err 
input=$(grep "Total" 1_Input-SRR400044-Heart-control_Rep1.filter.adaF.final.rmDup.bamStats|awk '{print 10000000/$3}')
tbx3=$(grep "Total" 1_SRR400043-Heart-Tbx3_Rep1.filter.final.rmDup.bamStats|awk '{print 10000000/$3}')
nkx25=$(grep "Total" 1_SRR400046-Heart-Nkx2-5_Rep1.filter.final.rmDup.bamStats|awk '{print 10000000/$3}')
/mnt/share/liym/tools/deepTools-2.0.0/bin/bamCompare -b1 1_SRR400043-Heart-Tbx3_Rep1.filter.final.rmDup.bam -b2 1_Input-SRR400044-Heart-control_Rep1.filter.adaF.final.rmDup.bam --scaleFactors ${tbx3}:$input --ratio subtract -bs 20 -p 10 --extendReads 109 -o Tbx3.subtract.bw --outFileFormat bigwig >log/bamCompare.log 2>log/bamCompare.err
/mnt/share/liym/tools/deepTools-2.0.0/bin/bamCompare -b1 1_SRR400046-Heart-Nkx2-5_Rep1.filter.final.rmDup.bamStats -b2 1_Input-SRR400044-Heart-control_Rep1.filter.adaF.final.rmDup.bam --scaleFactors ${nkx25}:$input --ratio subtract -bs 20 -p 10 --extendReads 97 -o Nkx2-5.subtract.bw --outFileFormat bigwig >log/bamCompare.log 2>log/bamCompare.err
cd GSE52123-2014NC #Gata4 and pol2
perl /mnt/share/liym/bin/fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGC -i 1_Pol2-WT-VentricularApex-Adult-CM_Rep1.fastq -o 1_Pol2-WT-VentricularApex-Adult-CM_Rep1.filter.fastq -r 1_Pol2-WT-VentricularApex-Adult-CM_Rep1.fqAdaFilter.log
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGC,GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAAAATCTCGTATGC,CGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGT -i 1_Pol2-WT-VentricularApex-Adult-CM_Rep2.fastq -o 1_Pol2-WT-VentricularApex-Adult-CM_Rep2.filter.fastq -r 1_Pol2-WT-VentricularApex-Adult-CM_Rep2.fqAdaFilter.log
ls *.fastq|sed 's/.fastq//'|while read file;do 
    bwa aln -t 5 /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}.fastq >${file}.out.sai 2>log/${file}.aln.log;
    bwa samse /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}.out.sai ${file}.fastq 2>log/${file}.samse.log |samtools view -bSu -| /home/share/local/bin/samtools sort -@ 5 -o ${file}.out.sorted.bam - ;
    bamtools filter -in ${file}.out.sorted.bam -out ${file}.uniq.sorted.bam -script /mnt/share/liym/data/bwa/bwaAln.filter.json 2>log/${file}.bamtoolsFilter.log
    samtools rmdup -s ${file}.uniq.sorted.bam ${file}.final.rmDup.bam 2>log/${file}.rmDup.log
    done;
macs2 callpeak -t 1_Pol2-WT-VentricularApex-Adult-CM_Rep2.filter.final.rmDup.bam -c 1_Input-Adult-histone-VentricularApex-CM_Rep1.final.rmDup.bam -f BAM -g mm --outdir Pol2_Rep2 -n Pol2 -B --SPMR >log/1_Pol2-WT-VentricularApex-Adult-CM.Rep2.macs2.log 2>log/1_Pol2-WT-VentricularApex-Adult-CM.Rep2.macs2.err &
macs2 callpeak -t 1_Pol2-WT-VentricularApex-Adult-CM_Rep1.filter.final.rmDup.bam -c 1_Input-Adult-histone-VentricularApex-CM_Rep1.final.rmDup.bam -f BAM -g mm --outdir Pol2_Rep1 -n Pol2 -B --SPMR >log/1_Pol2-WT-VentricularApex-Adult-CM.Rep1.macs2.log 2>log/1_Pol2-WT-VentricularApex-Adult-CM.Rep1.macs2.err &
macs2 callpeak -t 2_GATA4-fb-Adult-SRR1025222-SRR1025223_Rep2.final.rmDup.bam -c 2_Input-fb-Ab-Adult-SRR1025224_Rep1.final.rmDup.bam -f BAM -g mm --outdir Gata4_Rep2 -n Gata4 -B --SPMR >log/2_GATA4-fb-Adult-SRR1025222-SRR1025223.Rep2.macs2.log 2>log/2_GATA4-fb-Adult-SRR1025222-SRR1025223.Rep2.macs2.err &
macs2 callpeak -t 2_GATA4-fb-Adult-SRR1025222-SRR1025223_Rep1.final.rmDup.bam -c 2_Input-fb-Ab-Adult-SRR1025224_Rep1.final.rmDup.bam -f BAM -g mm --outdir Gata4_Rep1 -n Gata4 -B --SPMR >log/2_GATA4-fb-Adult-SRR1025222-SRR1025223.Rep1.macs2.log 2>log/2_GATA4-fb-Adult-SRR1025222-SRR1025223.Rep1.macs2.err &
mkdir result && cd result
bedtools intersect -wo -a ../Pol2_Rep1/Pol2_peaks.narrowPeak -b ../Pol2_Rep2/Pol2_peaks.narrowPeak|awk -v OFS="\t" -v ORS="" '{print $1"\t";if($2<=$12){print $2"\t";}else{print $12"\t";}if($3<=$13){print $13"\n"}else{print $3"\n"}}'|sort|uniq|sort -k1,1 -k2,2n |bedtools merge -i stdin >pol2.intersect.merged.peaks.bed
bedtools intersect -wo -a ../Gata4_Rep1/Gata4_peaks.narrowPeak -b ../Gata4_Rep2/Gata4_peaks.narrowPeak| awk -v OFS="\t" -v ORS="" '{print $1"\t";if($2<=$12){print $2"\t";}else{print $12"\t";}if($3<=$13){print $13"\n"}else{print $3"\n"}}'|sort|uniq |grep -vE "M|Un|random" >Gata4.intersect.merged.peaks.v1.bed 
cd ENCODE #P300
fqAdapterFilter.pl -a GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTGAAA -i Heart-Input-adult-8wks_Rep1.SRR489733.fq -o Heart-Input-adult-8wks_Rep1.SRR489733.adaF.fq -r Heart-Input-adult-8wks_Rep1.SRR489733.adaF.log 
fqSeFilter.pl -f 0.1 -l 30 -q 0.5 -c 10 -a 20 Heart-Input-adult-8wks_Rep1.SRR489733.adaF.fq >Heart-Input-adult-8wks_Rep1.SRR489733.Filter.fq 2>Heart-Input-adult-8wks_Rep1.SRR489733.fqFilter.log 
fqSeFilter.pl -f 0.1 -l 30 -q 0.5 -c 10 -a 20 Heart-Input-adult-8wks_Rep2.SRR489734.fq >Heart-Input-adult-8wks_Rep2.SRR489734.Filter.fq 2>Heart-Input-adult-8wks_Rep2.SRR489734.fqFilter.log
fqSeFilter.pl -f 0.1 -l 30 -q 0.5 -c 10 -a 20 -r 20 Heart-P300-adult-8wks_Rep1.SRR489717.fq >Heart-P300-adult-8wks_Rep1.SRR489717.Filter.fq 2>Heart-P300-adult-8wks_Rep1.SRR489717.fqFilter.log
fqAdapterFilter.pl -a GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTGAAA -i Heart-P300-adult-8wks_Rep2.SRR489718.fq -o Heart-P300-adult-8wks_Rep2.SRR489718.adaF.fq -r Heart-P300-adult-8wks_Rep2.SRR489718.adaF.log
fqSeFilter.pl -f 0.1 -l 30 -q 0.5 -c 10 -a 20 Heart-P300-adult-8wks_Rep2.SRR489718.adaF.fq >Heart-P300-adult-8wks_Rep2.SRR489718.Filter.fq 2>Heart-P300-adult-8wks_Rep2.SRR489718.fqFilter.log
ls *.fq|sed 's/.fq//'|while read file;do 
    bwa aln -t 5 /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}.fq >${file}.out.sai 2>${file}.aln.log;
    bwa samse /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}.out.sai ${file}.fq 2>${file}.samse.log |samtools view -bSu -| /home/share/local/bin/samtools sort -@ 3 -o ${file}.out.sorted.bam - ;
    bamtools filter -in ${file}.out.sorted.bam -out ${file}.uniq.sorted.bam -script /mnt/share/liym/data/bwa/bwaAln.filter.json 2>${file}.bamtoolsFilter.log
    samtools rmdup -s ${file}.uniq.sorted.bam ${file}.final.rmDup.bam 2>${file}.rmDup.log
done;
macs2 callpeak -t Heart-P300-adult-8wks_Rep1.SRR489717.Filter.final.rmDup.bam -c Heart-Input-adult-8wks_Rep1.SRR489733.Filter.final.rmDup.bam -f BAM -g mm --outdir p300_rep1 -n p300 >Heart-P300-adult-8wks_Rep1.SRR489717.macs2.log 2>Heart-P300-adult-8wks_Rep1.SRR489717.macs2.err
macs2 callpeak -t Heart-P300-adult-8wks_Rep2.SRR489718.Filter.final.rmDup.bam -c Heart-Input-adult-8wks_Rep2.SRR489734.Filter.final.rmDup.bam -f BAM -g mm --outdir p300_rep2 -n p300 >Heart-P300-adult-8wks_Rep2.SRR489718.log 2>Heart-P300-adult-8wks_Rep2.SRR489718.err
###1.1.2.2 Histone modification
cd histone
cd EEDproject
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGC -i 1_H3K27me3-EEDko-AdultCM_Rep1.fastq -o 1_H3K27me3-EEDko-AdultCM_Rep1.adaF.fastq -r log/1_H3K27me3-EEDko-AdultCM_Rep1.adaF.log &
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGC -i 1_H3K27me3-EEDko-AdultCM_Rep2.fastq -o 1_H3K27me3-EEDko-AdultCM_Rep2.adaF.fastq -r log/1_H3K27me3-EEDko-AdultCM_Rep2.adaF.log &
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGC -i 1_H3K27me3-WT-AdultCM_Rep1.fastq -o 1_H3K27me3-WT-AdultCM_Rep1.adaF.fastq -r log/1_H3K27me3-WT-AdultCM_Rep1.adaF.log &
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGC -i 1_H3K27me3-WT-AdultCM_Rep2.fastq -o 1_H3K27me3-WT-AdultCM_Rep2.adaF.fastq -r log/1_H3K27me3-WT-AdultCM_Rep2.adaF.log &
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGC -i 1_H3K27me3-WT-AdultCM_Rep3.fastq -o 1_H3K27me3-WT-AdultCM_Rep3.adaF.fastq -r log/1_H3K27me3-WT-AdultCM_Rep3.adaF.log &
ls *.fastq|sed 's/.fastq//'|while read file;do 
    bwa aln -t 10 /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}.fastq >${file}.out.sai 2>log/${file}.aln.log;
    bwa samse /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}.out.sai ${file}.fastq 2>log/${file}.samse.log |samtools view -bSu -| /home/share/local/bin/samtools sort -@ 5 -o ${file}.out.sorted.bam - ;
    bamtools filter -in ${file}.out.sorted.bam -out ${file}.uniq.sorted.bam -script /mnt/share/liym/data/bwa/bwaAln.filter.json 2>log/${file}.bamtoolsFilter.log
    samtools rmdup -s ${file}.uniq.sorted.bam ${file}.final.rmDup.bam 2>log/${file}.rmDup.log
done;
ls 1_H3K27ac*final.rmDup.bam|while read file;do
        prefix=$(echo $file|sed 's/.final.rmDup.bam//');
    macs2 callpeak -t $file -c 1_Input-H3K27ac-AdultCM_Rep1.final.rmDup.bam -f BAM -g mm --broad --outdir $prefix -n $prefix  >log/${prefix}.log 2>log/${prefix}.err
done;
ls 1_H3K27me3*final.rmDup.bam|while read file;do
        prefix=$(echo $file|sed 's/.final.rmDup.bam//');
    macs2 callpeak -t $file -c 1_Input-AdultCM_Rep1.final.rmDup.bam -f BAM -g mm --broad --outdir ${prefix}_broad -n $prefix -B --SPMR >log/${prefix}_broad.log 2>log/${prefix}_broad.err
done;
####generate signals
ls *final.rmDup.bam|while read file;do prefix=$(echo $file|sed 's/.final.rmDup.bam//');Rscript /mnt/share/liym/tools/phantompeakqualtools/run_spp.R -c=$file -s=0:5:1500 -savp=${prefix}.estFrag.pdf -out=${prefix}.spp.estFrag.rst & done;
ls *.final.rmDup.bam|while read file;do samtools index $file;done;
ls *H3K27me3*final.rmDup.bam|while read file;do
    prefix=$(echo $file|sed 's/.final.rmDup.bam//')
    input=$(grep "Total" 1_Input-AdultCM_Rep1.final.rmDup.bamStats|awk 'BEGIN{sum=0}{sum+=$3}END{print 10000000/sum}')
    factor=$(grep "Total" ${file}Stats|awk 'BEGIN{sum=0}{sum+=$3}END{print 10000000/sum}')
    /mnt/share/liym/tools/deepTools-2.0.0/bin/bamCompare -b1 $file -b2 1_Input-AdultCM_Rep1.final.rmDup.bam --scaleFactors ${factor}:$input --ratio subtract -bs 20 -p 15 --extendReads 140 -o ${prefix}.subtract.bw --outFileFormat bigwig >log/bamCompare.log 2>log/bamCompare.err
    done;
ls *H3K27ac*final.rmDup.bam|while read file;do
    prefix=$(echo $file|sed 's/.final.rmDup.bam//')
    input=$(grep "Total" 1_Input-H3K27ac-AdultCM_Rep1.final.rmDup.bamStats|awk 'BEGIN{sum=0}{sum+=$3}END{print 10000000/sum}')
    factor=$(grep "Total" ${file}Stats|awk 'BEGIN{sum=0}{sum+=$3}END{print 10000000/sum}')
    /mnt/share/liym/tools/deepTools-2.0.0/bin/bamCompare -b1 $file -b2 1_Input-H3K27ac-AdultCM_Rep1.final.rmDup.bam --scaleFactors ${factor}:$input --ratio subtract -bs 20 -p 15 --extendReads 140 -o ${prefix}.subtract.bw --outFileFormat bigwig >log/bamCompare.log 2>log/bamCompare.err
    done;
mkdir result && cd result
####Final peaks
bedtools intersect -a ../1_H3K27ac-EEDko-AdultCM_Rep1/1_H3K27ac-EEDko-AdultCM_Rep1_peaks.broadPeak -b ../1_H3K27ac-EEDko-AdultCM_Rep2/1_H3K27ac-EEDko-AdultCM_Rep2_peaks.broadPeak -wo |awk -v OFS="\t" -v ORS="" '{print $1"\t";if($2<=$11){print $2"\t";}else{print $11"\t";}if($3<=$12){print $12"\n"}else{print $3"\n"}}'|sort|uniq |sort -k1,1 -k2,2n|bedtools merge >H3K27ac-EEDko.intersect.broadPeaks.bed 
bedtools intersect -a ../1_H3K27ac-WT-AdultCM_Rep1/1_H3K27ac-WT-AdultCM_Rep1_peaks.broadPeak -b ../1_H3K27ac-WT-AdultCM_Rep2/1_H3K27ac-WT-AdultCM_Rep2_peaks.broadPeak -wo |awk -v OFS="\t" -v ORS="" '{print $1"\t";if($2<=$11){print $2"\t";}else{print $11"\t";}if($3<=$12){print $12"\n"}else{print $3"\n"}}'|sort|uniq |sort -k1,1 -k2,2n|bedtools merge >H3K27ac-WT.intersect.broadPeaks.bed
bedtools intersect -a ../1_H3K27me3-EEDko-AdultCM_Rep1.adaF_broad/1_H3K27me3-EEDko-AdultCM_Rep1.adaF_peaks.broadPeak -b ../1_H3K27me3-EEDko-AdultCM_Rep2.adaF_broad/1_H3K27me3-EEDko-AdultCM_Rep2.adaF_peaks.broadPeak -wo |awk -v OFS="\t" -v ORS="" '{print $1"\t";if($2<=$11){print $2"\t";}else{print $11"\t";}if($3<=$12){print $12"\n"}else{print $3"\n"}}'|sort|uniq |sort -k1,1 -k2,2n|bedtools merge >H3K27me3-EEDko.intersect.broadPeaks.bed
bedtools intersect -wo -a ../1_H3K27me3-WT-AdultCM_Rep1.adaF_broad/1_H3K27me3-WT-AdultCM_Rep1.adaF_peaks.broadPeak -b ../1_H3K27me3-WT-AdultCM_Rep2.adaF_broad/1_H3K27me3-WT-AdultCM_Rep2.adaF_peaks.broadPeak|awk -v OFS="\t" -v ORS="" '{print $1"\t";if($2<=$11){print $2"\t";}else{print $11"\t";}if($3<=$12){print $12"\n"}else{print $3"\n"}}'|sort|uniq|sort -k1,1 -k2,2n |bedtools merge|bedtools intersect -wo -a stdin -b ../1_H3K27me3-WT-AdultCM_Rep3.adaF_broad/1_H3K27me3-WT-AdultCM_Rep3.adaF_peaks.broadPeak|awk -v OFS="\t" -v ORS="" '{print $1"\t";if($2<=$5){print $2"\t";}else{print $5"\t";}if($3<=$6){print $6"\n"}else{print $3"\n"}}'|sort -k1,1 -k2,2n|uniq |bedtools merge >H3K27me3-WT.intersect.broadPeaks.bed
ls *broadPeaks.bed |while read file;do
    awk '$1!~"random|GL|Un|chrM"' $file >tmp;
    mv tmp $file;
    done;
cd ../banding-sham
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGC -i 1_CM-banding-H3K4me3-CHS00010885-GCCAAT_Rep1.fastq -o 1_CM-banding-H3K4me3-CHS00010885-GCCAAT_Rep1.filter.fastq -r log/1_CM-banding-H3K4me3-CHS00010885-GCCAAT_Rep1.fqFilter.log 
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGC -i 1_CM-banding-H3K9me3-CHS00010887-ACTTGA_Rep1.fastq -o 1_CM-banding-H3K9me3-CHS00010887-ACTTGA_Rep1.filter.fastq -r log/1_CM-banding-H3K9me3-CHS00010887-ACTTGA_Rep1.log
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGC -i 1_CM-banding-Pol2-CHS00010886-CAGATC_Rep1.fastq -o 1_CM-banding-Pol2-CHS00010886-CAGATC_Rep1.filter.fastq -r log/1_CM-banding-Pol2-CHS00010886-CAGATC_Rep1.fqFilter.log
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGC -i 1_CM-sham-H3K4me1-CHS00010890-CTTGTA_Rep1.fastq -o 1_CM-sham-H3K4me1-CHS00010890-CTTGTA_Rep1.filter.fastq -r log/1_CM-sham-H3K4me1-CHS00010890-CTTGTA_Rep1.fqFilter.log 
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGC -i 1_CM-sham-H3K4me3-CHS00010891-TGACCA_Rep1.fastq -o 1_CM-sham-H3K4me3-CHS00010891-TGACCA_Rep1.filter.fastq -r log/1_CM-sham-H3K4me3-CHS00010891-TGACCA_Rep1.fqFilter.log
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGC -i 1_CM-sham-H3K9me3-CHS00010893-GCCAAT_Rep1.fastq -o 1_CM-sham-H3K9me3-CHS00010893-GCCAAT_Rep1.filter.fastq -r log/1_CM-sham-H3K9me3-CHS00010893-GCCAAT_Rep1.fqFilter.log
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGC -i 1_CM-sham-H3K27ac-CHS00010889-TAGCTT_Rep1.fastq -o 1_CM-sham-H3K27ac-CHS00010889-TAGCTT_Rep1.filter.fastq -r log/1_CM-sham-H3K27ac-CHS00010889-TAGCTT_Rep1.fqFilter.log
fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGC -i 1_CM-sham-Pol2-CHS00010892-ACAGTG_Rep1.fastq -o 1_CM-sham-Pol2-CHS00010892-ACAGTG_Rep1.filter.fastq -r log/1_CM-sham-Pol2-CHS00010892-ACAGTG_Rep1.fqFilter.log
fqSeFilter.pl -f 0.1 -b N -q 0.5 -c 10 -a 20 H3K27me3_Heart_sham_ERR231658_Rep1.fastq >H3K27me3_Heart_sham_ERR231658_Rep1.filter.fastq 2>log/H3K27me3_Heart_sham_ERR231658_Rep1.filter.log 
fqSeFilter.pl -f 0.1 -b N -q 0.5 -c 10 -a 20 H3K27me3_Heart_TAC_ERR231655_Rep1.fastq >H3K27me3_Heart_TAC_ERR231655_Rep1.filter.fastq 2>log/H3K27me3_Heart_TAC_ERR231655_Rep1.filter.log 
fqSeFilter.pl -f 0.1 -b N -q 0.5 -c 10 -a 20 Input_Heart_sham_ERR231653.fastq >Input_Heart_sham_ERR231653.filter.fastq 2>log/Input_Heart_sham_ERR231653.filter.log 
fqSeFilter.pl -f 0.1 -b N -q 0.5 -c 10 -a 20 Input_Heart_TAC_ERR231657.fastq >Input_Heart_TAC_ERR231657.filter.fastq 2>log/Input_Heart_TAC_ERR231657.filter.log
#####mapping
ls *.fastq|sed 's/.fastq//'|while read file;do 
    bwa aln -t 10 /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}.fastq >${file}.out.sai 2>log/${file}.aln.log;
    bwa samse /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}.out.sai ${file}.fastq 2>log/${file}.samse.log |samtools view -bSu -| /home/share/local/bin/samtools sort -@ 3 -o ${file}.out.sorted.bam - ;
    bamtools filter -in ${file}.out.sorted.bam -out ${file}.uniq.sorted.bam -script /mnt/share/liym/data/bwa/bwaAln.filter.json 2>log/${file}.bamtoolsFilter.log
    samtools rmdup -s ${file}.uniq.sorted.bam ${file}.final.rmDup.bam 2>log/${file}.rmDup.log
done;
ls 1_CM-banding*final.rmDup.bam|while read file;do
    prefix=$(echo $file|cut -f1-3 -d '-')
    macs2 callpeak -t $file -c 1_Input-CM-banding-CHS00010888-GATCAG_Rep1.final.rmDup.bam -f BAM -g mm --outdir ${prefix} -n $prefix -B --SPMR >log/${prefix}_macs2.log 2>log/${prefix}_macs2.err
done;
ls 1_CM-sham*final.rmDup.bam|while read file;do
        name=$(echo $file|cut -f1-3 -d '-')
    macs2 callpeak -t $file -c 1_Input-CM-sham-CHS00010894-CAGATC_Rep1.final.rmDup.bam -f BAM -g mm --outdir $name -n $name -B --SPMR >log/$name_macs2.log 2>log/$name_macs2.err
done;   
#####H3K27ac broad peaks 
macs2 callpeak -t 1_CM-banding-H3K27ac-CHS00010883-TGACCA_Rep1.final.rmDup.bam -c 1_Input-CM-banding-CHS00010888-GATCAG_Rep1.final.rmDup.bam -f BAM -g mm --outdir 1_CM-banding-H3K27ac_broad --broad >log/1_CM-banding-H3K27ac_broad.macs2.log 2>log/1_CM-banding-H3K27acbroad._macs2.err
macs2 callpeak -t 1_CM-sham-H3K27ac-CHS00010889-TAGCTT_Rep1.filter.final.rmDup.bam -c 1_Input-CM-sham-CHS00010894-CAGATC_Rep1.final.rmDup.bam -f BAM -g mm --outdir 1_CM-sham-H3K27ac_broad --broad >log/1_CM-sham-H3K27ac_broad.macs2.log 2>log/1_CM-sham-H3K27ac_broad.macs2.err
#####H3K9me3 broad peaks 
macs2 callpeak -t 1_CM-sham-H3K9me3-CHS00010893-GCCAAT_Rep1.filter.final.rmDup.bam -c 1_Input-CM-sham-CHS00010894-CAGATC_Rep1.final.rmDup.bam -f BAM -g mm --outdir 1_CM-sham-H3K9me3 -n 1_CM-sham-H3K9me3 -B --SPMR --qvalue 0.05 --broad >log/1_CM-sham-H3K9me3_macs2.log 2>log/1_CM-sham-H3K9me3_macs2.err
macs2 callpeak -t 1_CM-banding-H3K9me3-CHS00010887-ACTTGA_Rep1.filter.final.rmDup.bam -c 1_Input-CM-banding-CHS00010888-GATCAG_Rep1.final.rmDup.bam -f BAM -g mm --outdir 1_CM-banding-H3K9me3 -n 1_CM-banding-H3K9me3 -B --SPMR --qvalue 0.05 --broad >log/1_CM-banding-H3K9me3_macs2.log 2>log/1_CM-banding-H3K9me3_macs2.err
macs2 callpeak -t 1_CM-sham-H3K9me3-CHS00010893-GCCAAT_Rep1.filter.final.rmDup.bam -c 1_Input-CM-sham-CHS00010894-CAGATC_Rep1.final.rmDup.bam -f BAM -g mm --outdir 1_CM-sham-H3K9me3_v2 -n 1_CM-sham-H3K9me3 --pvalue 1e-5 --broad-cutoff 1e-3 --broad >log/1_CM-sham-H3K9me3_macs2_v2.log 2>log/1_CM-sham-H3K9me3_macs2_v2.err
macs2 callpeak -t 1_CM-banding-H3K9me3-CHS00010887-ACTTGA_Rep1.filter.final.rmDup.bam -c 1_Input-CM-banding-CHS00010888-GATCAG_Rep1.final.rmDup.bam -f BAM -g mm --outdir 1_CM-banding-H3K9me3_v2 -n 1_CM-banding-H3K9me3 --pvalue 1e-5 --broad-cutoff 1e-3 --broad >log/1_CM-banding-H3K9me3_macs2_v2.log 2>log/1_CM-banding-H3K9me3_macs2_v2.err
#####H3K27me3 broad peaks
macs2 callpeak -t H3K27me3_Heart_TAC_ERR231655.final.rmDup.bam -c Input_Heart_TAC_ERR231657.final.rmDup.bam --broad-cutoff 1e-3 -f BAM -g mm --outdir H3K27me3_Heart_TAC -n H3K27me3_Heart_TAC --broad >log/H3K27me3_Heart_TAC_broad.macs2.log 2>log/H3K27me3_Heart_TAC_broad.macs2.err
macs2 callpeak -t H3K27me3_Heart_sham_ERR231658.final.rmDup.bam -c Input_Heart_sham_ERR231653.final.rmDup.bam --broad-cutoff 1e-3 -f BAM -g mm --outdir H3K27me3_Heart_sham -nH3K27me3_Heart_sham --broad >log/H3K27me3_Heart_sham_broad.macs2.log 2>log/H3K27me3_Heart_sham_broad.macs2.err
/home/liym/.local/bin/macs14 callpeak -t ../H3K27me3_Heart_TAC_ERR231655.final.rmDup.bam -c ../Input_Heart_TAC_ERR231657.final.rmDup.bam -p 1e-3 -f BAM -g mm -n H3K27me3_Heart_TAC >../log/H3K27me3_Heart_TAC_macs14.log 2>../log/H3K27me3_Heart_TAC_macs14.err
/home/liym/.local/bin/macs14 callpeak -t ../H3K27me3_Heart_sham_ERR231658.final.rmDup.bam -c ../Input_Heart_sham_ERR231653.final.rmDup.bam -p 1e-3 -f BAM -g mm -n H3K27me3_Heart_sham >../log/H3K27me3_Heart_sham_macs14.log 2>../log/H3K27me3_Heart_sham_macs14.err
mkdir result && cd result
#####H3K4me1 & H3K4me3
myPeak(){
    bedtools intersect -wao -a ../1_CM-banding-${prefix}/1_CM-banding-${prefix}_peaks.narrowPeak -b ../1_CM-sham-${prefix}/1_CM-sham-${prefix}_peaks.narrowPeak >banding_sham.${prefix}.intersect.bed 
    bedtools intersect -wao -a ../1_CM-sham-${prefix}/1_CM-sham-${prefix}_peaks.narrowPeak -b ../1_CM-banding-${prefix}/1_CM-banding-${prefix}_peaks.narrowPeak >sham_banding.${prefix}.intersect.bed; 
    awk '$12<0' banding_sham.${prefix}.intersect.bed |cut -f1-10 >banding_specific.${prefix}.bed;
    awk '$12<0' sham_banding.${prefix}.intersect.bed |cut -f1-10 >sham_specific.${prefix}.bed;
    awk '$12>0' banding_sham.${prefix}.intersect.bed |awk '$7/$17>1.5 && $3-$2>=20'|awk -v OFS="\t" -v ORS="" 'BEGIN{i=1}{print $1"\t";if($2<=$12){print $12"\t"}else{print $2"\t"}if($3<=$13){print $3"\t"}else{print $13"\t"}print "banding_up_"i"\n";i+=1;}' >banding_up.${prefix}.bed
    awk '$12>0' banding_sham.${prefix}.intersect.bed |awk '$7/$17<0.5 && $3-$2>=20'|awk -v OFS="\t" -v ORS="" 'BEGIN{i=1}{print $1"\t";if($2<=$12){print $12"\t"}else{print $2"\t"}if($3<=$13){print $3"\t"}else{print $13"\t"}print "banding_down_"i"\n";i+=1;}' >banding_down.${prefix}.bed
    awk '$12>0' banding_sham.${prefix}.intersect.bed |awk '$7/$17>=0.5 && $7/$17<=1.5 && $3-$2>=20'|awk -v OFS="\t" -v ORS="" 'BEGIN{i=1}{print $1"\t";if($2<=$12){print $12"\t"}else{print $2"\t"}if($3<=$13){print $3"\t"}else{print $13"\t"}print "unchange_"i"\n";i+=1;}' >banding_sham.unchange.${prefix}.bed
}
prefix=H3K4me1
myPeak
prefix=H3K4me3
myPeak
#####H3K9me3, H3K27ac and H3K27me3 broad peaks
myBroadPeak(){
    awk '$11<0' banding_sham.${prefix}.intersect.bed |cut -f1-9>banding_specific.${prefix}.bed;
    awk '$11<0' sham_banding.${prefix}.intersect.bed |cut -f1-9 >sham_specific.${prefix}.bed;
    awk '$11>0' banding_sham.${prefix}.intersect.bed |awk '$7/$16>1.5 && $3-$2>=20'|awk -v OFS="\t" -v ORS="" 'BEGIN{i=1}{print $1"\t";if($2<=$11){print $11"\t"}else{print $2"\t"}if($3<=$12){print $3"\t"}else{print $12"\t"}print "banding_up_"i"\n";i+=1;}' >banding_up.${prefix}.bed
    awk '$11>0' banding_sham.${prefix}.intersect.bed |awk '$7/$16<0.5 && $3-$2>=20'|awk -v OFS="\t" -v ORS="" 'BEGIN{i=1}{print $1"\t";if($2<=$11){print $11"\t"}else{print $2"\t"}if($3<=$12){print $3"\t"}else{print $12"\t"}print "banding_down_"i"\n";i+=1;}' >banding_down.${prefix}.bed
    awk '$11>0' banding_sham.${prefix}.intersect.bed |awk '$7/$16>=0.5 && $7/$16<=1.5 && $3-$2>=20'|awk -v OFS="\t" -v ORS="" 'BEGIN{i=1}{print $1"\t";if($2<=$11){print $11"\t"}else{print $2"\t"}if($3<=$12){print $3"\t"}else{print $12"\t"}print "unchange_"i"\n";i+=1;}' >banding_sham.unchange.${prefix}.bed
}
prefix=H3K9me3
bedtools intersect -wao -a ../1_CM-banding-${prefix}_v2/1_CM-banding-${prefix}_peaks.broadPeak -b ../1_CM-sham-${prefix}/1_CM-sham-${prefix}_peaks.broadPeak >banding_sham.${prefix}.intersect.bed 
bedtools intersect -wao -a ../1_CM-sham-${prefix}_v2/1_CM-sham-${prefix}_peaks.broadPeak -b ../1_CM-banding-${prefix}/1_CM-banding-${prefix}_peaks.broadPeak >sham_banding.${prefix}.intersect.bed    
myBroadPeak
prefix=H3K27ac
bedtools intersect -wao -a ../1_CM-banding-${prefix}_broad/*_peaks.broadPeak -b ../1_CM-sham-${prefix}_broad/*_peaks.broadPeak >banding_sham.${prefix}.intersect.bed
bedtools intersect -wao -a ../1_CM-sham-${prefix}_broad/*_peaks.broadPeak -b ../1_CM-banding-${prefix}_broad/*_peaks.broadPeak >sham_banding.${prefix}.intersect.bed
myBroadPeak
bedtools intersect -wao -a ../H3K27me3_Heart_TAC_macs14/H3K27me3_Heart_TAC_peaks.FE.bed -b ../H3K27me3_Heart_sham_macs14/H3K27me3_Heart_sham_peaks.FE.bed >banding_sham.H3K27me3.intersect.bed
bedtools intersect -wao -a ../H3K27me3_Heart_sham_macs14/H3K27me3_Heart_sham_peaks.FE.bed -b ../H3K27me3_Heart_TAC_macs14/H3K27me3_Heart_TAC_peaks.FE.bed >sham_banding.H3K27me3.intersect.bed
prefix=H3K27me3
awk '$7<0' banding_sham.${prefix}.intersect.bed |cut -f1-3 >banding_specific.${prefix}.bed;
awk '$7<0' sham_banding.${prefix}.intersect.bed |cut -f1-3 >sham_specific.${prefix}.bed;
awk '$7>0' banding_sham.${prefix}.intersect.bed |awk '$5/$10>1.5 && $3-$2>=20'|awk -v OFS="\t" -v ORS="" 'BEGIN{i=1}{print $1"\t";if($2<=$7){print $7"\t"}else{print $2"\t"}if($3<=$8){print $3"\t"}else{print $8"\t"}print "banding_up_"i"\n";i+=1;}' >banding_up.${prefix}.bed
awk '$7>0' banding_sham.${prefix}.intersect.bed |awk '$5/$10<0.5 && $3-$2>=20'|awk -v OFS="\t" -v ORS="" 'BEGIN{i=1}{print $1"\t";if($2<=$7){print $7"\t"}else{print $2"\t"}if($3<=$8){print $3"\t"}else{print $8"\t"}print "banding_down_"i"\n";i+=1;}' >banding_down.${prefix}.bed
awk '$7>0' banding_sham.${prefix}.intersect.bed |awk '$5/$10>=0.5 && $5/$10<=1.5 && $3-$2>=20'|awk -v OFS="\t" -v ORS="" 'BEGIN{i=1}{print $1"\t";if($2<=$7){print $7"\t"}else{print $2"\t"}if($3<=$8){print $3"\t"}else{print $8"\t"}print "unchange_"i"\n";i+=1;}' >banding_sham.unchange.${prefix}.bed

###1.1.2.3 others
####SUZ12
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
macs2 callpeak -t SUZ12-ChIP-heart-WT-SRR1297213.final.rmDup.bam -c Input-ChIP-heart-WT-SRR1297212.final.rmDup.bam -f BAM -g mm --outdir suz12_v0 -n suz12 >SUZ12-ChIP-heart-WT-SRR1297213.macs2.log 2>SUZ12-ChIP-heart-WT-SRR1297213.macs2.err
macs2 callpeak -t SUZ12-ChIP-heart-WT-SRR1297213.final.rmDup.bam -c Input-ChIP-heart-WT-SRR1297212.final.rmDup.bam -f BAM -g mm --pvalue 1e-5 --broad-cutoff 1e-3 --broad --outdir suz12_broad_v1 -n suz12 >SUZ12-ChIP-heart-WT-SRR1297213.macs2_broad.log 2>SUZ12-ChIP-heart-WT-SRR1297213.macs2_broad.err
macs2 callpeak -t SUZ12-ChIP-heart-WT-SRR1297213.final.rmDup.bam -c Input-ChIP-heart-WT-SRR1297212.final.rmDup.bam -f BAM -g mm --pvalue 1e-5 --broad-cutoff 0.05 --broad --outdir suz12_broad_v2 -n suz12 > SUZ12-ChIP-heart-WT-SRR1297213.macs2_broad_3.log 2> SUZ12-ChIP-heart-WT-SRR1297213.macs2_broad_3.err
macs2 callpeak -t SUZ12-ChIP-heart-WT-SRR1297213.final.rmDup.bam -c Input-ChIP-heart-WT-SRR1297212.final.rmDup.bam -f BAM -g mm --pvalue 1e-5 --broad-cutoff 0.01 --broad --outdir suz12_broad_v3 -n suz12 > SUZ12-ChIP-heart-WT-SRR1297213.macs2_broad_4.log 2> SUZ12-ChIP-heart-WT-SRR1297213.macs2_broad_4.err
####EED,HDAC1,HDAC2
ls *_1.fastq |sed 's/_1.fastq//'|while read file;do
    cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -o ${file}_1.cutAda.fq -p ${file}_2.cutAda.fq -O 10 ${file}_1.fastq ${file}_2.fastq >${file}.cutAda.log 2>${file}.cutAda.err 
done;

ls *_1.cutAda.fq |sed 's/_1.cutAda.fq//'|while read file;do 
    bwa mem -t 15 /mnt/share/liym/data/bwa/mm10/mm10.fa ${file}_1.cutAda.fq ${file}_2.cutAda.fq 2>${file}.bwa.log |samtools view -@ 5 -bSu - |samtools sort -@ 5 -o ${file}.out.sorted.bam - 
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
cd EED
macs2 callpeak -t EED_ChIP_rep1.SRR4241119.filtered.rmdup.bam -c EED_ChIP_input.SRR4241121.filtered.rmdup.bam -f BAMPE -g mm --outdir EED_rep1 -n EED_rep1 > EED_ChIP_rep1.SRR4241119.macs2.log 2> EED_ChIP_rep1.SRR4241119.macs2.err
macs2 callpeak -t EED_ChIP_rep2.SRR4241120.filtered.rmdup.bam -c EED_ChIP_input.SRR4241121.filtered.rmdup.bam -f BAMPE -g mm --outdir EED_rep2 -n EED_rep2 > EED_ChIP_rep2.SRR4241120.macs2.log 2> EED_ChIP_rep2.SRR4241120.macs2.err
macs2 callpeak -t EED_ChIP_rep1.SRR4241119.filtered.rmdup.bam -c EED_ChIP_input.SRR4241121.filtered.rmdup.bam -f BAMPE -g mm --pvalue 1e-5 --broad-cutoff 1e-3 --broad --outdir EED_rep1_broad_2 -n EED_rep1 > EED_ChIP_rep1.SRR4241119.macs2_broad_2.log 2> EED_ChIP_rep1.SRR4241119.macs2_broad_2.err
macs2 callpeak -t EED_ChIP_rep2.SRR4241120.filtered.rmdup.bam -c EED_ChIP_input.SRR4241121.filtered.rmdup.bam -f BAMPE -g mm --pvalue 1e-5 --broad-cutoff 1e-3 --broad --outdir EED_rep2_broad_2 -n EED_rep2 > EED_ChIP_rep2.SRR4241120.macs2_broad_2.log 2> EED_ChIP_rep2.SRR4241120.macs2_broad_2.err
cd ../HDAC1
macs2 callpeak -t HDAC1_ChIP_EEDhete.lib1643_HVLCNCCXX_L4.filtered.rmdup.bam -c ../HDAC2/HDAC2_input_EEDhete.SRR4241137.filtered.rmdup.bam -f BAMPE -g mm --outdir EEDhete -n HDAC1 >HDAC1_ChIP_EEDhete.lib1643.macs2.log 2>HDAC1_ChIP_EEDhete.lib1643.macs2.err
macs2 callpeak -t HDAC1_ChIP_EEDhomo.lib1641_HVLCNCCXX_L4.filtered.rmdup.bam -c ../HDAC2/HDAC2_input_EEDhete.SRR4241137.filtered.rmdup.bam -f BAMPE -g mm --outdir EEDhomo -n HDAC1 >HDAC1_ChIP_EEDhomo.lib1641.macs2.log 2>HDAC1_ChIP_EEDhomo.lib1641.macs2.err
macs2 callpeak -t HDAC1_ChIP_EEDhete.lib1643_HVLCNCCXX_L4.filtered.rmdup.bam -c ../HDAC2/HDAC2_input_EEDhete.SRR4241137.filtered.rmdup.bam -f BAMPE -g mm --pvalue 1e-5 --broad-cutoff 1e-3 --broad --outdir EEDhete_broad_2 -n HDAC1 >HDAC1_ChIP_EEDhete.lib1643.macs2_broad_2.log 2>HDAC1_ChIP_EEDhete.lib1643.macs2_broad_2.err
macs2 callpeak -t HDAC1_ChIP_EEDhomo.lib1641_HVLCNCCXX_L4.filtered.rmdup.bam -c ../HDAC2/HDAC2_input_EEDhete.SRR4241137.filtered.rmdup.bam -f BAMPE -g mm --pvalue 1e-5 --broad-cutoff 1e-3 --broad --outdir EEDhomo_broad_2 -n HDAC1 >HDAC1_ChIP_EEDhomo.lib1641.macs2_broad_2.log 2>HDAC1_ChIP_EEDhomo.lib1641.macs2_broad_2.err
cd ../HDAC2
macs2 callpeak -t HDAC2_ChIP_EEDhete.SRR4241136.filtered.rmdup.bam -c HDAC2_input_EEDhete.SRR4241137.filtered.rmdup.bam -f BAMPE -g mm --outdir EEDhete -n HDAC2 >HDAC2_ChIP_EEDhete.SRR4241136.macs2.log 2>HDAC2_ChIP_EEDhete.SRR4241136.macs2.err
macs2 callpeak -t HDAC2_ChIP_EEDhomo.lib1642.filtered.rmdup.bam -c HDAC2_input_EEDhete.SRR4241137.filtered.rmdup.bam -f BAMPE -g mm --outdir EEDhomo -n HDAC2 >HDAC2_ChIP_EEDhomo.lib1642.macs2.log 2>HDAC2_ChIP_EEDhomo.lib1642.macs2.err
macs2 callpeak -t HDAC2_ChIP_EEDhete.SRR4241136.filtered.rmdup.bam -c HDAC2_input_EEDhete.SRR4241137.filtered.rmdup.bam -f BAMPE -g mm --pvalue 1e-5 --broad-cutoff 1e-3 --broad --outdir EEDhete_broad_2 -n HDAC2 >HDAC2_ChIP_EEDhete.SRR4241136.macs2_broad_2.log 2>HDAC2_ChIP_EEDhete.SRR4241136.macs2_broad_2.err
macs2 callpeak -t HDAC2_ChIP_EEDhomo.lib1642.filtered.rmdup.bam -c HDAC2_input_EEDhete.SRR4241137.filtered.rmdup.bam -f BAMPE -g mm --pvalue 1e-5 --broad-cutoff 1e-3 --broad --outdir EEDhomo_broad_2 -n HDAC2 >HDAC2_ChIP_EEDhomo.lib1642.macs2_broad_2.log 2>HDAC2_ChIP_EEDhomo.lib1642.macs2_broad_2.err

##1.2 RNA-seq data
cd /rd1/user/liym/nucleosomeTurnover/RNA-seq
###1.2.1 Pre-mapping filter
cd EEDheto_2
fqAdapterFilter.pl -a GAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCT,AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATG -i read1_filter.fq -o read1_adaF_filter.fq --a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCC --i2 read2_filter.fq --o2 read2_adaF_filter.fq -r fqAdaFilter.log &
cd ../EEDko_1
fqAdapterFilter.pl -a GAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCT,AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATG -i read1_filter.fq -o read1_adaF_filter.fq --i2 read2_filter.fq --o2 read2_adaF_filter.fq -r fqAdaFilter.log &
cd ../EEDko_2
fqAdapterFilter.pl -a GAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCT,AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATG,GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGC -i read1_filter.fq -o read1_adaF_filter.fq --a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCC,GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCG,GAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCA --i2 read2_filter.fq --o2 read2_adaF_filter.fq -r fqAdaFilter.log &
###1.2.2 Mapping and expression calculation
ls -d *_*/|sed 's;/;;'|while read dir;do
    cd $dir;
    /rd1/user/liym/tools/tophat-2.1.1.Linux_x86_64/tophat2 --read-mismatches 4 --read-gap-length 3 --read-edit-dist 4 --min-anchor 4 --splice-mismatches 0 --num-threads 10 --no-coverage-search --segment-length 25 --segment-mismatches 2 --library-type fr-unstranded /mnt/share/liym/data/bowtie2/mm10/mm10 read1_adaF_filter.fq read2_adaF_filter.fq > tophat2.log 2> tophat2.err;
    samtools view -bu -q 50 -@ 5 tophat_out/accepted_hits.bam | samtools sort -@ 5 -m 10G -o ${dir}.uniq.sorted.bam - 2>samtoolsSort.log
    cd ../;
    done;
ls -d */|grep -v "data"|sed 's;/;;'|while read dir;do cd $dir;cufflinks -p 5 -o cufflinks -G /rd1/user/liym/nucleosomeTurnover/data/mm10.refGene.filterRandom.gtf --max-bundle-frags 10000000 --no-update-check ${dir}.uniq.sorted.bam >cufflinks.log 2>cufflinks.err;cd ../;done;

##1.3 MNase-seq data
###1.3.1 Pre-mapping filter
cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -o 1_MNase-EEDhetoMF-lib1612_Rep2_1.cutAda.fq -p 1_MNase-EEDhetoMF-lib1612_Rep2_2.cutAda.fq -O 10 -u 1 -U 1 1_MNase-EEDhetoMF-lib1612_Rep2_1.fastq 1_MNase-EEDhetoMF-lib1612_Rep2_2.fastq > 1_MNase-EEDhetoMF-lib1612.cutAda.log 2> 1_MNase-EEDhetoMF-lib1612.cutAda.err
cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -o 1_MNase-EEDkoMF-lib1613_Rep2_1.cutAda.fq -p 1_MNase-EEDkoMF-lib1613_Rep2_2.cutAda.fq -O 10 -u 1 -U 1 1_MNase-EEDkoMF-lib1613_Rep2_1.fastq 1_MNase-EEDkoMF-lib1613_Rep2_2.fastq > 1_MNase-EEDkoMF-lib1613.cutAda.log 2> 1_MNase-EEDkoMF-lib1613.cutAda.err
cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -o 1_MNase-EEDhetoMF-lib1610_Rep1_1.cutAda.fq -p 1_MNase-EEDhetoMF-lib1610_Rep1_2.cutAda.fq -O 10 -u 1 -U 1 1_MNase-EEDhetoMF-lib1610_Rep1_1.fastq 1_MNase-EEDhetoMF-lib1610_Rep1_2.fastq >1_MNase-EEDhetoMF-lib1610.cutAda.log 2>1_MNase-EEDhetoMF-lib1610.cutAda.err
cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -o 1_MNase-EEDkoMF-lib1611_Rep1_1.cutAda.fq -p 1_MNase-EEDkoMF-lib1611_Rep1_2.cutAda.fq -O 10 -u 1 -U 1 1_MNase-EEDkoMF-lib1611_Rep1_1.fastq 1_MNase-EEDkoMF-lib1611_Rep1_2.fastq >1_MNase-EEDkoMF-lib1611.cutAda.log 2>1_MNase-EEDkoMF-lib1611.cutAda.err
###1.3.2 Mapping
bwa mem -t 8 /mnt/share/liym/data/bwa/mm10/mm10.fa 1_MNase-EEDhetoMF-lib1612_Rep2_1.cutAda.fq 1_MNase-EEDhetoMF-lib1612_Rep2_2.cutAda.fq 2> 1_MNase-EEDhetoMF-lib1612.bwa.log | samtools view -bSu - | samtools sort - 1_MNase-EEDhetoMF-lib1612.out.sorted
bwa mem -t 8 /mnt/share/liym/data/bwa/mm10/mm10.fa 1_MNase-EEDkoMF-lib1613_Rep2_1.cutAda.fq 1_MNase-EEDkoMF-lib1613_Rep2_2.cutAda.fq 2> 1_MNase-EEDkoMF-lib1613.log | samtools view -bSu - | samtools sort - 1_MNase-EEDkoMF-lib1613.out.sorted
bwa mem -t 8 /mnt/share/liym/data/bwa/mm10/mm10.fa 1_MNase-EEDhetoMF-lib1610_Rep1_1.cutAda.fq 1_MNase-EEDhetoMF-lib1610_Rep1_2.cutAda.fq 2> 1_MNase-EEDhetoMF-lib1610.bwa.log |samtools view -bSu - | samtools sort - 1_MNase-EEDhetoMF-lib1610.out.sorted 
bwa mem -t 8 /mnt/share/liym/data/bwa/mm10/mm10.fa 1_MNase-EEDkoMF-lib1611_Rep1_1.cutAda.fq 1_MNase-EEDkoMF-lib1611_Rep1_2.cutAda.fq 2> 1_MNase-EEDkoMF-lib1611.bwa.log |samtools view -bSu - | samtools sort - 1_MNase-EEDkoMF-lib1611.out.sorted
###1.3.3 Post-mapping filter
ls *.out.sorted.bam|while read file;do
    prefix=$(echo $file|sed 's/.out.sorted.bam//')
    bamtools filter -in $file -out ${prefix}.myFiltered.sorted.bam -script /mnt/share/liym/data/bwa/bwaMem.filter.json 2>log/bamtoolsFilter.log;
    bamtools stats -in ${prefix}.myFiltered.sorted.bam >${prefix}.myFiltered.sorted.bamStats & 
    samtools index ${prefix}.myFiltered.sorted.bam;
    samtools view -@ 5 -bu ${prefix}.myFiltered.sorted.bam chr1 chr2 chrX chr3 chr4 chr5 chr6 chr7 chr10 chr8 chr14 chr9 chr11 chr13 chr12 chr15 chr16 chr17 chrY chr18 chr19 chrM |samtools sort -@ 5 -o ${prefix}.final - 2>log/samtoolsSort.log;
    samtools rmdup ${prefix}.final.bam ${prefix}.final.rmDup.bam 2>${prefix}.rmDup.log 
    bamtools stats -in ${prefix}.final.bam >${prefix}.final.bamStats & 
done;
ls *.final.rmDup.bam|sed 's/.bam//'|while read file;do bamtools filter -in ${file}.bam -out ${file}.SizeF.bam -script bwaMem.insertSize.json ; bamtools stats -insert -in ${file}.SizeF.bam >${file}.SizeF.bamStats & done;
###1.3.4 Generate nucleosome occupancy
danpos.py dpos 1_MNase-EEDhetoMF-lib1610.final.rmDup.SizeF.bam,1_MNase-EEDhetoMF-lib1612.final.rmDup.SizeF.bam,1_MNase-EEDkoMF-lib1611.final.rmDup.SizeF.bam,1_MNase-EEDkoMF-lib1613.final.rmDup.SizeF.bam -jd 147 --extend 74 -m 1 -c 10000000 -o danpos_rst > danpos.log 2> danpos.err
ls danpos_rst/pooled/*.wig|while read file;do prefix=$(echo $file|sed 's/.final.rmDup.SizeF.Fnor.smooth.wig//');wigToBigWig -clip $file /rd1/user/liym/nucleosomeTurnover/data/mm10.filterRandom.size ${prefix}.bw 2>wigTobw.log ;done;

#2. Final figures
cd /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/figures
##2.1 Figure1
###Generate H2BGFP signals for 6 time points
toolRunner.sh wigmath.Average -f --step 10 -p 2 -o week0.wig /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/1/danposRst/result/pooled/1_H2bGFP-day1_Rep1.final.bgsub.Fnor.smooth.wig /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/1/danposRst/result/pooled/1_H2bGFP-day1_Rep2.final.bgsub.Fnor.smooth.wig > week0.wigmath.log
toolRunner.sh wigmath.Average -f --step 10 -p 2 -o week1.wig /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/1/danposRst/result/pooled/1_H2bGFP-week1_Rep1.final.bgsub.Fnor.smooth.wig /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/1/danposRst/result/pooled/1_H2bGFP-week1_Rep2.final.bgsub.Fnor.smooth.wig > week1.wigmath.log
toolRunner.sh wigmath.Average -f --step 10 -p 2 -o week2.wig /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/1/danposRst/result/pooled/1_H2bGFP-week2_Rep1.final.bgsub.Fnor.smooth.wig /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/1/danposRst/result/pooled/1_H2bGFP-week2_Rep2.final.bgsub.Fnor.smooth.wig > week2.wigmath.log
toolRunner.sh wigmath.Average -f --step 10 -p 2 -o week4.wig /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/1/danposRst/result/pooled/1_H2bGFP-week4_Rep1.final.bgsub.Fnor.smooth.wig /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/2/danposRst/result/pooled/2_H2BGFP-CM-week4-2015_Rep1.final.bgsub.Fnor.smooth.wig /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/2/danposRst/result/pooled/2_H2BGFP-CM-week4-2015_Rep2.final.bgsub.Fnor.smooth.wig > week4.wigmath.log
ls *.wig |sed 's/.wig//'|while read file;do wigToBigWig -clip ${file}.wig /mnt/share/liym/data/chr.size/mm10.filterRandom.size ${file}.bw 2>log/${file}.wigTobw.log;done;
toolRunner.sh wigmath.Average -f --step 10 -p 2 -o week0-new.wig /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/1/danposRst/day1_input2/pooled/1_H2bGFP-day1_Rep1.final.bgsub.Fnor.smooth.wig /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/1/danposRst/day1_input2/pooled/1_H2bGFP-day1_Rep2.final.bgsub.Fnor.smooth.wig > log/week0-new.wigmath.log
wigToBigWig -clip week0-new.wig /rd1/user/liym/nucleosomeTurnover/data/mm10.filterRandom.size week0-new.bw
###H2BGFP signals for 1kb genomic intervals
ls *.bw |sed 's/.bw//'|while read file;do bwtool summary /rd1/user/liym/nucleosomeTurnover/data/mm10.1kbIntervals.filterGap.bed ${file}.bw ${file}.1kb.tsv;done;
paste <(cut -f8 H3.1kb.tsv) <(cut -f8 week0-new.1kb.tsv) <(cut -f8 week1.1kb.tsv) <(cut -f8 week2.1kb.tsv) <(cut -f8 week4.tsv) <(cut -f8 week6.1kb.tsv) <(cut -f8 week8.1kb.tsv) >H3.H2BGFP.1kb.sum.tsv
####Start in R
R
library("ggplot2")
library("reshape2")
data<-read.delim(file="H3.H2BGFP.1kb.sum.tsv",header=F)
colnames(data)<-c("H3","0w","1w","2w","4w","6w","8w")
quantile(data$H3,probs=seq(0,1,0.1),na.rm=T)
#  0%  10%  20%  30%  40%  50%  60%  70%  80%  90% 100% 
#  0.00 0.00 0.01 0.04 0.10 0.17 0.26 0.36 0.51 0.73 6.58 
quantile(data$'0w',probs=seq(0,1,0.1),na.rm=T)
#    0%   10%   20%   30%   40%   50%   60%   70%   80%   90%  100% 
#     0.00  0.00  0.02  0.08  0.14  0.22  0.30  0.39  0.50  0.66 43.39 
subset<-data[data$H3>0.10 & data$'0w'>0.14,]
subset<-subset[,-1]
subset<-melt(subset)
pdf(file="Fig1E.boxplot.final.pdf")
ggplot(subset,aes(x=variable, y=value, fill=variable))+geom_violin(trim=F)+geom_boxplot(width=0.2, fill="white")+scale_fill_brewer(palette="Blues",direction = -1) + theme_classic()+ylim(0,1.2)+labs(y="H2BGFP signals",x="Time")+stat_summary(fun.y=mean, geom="point", shape=23, size=2)
dev.off()
wilcox.test(subset[subset$variable=="0w",2],subset[subset$variable=="1w",2],alternative="g")#p-value <2.2e-16
wilcox.test(subset[subset$variable=="1w",2],subset[subset$variable=="2w",2],alternative="g")#p-value <2.2e-16
wilcox.test(subset[subset$variable=="2w",2],subset[subset$variable=="4w",2],alternative="g")#p-value < 2.2e-16
wilcox.test(subset[subset$variable=="4w",2],subset[subset$variable=="6w",2],alternative="g")#p-value < 2.2e-16
wilcox.test(subset[subset$variable=="6w",2],subset[subset$variable=="8w",2],alternative="g")#p-value < 2.2e-16
####End in R
###Gene profiles and heatmap
bwtool agg 5000:10000:5000 /rd1/user/liym/nucleosomeTurnover/data/mm10.refGene.bed6 H3.bw,week0.bw,week1.bw,week2.bw,week4.bw,week6.bw,week8.bw gene.profile.txt -header 
Rscript /mnt/share/liym/bin/smooth_lines.R -d=1000 -i=gene.profile.txt -x="Gene relative position" -y="Normalized H3/H2BGFP depth" -c="#FF0000FF,#FFDB00FF,#49FF00FF,#00FF92FF,#0092FFFF,#4900FFFF,#FF00DBFF" -l=H3,H2BGFP_0w,H2BGFP_1w,H2BGFP_2w,H2BGFP_4w,H2BGFP_6w,H2BGFP_8w -o=gene.profile.pdf
bwtool agg 5000:10000:5000 /rd1/user/liym/nucleosomeTurnover/data/mm10.refGene.bed6 week0-new.bw gene.profile.week0-new.txt
paste <(cut -f1,2,4- gene.profile.txt) <(cut -f2 gene.profile.week0-new.txt) |awk -v OFS="\t" '{print $1,$2,$8,$3,$4,$5,$6,$7}' >gene.profile.final.txt
Rscript /mnt/share/liym/bin/smooth_lines.R -d=1000 -i=gene.profile.final.txt -x="Gene relative position" -y="Normalized H3/H2BGFP depth" -c="#FF0000FF,#FFDB00FF,#49FF00FF,#00FF92FF,#0092FFFF,#4900FFFF,#FF00DBFF" -l=H3,H2BGFP_0w,H2BGFP_1w,H2BGFP_2w,H2BGFP_4w,H2BGFP_6w,H2BGFP_8w -o=gene.profile.final.pdf
mkdir geneReion && cd geneReion
/mnt/share/liym/tools/deepTools-2.0.0/bin/computeMatrix scale-regions -R ../../fig2/EEDhete.Meanfpkm.bed6 -S ../H3.bw ../week0-new.bw ../week1.bw ../week2.bw ../week4.bw ../week6.bw ../week8.bw -out geneRegions.matrix.gz --outFileNameMatrix geneRegions.matrix.tsv --regionBodyLength 10000 --startLabel TSS --endLabel TTS -b 5000 -a 5000 -bs 40 --skipZeros -p 15 >log/computeMatrix.log 2>log/computeMatrix.err
gunzip geneRegions.matrix.gz
####Start in R
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
####End in R
#####Row mean sort
paste (grep -v "@" geneRegions.matrix|cut -f1-6) final.matrix.v2 >tmp;
cat header tmp >final.matrix.v2
gzip -k final.matrix.v2
/mnt/share/liym/tools/deepTools-2.0.0/bin/plotHeatmap -m final.matrix.v2.gz -out final.matrix.rowMean.sort.heatmap.pdf --sortRegions descend --sortUsing mean --colorList royalblue,white,red --colorNumber 5 --whatToShow "heatmap and colorbar" --startLabel TSS --endLabel TTS -x "" -y "" --samplesLabel H3 0w 1w 2w 4w 6w 8w --heatmapHeight 16 -max 1
#####H3 TSS signals sort
grep -v "^@" final.matrix.v2|Rscript /mnt/share/liym/bin/sortByMultiColumns.R -s=100 -e=150 -d=FALSE -o=final.matrix.H3.TSS.sort
cat header final.matrix.H3.TSS.sort >tmp;mv tmp final.matrix.H3.TSS.sort
gzip final.matrix.H3.TSS.sort
/mnt/share/liym/tools/deepTools-2.0.0/bin/plotHeatmap -m final.matrix.H3.TSS.sort.gz -out final.matrix.H3.TSS.sort.heatmap.pdf --sortRegions no --colorList royalblue,white,red --colorNumber 5 --whatToShow "heatmap and colorbar" --startLabel TSS --endLabel TTS -x "" -y "" --samplesLabel H3 0w 1w 2w 4w 6w 8w --heatmapHeight 16 -max 1 
#####H3 gene region signals sort
grep -v "^@" final.matrix.v2|Rscript /mnt/share/liym/bin/sortByMultiColumns.R -s=125 -e=375 -d=FALSE -o=final.matrix.H3.geneRegion.sort
cat header final.matrix.H3.geneRegion.sort >tmp;mv tmp final.matrix.H3.geneRegion.sort
gzip final.matrix.H3.geneRegion.sort
/mnt/share/liym/tools/deepTools-2.0.0/bin/plotHeatmap -m final.matrix.H3.geneRegion.sort.gz -out final.matrix.H3.geneRegion.sort.heatmap.pdf --sortRegions no --colorList royalblue,white,red --colorNumber 5 --whatToShow "heatmap and colorbar" --startLabel TSS --endLabel TTS -x "" -y "" --samplesLabel H3 0w 1w 2w 4w 6w 8w --heatmapHeight 16 -max 1 
#####Gene FPKM sort
sort -k5,5n final.matrix.v2 >final.fpkmSort.matrix
cat header final.fpkmSort.matrix >tmp;mv tmp final.fpkmSort.matrix
#####Pol2 signals sort
cut -f1-4 ../../fig2/EEDhete.Meanfpkm.TSSfl500.bed6 | bwtool summary -skip-median stdin /rd1/user/liym/nucleosomeTurnover/ChIP-seq/TF/GSE52123-2014NC/Pol2.subtract.bw TSSfl500.pol2.signals.tsv
join.pl -i2 final.matrix.v2 -f2 4 -i1 TSSfl500.pol2.signals.bed6+ -f1 4|cut -f7-|sort -k1,1n|cut -f2- >final.pol2Sort.matrix.v2
cat header final.pol2Sort.matrix >tmp;mv tmp final.pol2Sort.matrix.v2
gzip final.pol2Sort.matrix.v2
/mnt/share/liym/tools/deepTools-2.0.0/bin/plotHeatmap -m final.matrix.v2.pol2Sort.gz  -out final.heattmap.pol2Sort.pdf --sortRegions no --colorList royalblue,white,red --colorNumber 5 --whatToShow "heatmap and colorbar" --startLabel TSS --endLabel TTS -x "" -y "" --samplesLabel H3 0w 1w 2w 4w 6w 8w --heatmapHeight 16 -max 1 

##2.2 Figure2
###2.2.1 1kb genomic intervals
cut -f2- ../fig1/H3.H2BGFP.1kb.sum.tsv >H2BGFP.1kb.sum.tsv 
grep -v "NA" H2BGFP.1kb.sum.tsv>tmp;mv tmp H2BGFP.1kb.sum.tsv
mkdir split plus1e-4 plus1e-5 plus1e-15
cd split;split -l 10000 -d ../H2BGFP.1kb.sum.tsv
ls split/*|while read file;do prefix=$(basename $file);Rscript /rd1/user/liym/nucleosomeTurnover/scripts/NTR.R -i=$file -s=1 -e=6 -p=0.0001 -o=plus1e-4/${prefix}.NTR.tsv ;done;
#ls split/*|while read file;do prefix=$(basename $file);Rscript /rd1/user/liym/nucleosomeTurnover/scripts/NTR.R -i=$file -s=1 -e=6 -p=0.00001 -o=plus1e-5/${prefix}.NTR.tsv ;done;
#ls split/*|while read file;do prefix=$(basename $file);Rscript /rd1/user/liym/nucleosomeTurnover/scripts/NTR.R -i=$file -s=1 -e=6 -p=1e-15 -o=plus1e-15/${prefix}.NTR.tsv ;done;
cd plus1e-4
paste ../intervals_1k.bed3 NTR.summary.tsv >regions.NTR.p.tsv
sort -k10 -g -r regions.NTR.p.tsv >regions.NTR.p.sorted.tsv
awk '$11<0.05' regions.NTR.p.sorted.tsv |head -n10000 |cut -f1-3 >top10k.p0.05.bed3
Rscript /mnt/share/liym/bin/PeakAnnoGO.R top10k.p0.05.bed3 top10k.p0.05.Anno.tsv top10k.p0.05.annoPie.pdf top10k.p0.05.GO.dotplot.pdf
#####Start in R
ggplot(subset20,aes(x=p.adjust,y=Description,color=GeneRatio,size=Count))+geom_point(stat='identity')+scale_colour_gradientn(colours = c("blue","yellow","red"))+xlim(9,11.5)+geom_text(aes(label = Description), size=3, colour="black")+theme(axis.text.y = element_blank(),panel.background = element_rect(fill="white"),panel.grid.major = element_line(colour="gray"), panel.grid.minor = element_line(colour="gray"),panel.border = element_rect(colour="black",fill=NA))
#####End in R
###2.2.2 NTR binnd gene expression
awk -v OFS="\t" '{if($1!="tracking_id"){print $1":"$7,$10}}' /rd1/user/liym/nucleosomeTurnover/RNA-seq/EEDheto_1/cufflinks/genes.fpkm_tracking >EEDheto_1.fpkm.tsv
awk -v OFS="\t" '{if($1!="tracking_id"){print $1":"$7,$10}}' /rd1/user/liym/nucleosomeTurnover/RNA-seq/EEDheto_2/cufflinks/genes.fpkm_tracking >EEDheto_2.fpkm.tsv
cut -f2,3,4,12 /rd1/user/liym/nucleosomeTurnover/data/mm10.refGene.filterRandom.gpe|sort|uniq|awk -v OFS="\t" '{print $4":"$1":"$3,$2}' >mm10.refGene.start.tsv 
join.pl -i1 EEDheto_1.fpkm.tsv -i2 EEDheto_2.fpkm.tsv |awk -v OFS="\t" '{split($1,a,":");split(a[3],b,"-");print a[2],b[1],b[2],a[1]":"a[2]":"b[1],($2+$4)/2}'|join.pl -f1 4 -i2 mm10.refGene.start.tsv -f2 1|awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$7}'|sort|uniq >EEDhete.Meanfpkm.bed6 
awk -v OFS="\t" '{if($6=="+"){print $1,$2-500,$2+500,$4,$5,$6}else{print $1,$3-500,$3+500,$4,$5,$6}}' EEDhete.Meanfpkm.bed6 >EEDhete.Meanfpkm.TSSfl500.bed6 
ls ../fig1/week*.bw |grep -v "week0.bw"|while read file;do prefix=$(basename $file|sed 's/.bw//'|sed 's/-new//');cut -f1-4 EEDhete.Meanfpkm.TSSfl500.bed6 |bigWigAverageOverBed $file stdin ${prefix}.TSSfl500.H2BGFP.tsv;done;
paste <(cut -f1,5 week0.TSSfl500.H2BGFP.tsv) <(cut -f5 week1.TSSfl500.H2BGFP.tsv) <(cut -f5 week2.TSSfl500.H2BGFP.tsv) <(cut -f5 week4.TSSfl500.H2BGFP.tsv) <(cut -f5 week6.TSSfl500.H2BGFP.tsv) <(cut -f5 week8.TSSfl500.H2BGFP.tsv)|join.pl -f1 1 -i2 EEDhete.Meanfpkm.TSSfl500.bed6 -f2 4|cut -f1-7,12|sort -k8 -gr >TSS.fl500.H2BGFP.fpkmSort.tsv
Rscript /rd1/user/liym/nucleosomeTurnover/scripts/NTR.R -i=TSS.fl500.H2BGFP.fpkmSort.tsv -p=1e-4 -s=2 -e=7 -o=TSS.fl500.H2BGFP.fpkmSort.NTR.tsv 
###2.2.3 NTR of mutiple histone modification regions
#####Agg lines
ls /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/*.bed |grep -vE "H3K27ac|H3K27me3"|while read file;do prefix=$(basename $file|sed 's/.ChIP.8w.mouseHeart.narrowPeak.bed//');computeMatrix reference-point -R $file -S ../fig1/week0-new.bw ../fig1/week1.bw ../fig1/week2.bw ../fig1/week4.bw ../fig1/week6.bw ../fig1/week8.bw -out ${prefix}.H2BGFP.gz --referencePoint center -b 5000 -a 5000 -bs 20 --sortRegions no -p 20 >computeMatrix.log 2>computeMatrix.err;done;
computeMatrix reference-point -R /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/result/H3K27ac-WT.intersect.broadPeaks.bed -S ../fig1/week0-new.bw ../fig1/week1.bw ../fig1/week2.bw ../fig1/week4.bw ../fig1/week6.bw ../fig1/week8.bw -out H3K27ac.H2BGFP.gz --referencePoint center -b 5000 -a 5000 -bs 20 --sortRegions no -p 20 >computeMatrix.log 2>computeMatrix.err
computeMatrix reference-point -R /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/result/H3K27me3-WT.intersect.broadPeaks.bed -S ../fig1/week0-new.bw ../fig1/week1.bw ../fig1/week2.bw ../fig1/week4.bw ../fig1/week6.bw ../fig1/week8.bw -out H3K27me3.H2BGFP.gz --referencePoint center -b 5000 -a 5000 -bs 20 --sortRegions no -p 20 >computeMatrix.log 2>computeMatrix.err
ls *.H2BGFP.gz|sed 's/.H2BGFP.gz//'|while read file;do 
        mkdir $file;
        mv ${file}.H2BGFP.gz $file;
        cd $file;
        gunzip ${file}.H2BGFP.gz;
        /rd1/user/liym/nucleosomeTurnover/scripts/NTR.for.lines.sh ${file}.H2BGFP $file 
        Rscript /mnt/share/liym/bin/columnMean.R -i=${file}.summary.NTR.tsv -s=7 -e=506 -o=colMean.NTR.tsv
        cd ../;
done;
paste H3*/colMean.NTR.tsv >HM.lines.NTR.tsv
#####Violion plot
ls /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/*.bed|grep -vE "H3K27ac|H3K27me3"|while read file;do 
    prefix=$(basename $file|sed 's/.ChIP.8w.mouseHeart.narrowPeak.bed//');
    #/rd1/user/liym/nucleosomeTurnover/scripts/FileForNTR.sh $file ../fig1/week0-new.bw ../fig1/week1.bw ../fig1/week2.bw ../fig1/week4.bw ../fig1/week6.bw ../fig1/week8.bw ${prefix}.H2BGFP.tsv;
    Rscript /rd1/user/liym/nucleosomeTurnover/scripts/NTR.R -i=${prefix}.H2BGFP.tsv -p=0.0001 -s=4 -e=9 -o=${prefix}.NTR.tsv 
    done;
/rd1/user/liym/nucleosomeTurnover/scripts/FileForNTR.sh /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/result/H3K27me3-WT.intersect.broadPeaks.bed ../fig1/week0-new.bw ../fig1/week1.bw ../fig1/week2.bw ../fig1/week4.bw ../fig1/week6.bw ../fig1/week8.bw H3K27me3.peaks.H2BGFP.tsv 
Rscript /rd1/user/liym/nucleosomeTurnover/scripts/NTR.R -i=H3K27me3.peaks.H2BGFP.tsv -p=0.0001 -s=4 -e=9 -o=H3K27me3.peaks.NTR.tsv 
###2.2.4 NTR of eRNAs
awk -v OFS="\t" '{if($1!="tracking_id"){print $7,$10}}' /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/result/cufflinks_rep2/genes.fpkm_tracking >enhancers.fpkm.rep2.tsv 
awk -v OFS="\t" '{if($1!="tracking_id"){print $7,$10}}' /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/result/cufflinks_rep1/genes.fpkm_tracking >enhancers.fpkm.rep1.tsv 
join.pl -i1 enhancers.fpkm.rep1.tsv -f1 1 -i2 enhancers.fpkm.rep2.tsv -f2 1 |awk -v OFS="\t" '{split($1,a,":");split(a[2],b,"-");print a[1],b[1],b[2],($2+$4)/2}' >enhancers.Mean.fpkm.tsv 
bedtools intersect -f 1.0 -r -wao -a enhancers.Mean.fpkm.tsv -b enhancers.H2BGFP.NTR.tsv |cut -f1-4,14 |sort -k4,4n >enhancers.fpkm.NTR.tsv
bedtools subtract -a enhancers.fpkm.NTR.tsv -b /rd1/user/liym/nucleosomeTurnover/data/refGene.TSS.fl2k.bed6|bedtools subtract -a stdin -b /rd1/user/liym/nucleosomeTurnover/data/refGene.exons.bed6|awk '$3-$2>100' |cut -f1-3 >eRNA.regions.bed3
perl /mnt/share/zhangsj/bin/regionRPKM_speed.pl -b eRNA.regions.bed3 /rd1/user/liym/nucleosomeTurnover/RNA-seq/EEDheto_2/EEDheto_2.uniq.sorted.bam >eRNA.regions.EEDhete_rep2.RPKM.bed 2>eRNA.regions.EEDhete_rep2.RPKM.log 
perl /mnt/share/zhangsj/bin/regionRPKM_speed.pl -b eRNA.regions.bed3 /rd1/user/liym/nucleosomeTurnover/RNA-seq/EEDheto_1/EEDheto_1.uniq.sorted.bam > eRNA.regions.EEDhete_rep1.RPKM.bed 2> eRNA.regions.EEDhete_rep1.RPKM.log
/rd1/user/liym/nucleosomeTurnover/scripts/FileForNTR.sh eRNA.regions.bed3 ../fig1/week0-new.bw ../fig1/week1.bw ../fig1/week2.bw ../fig1/week4.bw ../fig1/week6.bw ../fig1/week8.bw eRNA.regions.H2BGFP.tsv
Rscript /rd1/user/liym/nucleosomeTurnover/scripts/NTR.R -i=eRNA.regions.H2BGFP.tsv -p=1e-4 -s=4 -e=9 -o=eRNA.regions.NTR.tsv
bedtools intersect -f 1.0 -r -wao -a eRNA.regions.EEDhete_rep1.RPKM.bed -b eRNA.regions.EEDhete_rep2.RPKM.bed|awk -v OFS="\t" '{print $1,$2,$3,($5+$10)/2}'|bedtools intersect -f 1.0 -r -wao -a stdin -b eRNA.regions.NTR.tsv |awk -v OFS="\t" '{print $1,$2,$3,$4,$14,$15}' >eRNA.regions.RPKM.NTR.tsv 
###2.2.5 Promoters and enhaners vs. pol2 binding
bedtools intersect -wao -a /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/result/H3K27ac-WT.intersect.broadPeaks.bed -b /rd1/user/liym/nucleosomeTurnover/ChIP-seq/TF/GSE52123-2014NC/result/pol2.intersect.merged.peaks.bed >enhancers.pol2.intersect.bed
awk '$5>0' enhancers.pol2.intersect.bed |cut -f1-3|sort|uniq|awk -v OFS="\t" '{print $1,$2,$3,"wt"}' >enhancers.pol2.status.bed6
awk '$5<0' enhancers.pol2.intersect.bed |cut -f1-3|sort|uniq|awk -v OFS="\t" '{print $1,$2,$3,"wo"}' >>enhancers.pol2.status.bed6
awk -v OFS="" -v ORS="" '{if($5>0){print $1"\t";if($2<$5){print $5"\t"}else{print $2"\t"}if($3<$6){print $3"\n"}else{print $6"\n"}}}' enhancers.pol2.intersect.bed >enhancers.wt.pol2.bed3 
/rd1/user/liym/nucleosomeTurnover/scripts/FileForNTR.sh enhancers.wt.pol2.bed3 ../fig1/week0-new.bw ../fig1/week1.bw ../fig1/week2.bw ../fig1/week4.bw ../fig1/week6.bw ../fig1/week8.bw enhancers.wt.pol2.H2BGFP.tsv 
Rscript /rd1/user/liym/nucleosomeTurnover/scripts/NTR.R -i=enhancers.wt.pol2.H2BGFP.tsv -p=0.0001 -s=4 -e=9 -o=enhancers.wt.pol2.NTR.tsv 
cat <(awk '$4=="wo"' enhancers.pol2.status.NTR.tsv) <(awk -v OFS="\t" '{print $1,$2,$3,"wt",$10}' enhancers.wt.pol2.NTR.tsv) >enhancers.pol2.status.NTR.final.tsv

##2.3 Figure3
ln -s ../fig1/week0-new.bw week0.bw
ln -s ../fig1/week1.bw week1.bw
ln -s ../fig2/week1.bw week2.bw
ln -s ../fig4/week1.bw week4.bw
ln -s ../fig6/week1.bw week6.bw
ln -s ../fig8/week1.bw week8.bw
ln -s /rd1/user/liym/nucleosomeTurnover/ChIP-seq/TF/GSE52123-2014NC/result/Gata4.intersect.merged.peaks.v1.bed Gata4.narrowPeak.bed
ln -s /rd1/user/liym/nucleosomeTurnover/ChIP-seq/TF/GSE35151_adultHeart/Nkx2-5/Nkx2-5_peaks.narrowPeak Nkx2-5.narrowPeak.bed
ln -s /rd1/user/liym/nucleosomeTurnover/ChIP-seq/TF/GSE29636_JCI2011_SE36/Tbx20_combined.peaks.bed3 Tbx20.narrowPeak.bed
ln -s /rd1/user/liym/nucleosomeTurnover/ChIP-seq/TF/GSE35151_adultHeart/Tbx3/Tbx3_peaks.narrowPeak Tbx3.narrowPeak.bed
ln -s /rd1/user/liym/nucleosomeTurnover/ChIP-seq/TF/GSE52123-2014NC/Gata4.subtract.bw Gata4.signal.bw
ln -s /rd1/user/liym/nucleosomeTurnover/ChIP-seq/TF/GSE35151_adultHeart/Nkx2-5.subtract.bw Nkx2-5.signal.bw
ln -s /rd1/user/liym/nucleosomeTurnover/ChIP-seq/TF/GSE29636_JCI2011_SE36/Tbx20.subtract.bw Tbx20.signal.bw
ln -s /rd1/user/liym/nucleosomeTurnover/ChIP-seq/TF/GSE35151_adultHeart/Tbx3.subtract.bw Tbx3.signal.bw
###2.3.1 enhancer NTR vs. TF binidng status
cat *narrowPeak.bed |cut -f1-3 >TFs.sum.bed3 
awk '$7=="wo"' ../fig2/promoters.pol2.status.NTR.tsv|bedtools intersect -wao -a stdin -b TFs.sum.bed3 >promoters.TFs.intersect.bed
awk '$10>0' promoters.TFs.intersect.bed |cut -f1-6,8|sort|uniq |awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,"wt",$7}' >promoters.TFs.status.NTR.tsv
awk '$10<0' promoters.TFs.intersect.bed |cut -f1-6,8|sort|uniq |awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,"wo",$7}' >>promoters.TFs.status.NTR.tsv
awk '$4=="wo"' ../fig2/enhancers.pol2.status.NTR.tsv|bedtools intersect -wao -a stdin -b TFs.sum.bed3 >enhancers.TFs.intersect.bed
awk '$7>0' enhancers.TFs.intersect.bed |awk -v OFS="" -v ORS="" '{print $1"\t";if($2<$7){print $7"\t";}else{print $2"\t"}if($3<$8){print $3"\n"}else{print $8"\n"}}'|awk -v OFS="\t" '{print $1,$2,$3,"wt"}' |sort -k1,1 -k2,2n |bedtools merge >enhancers.TFs.status.bed3+
/rd1/user/liym/nucleosomeTurnover/scripts/FileForNTR.sh enhancers.TFs.status.bed3+ week0.bw week1.bw week2.bw week4.bw week6.bw week8.bw enhancers.TFs.wt.H2BGFP.tsv 
Rscript /rd1/user/liym/nucleosomeTurnover/scripts/NTR.R -i=enhancers.TFs.wt.H2BGFP.tsv -s=4 -e=9 -p=0.0001 -o=enhancers.TFs.wt.NTR.tsv 
awk '$7<0' enhancers.TFs.intersect.bed |cut -f1-3,5|sort|uniq|awk -v OFS="\t" '{print $1,$2,$3,"wo",$4}' >enhancers.TFs.status.NTR.tsv
awk -v OFS="\t" '{print $1,$2,$3,"wt",$10}' enhancers.TFs.wt.NTR.tsv >>enhancers.TFs.status.NTR.tsv
###2.3.2 enhancers NTR binned by TF binding number
bedtools multiinter -i <(awk '$4=="wo"' ../fig2/enhancers.pol2.status.NTR.tsv|sort -k1,1 -k2,2n) <(sort -k1,1 -k2,2n Gata4.narrowPeak.bed) <(sort -k1,1 -k2,2n Nkx2-5.narrowPeak.bed) <(sort -k1,1 -k2,2n Tbx20.narrowPeak.bed) <(sort -k1,1 -k2,2n Tbx3.narrowPeak.bed)|awk -v OFS="\t" '{if($5~"1" && $5~","){print $1,$2,$3,$4}}' >enhancers.TFnum.bed3+
/rd1/user/liym/nucleosomeTurnover/scripts/FileForNTR.sh enhancers.TFnum.bed3+ week0.bw week1.bw week2.bw week4.bw week6.bw week8.bw enhancers.TFnum.H2BGFP.tsv
Rscript /rd1/user/liym/nucleosomeTurnover/scripts/NTR.R -i=enhancers.TFnum.H2BGFP.tsv -s=4 -e=9 -p=0.0001 -o=enhancers.TFnum.NTR.tsv
paste enhancers.TFnum.bed3+ <(cut -f10 enhancers.TFnum.NTR.tsv) >enhancers.TFnum.NTR.finla.tsv
###2.3.3 enhancers NTR binned by total TF binding signals
ls *.signal.bw |sed 's/.signal.bw//'|while read file;do bwtool summary enhancers.TFnum.bed3+ ${file}.signal.bw enhancers.TFnum.${file}.signal.tsv;done;
paste *TFnum*signal.tsv |cut -f1-3,8,17,16,25|bedtools intersect -wo -a stdin -b enhancers.TFnum.NTR.finla.tsv|awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$11,$12}' >enhancers.TFnum.allTFsignal.NTR.tsv

##2.4 Figure4 #Mutiple chromatin regulators
cd fig4
ln -s ../fig1/week0-new.bw week0.bw
ln -s ../fig1/week1.bw week1.bw
ln -s ../fig2/week1.bw week2.bw
ln -s ../fig4/week1.bw week4.bw
ln -s ../fig6/week1.bw week6.bw
ln -s ../fig8/week1.bw week8.bw
###2.4.1 enhancers NTR binned by signals
cut -f1-3,10 ../fig2/enhancers.H2BGFP.NTR.tsv >enhancers.NTR.bed3+ 
awk -v OFS="\t" 'BEGIN{i=1}{print $1,$2,$3,i;i++;}' enhancers.NTR.bed3+|bigWigAverageOverBed /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/1_H3K27ac-WT-AdultCM_Rep1.subtract.bw stdin enhancers.NTR.signal.rep1.tsv 
awk -v OFS="\t" 'BEGIN{i=1}{print $1,$2,$3,i;i++;}' enhancers.NTR.bed3+|bigWigAverageOverBed /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/1_H3K27ac-WT-AdultCM_Rep2.subtract.bw stdin enhancers.NTR.signal.rep2.tsv 
paste enhancers.NTR.bed3+ <(cut -f5 enhancers.NTR.signal.rep1.tsv) <(cut -f5 enhancers.NTR.signal.rep2.tsv)|awk -v OFS="\t" '{print $1,$2,$3,$4,($5+$6)/2}' |sort -k5,5n >enhancers.NTR.H3K27ac.signal.tsv 
###2.4.2 H2K27me3 intersect with EED/SUEZ12 peaks
bedtools intersect -wao -a /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/result/H3K27me3-WT.intersect.broadPeaks.bed -b eed/EED.merged.peaks.V1.bed |awk '$1!~"random|GL|Un"' >H3K27me3.EED.intersect.bed; 
bedtools intersect -v -a /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/result/H3K27me3-WT.intersect.broadPeaks.bed -b eed/EED.merged.peaks.V1.bed |awk '$1!~"random|GL|Un"' >H3K27me3.unique.bed3;
bedtools intersect -v -a eed/EED.merged.peaks.V1.bed -b /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/result/H3K27me3-WT.intersect.broadPeaks.bed |awk '$1!~"random|GL|Un"' >EED.unique.bed3;
awk '$5>0' H3K27me3.EED.intersect.bed |awk -v OFS="" -v ORS="" '{print $1"\t";if($2<=$5){print $5"\t";}else{print $2"\t";}if($3<=$6){print $3"\n"}else{print $6"\n"}}' >H3K27me3.EED.shareRegion.bed3; 
ls *.bed3|while read file;do prefix=$(echo $file|sed 's/.bed3//');/rd1/user/liym/nucleosomeTurnover/scripts/FileForNTR.sh $file week0.bw week1.bw week2.bw week4.bw week6.bw week8.bw ${prefix}.H2BGFP.tsv;Rscript /rd1/user/liym/nucleosomeTurnover/scripts/NTR.R -i=${prefix}.H2BGFP.tsv -p=0.0001 -s=4 -e=9 -o=${prefix}.NTR.tsv 2>NTR.err;done;
bedtools intersect -v -a /rd1/user/liym/nucleosomeTurnover/ChIP-seq/others/suz12/suz12_broad_v2/suz12_peaks.broadPeak -b ../fig2/H3K27me3.peaks.NTR.tsv >SUZ12.unique.vsK27me3.bed3 
bedtools intersect -wao -a ../fig2/H3K27me3.peaks.NTR.tsv -b /rd1/user/liym/nucleosomeTurnover/ChIP-seq/others/suz12/suz12_broad_v2/suz12_peaks.broadPeak|cut -f1-3,10,12- >H3K27me3.SUZ12.intersect.bed
awk '$6>0' H3K27me3.SUZ12.intersect.bed|awk -v OFS="" -v ORS="" '{print $1"\t";if($2<=$6){print $6"\t";}else{print $2"\t";}if($3<=$7){print $3"\n"}else{print $7"\n"}}' |awk '$1!~"random|GL|Un"' >H3K27me3.SUZ12.shareRegion.bed3 
/rd1/user/liym/nucleosomeTurnover/scripts/FileForNTR.sh H3K27me3.SUZ12.shareRegion.bed3 week0.bw week1.bw week2.bw week4.bw week6.bw week8.bw H3K27me3.SUZ12.shareRegion.H2BGFP.tsv
Rscript /rd1/user/liym/nucleosomeTurnover/scripts/NTR.R -i=H3K27me3.SUZ12.shareRegion.H2BGFP.tsv -p=0.0001 -s=4 -e=9 -o=H3K27me3.SUZ12.shareRegion.NTR.tsv 
/rd1/user/liym/nucleosomeTurnover/scripts/FileForNTR.sh SUZ12.unique.vsK27me3.bed3 week0.bw week1.bw week2.bw week4.bw week6.bw week8.bw SUZ12.unique.vsK27me3.H2BGFP.tsv
Rscript /rd1/user/liym/nucleosomeTurnover/scripts/NTR.R -i=SUZ12.unique.vsK27me3.H2BGFP.tsv -p=0.0001 -s=4 -e=9 -o=SUZ12.unique.vsK27me3.NTR.tsv
###2.4.3 enhancers intersect with mutiple regulators
myFunction(){
    bedtools intersect -wao -a ../enhancers.NTR.bed3+ -b ${prefix}.merged.peaks.${version}.bed >${prefix}.enhancers.intersect.${version}.bed; 
    awk '$6>0' ${prefix}.enhancers.intersect.${version}.bed |awk -v OFS="" -v ORS="" '{print $1"\t";if($2<=$6){print $6"\t";}else{print $2"\t";}if($3<=$7){print $3"\n"}else{print $7"\n"}}' |awk '$1!~"random|GL|Un"' >${prefix}.enhancers.shareRegion.${version}.bed3;
    awk '$6>0' ${prefix}.enhancers.intersect.${version}.bed |awk -v OFS="" -v ORS="" '{print $1"\t";if($2<=$6){print $2"\t";}else{print $6"\t";}if($3<=$7){print $7"\n"}else{print $3"\n"}}' |awk '$1!~"random|GL|Un"' >${prefix}.enhancers.mergeRegion.${version}.bed3;
    bedtools intersect -v -a ${prefix}.merged.peaks.${version}.bed -b ../enhancers.NTR.bed3+ |cut -f1-3| awk '$1!~"random|GL|Un"'>${prefix}.unique.${version}.bed3;
    ls *.bed3|while read file;do prefix=$(echo $file|sed 's/.bed3//');/rd1/user/liym/nucleosomeTurnover/scripts/FileForNTR.sh $file ../week0.bw ../week1.bw ../week2.bw ../week4.bw ../week6.bw ../week8.bw ${prefix}.H2BGFP.tsv;Rscript /rd1/user/liym/nucleosomeTurnover/scripts/NTR.R -i=${prefix}.H2BGFP.tsv -p=0.0001 -s=4 -e=9 -o=${prefix}.NTR.tsv 2>NTR.err;done;
}
mkdir eed hdac1 hdac2 p300 suz12
cd suz12
cut -f1-3 /rd1/user/liym/nucleosomeTurnover/ChIP-seq/others/suz12/suz12_broad_v2/suz12_peaks.broadPeak |bedtools intersect -wao -a enhancers.NTR.bed3+ -b stdin >suz12.enhancers.intersect.V2.bed
prefix=suz12
version=V2
myFunction
cd p300
/rd1/user/liym/nucleosomeTurnover/scripts/macs2Peaks.union.sh narrow /rd1/user/liym/nucleosomeTurnover/ChIP-seq/TF/ENCODE/p300_rep1/p300_peaks.narrowPeak /rd1/user/liym/nucleosomeTurnover/ChIP-seq/TF/ENCODE/p300_rep2/p300_peaks.narrowPeak p300.merged.peaks.V0.bed
prefix=p300
version=V0
myFunction
cd eed
/rd1/user/liym/nucleosomeTurnover/scripts/macs2Peaks.union.sh broad /rd1/user/liym/nucleosomeTurnover/ChIP-seq/others/EED/EED_rep1_broad_v1/EED_rep1_peaks.broadPeak /rd1/user/liym/nucleosomeTurnover/ChIP-seq/others/EED/EED_rep2_broad_v1/EED_rep2_peaks.broadPeak EED.merged.peaks.V1.bed 
prefix=EED
version=V1
myFunction
cd hdac1
cut -f1-3 /rd1/user/liym/nucleosomeTurnover/ChIP-seq/others/HDAC1/EEDhete_broad_v1/HDAC1_peaks.broadPeak >HDAC1.merged.peaks.V1.bed
prefix=HDAC1
version=V1
myFunction
cd hdac2
cut -f1-3 /rd1/user/liym/nucleosomeTurnover/ChIP-seq/others/HDAC2/EEDhete_broad_v1/HDAC2_peaks.broadPeak >HDAC2.merged.peaks.V1.bed
prefix=HDAC2
version=V1
myFunction
###2.4.4 enhancers vs regulators number
bedtools multiinter -i <(sort -k1,1 -k2,2n enhancers.NTR.bed3+) <(sort -k1,1 -k2,2n eed/EED.merged.peaks.V1.bed) <(sort -k1,1 -k2,2n hdac1/HDAC1.merged.peaks.V1.bed) <(sort -k1,1 -k2,2n hdac2/HDAC2.merged.peaks.V1.bed) <(sort -k1,1 -k2,2n p300/p300.merged.peaks.V0.bed) <(sort -k1,1 -k2,2n /rd1/user/liym/nucleosomeTurnover/ChIP-seq/others/suz12/suz12_broad_4/suz12_peaks.broadPeak|cut -f1-3) >enhancers.coBind.multiinter.bed 
awk '$5~"1" && $5~","' enhancers.coBind.multiinter.bed |cut -f1-4 >enhancers.coBind.bed4 
/rd1/user/liym/nucleosomeTurnover/scripts/FileForNTR.sh enhancers.coBind.bed4 week0.bw week1.bw week2.bw week4.bw week6.bw week8.bw enhancers.coBind.H2BGFP.tsv
grep -v "NA" enhancers.coBind.H2BGFP.tsv >tmp;mv tmp enhancers.coBind.H2BGFP.tsv
Rscript /rd1/user/liym/nucleosomeTurnover/scripts/NTR.R -i=enhancers.coBind.H2BGFP.tsv -p=0.0001 -s=4 -e=9 -o=enhancers.coBind.NTR.tsv 
bedtools intersect -wo -f 1.0 -r -a enhancers.coBind.bed4 -b enhancers.coBind.NTR.tsv|cut -f1-4,14 >tmp;mv tmp enhancers.coBind.NTR.tsv
###2.4.5 EEDKO vs. WT
####2.4.5.1 Prepare data
ls ../../3/danposRst/result/pooled/*.wig |grep -v "Input"|while read file;do prefix=$(basename $file|sed 's/.final.bgsub.Fnor.smooth.wig//');wigToBigWig -clip $file /mnt/share/liym/data/chr.size/mm10.chrom.sizes ${prefix}.bw 2>wigTobw.log & done;
ls ../../4/danposRst/result/pooled/*.wig|grep -v "Input"|while read file;do prefix=$(basename $file|sed 's/.final.bgsub.Fnor.smooth.wig//');wigToBigWig -clip $file /mnt/share/liym/data/chr.size/mm10.chrom.sizes ${prefix}.bw 2>wigTobw.log; done;
####2.4.5.2 H3K27ac
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
#####All peaks together
cat EEDko.specific.peaks.bed3 WT.EEDko.common.peaks.bed3 WT.specific.peaks.bed3 >H3K27ac/H3K27ac.allPeaks.bed3
cd H3K27ac
ls /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/*H3K27ac*.bw |while read file;do prefix=$(basename $file|sed 's/.subtract.bw//');bwtool summary H3K27ac.allPeaks.bed3 $file peaks.${prefix}.signal.tsv -skip-median; done;
ls ../*H2BGFP*bw |while read file;do prefix=$(basename $file|sed 's/.bw//');bwtool summary H3K27ac.allPeaks.bed3 $file peaks.${prefix}.H2BGFP.tsv -skip-median;done;
paste *.signal.tsv|cut -f1-3,8,16,24,32|awk -v OFS="\t" '{print $1,$2,$3,($4+$5)/2-($6+$7)/2}' >peaks.ko.wt.signals.diff.tsv 
paste *H2BGFP.tsv|cut -f1-3,8,16,24,32,40,48,56|awk -v OFS="\t" '{print $1,$2,$3,($4+$5)/2-($7+$8)/2,$6-($9+$10)/2}' |awk -v OFS="\t" '{print $1,$2,$3,$5-$4}' >peaks.ko.wt.NTR.diff.tsv 
bedtools intersect -wo -a peaks.ko.wt.signals.diff.tsv -b peaks.ko.wt.NTR.diff.tsv|cut -f1-4,8 >peaks.ko.wt.signal.NTR.diff.tsv
paste *.signal.tsv|cut -f1-3,8,16,24,32|awk -v OFS="\t" '{if($6+$7!=0){print $1,$2,$3,($4+$5)/($6+$7)}else{print $1,$2,$3,"inf"}}'|bedtools intersect -wo -a stdin -b peaks.ko.wt.NTR.diff.tsv|cut -f1-4,8 >peaks.ko.wt.signalFC.NTR.diff.tsv
paste *.signal.tsv|cut -f1-3,8,16,24,32|bedtools intersect -wo -a stdin -b peaks.ko.wt.NTR.diff.tsv|cut -f1-7,11 >peaks.ko.wt.repSignal.NTR.diff.tsv
ls /rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/*.bw |while read file;do prefix=$(basename $file|sed 's/.bw//');cut -f1-3 peaks.ko.wt.signalFC.NTR.diff.tsv|bwtool summary -skip-median stdin $file peaks.${prefix}.NC.tsv;done;
paste *NC.tsv|cut -f1-3,8,16,24,32 >peaks.wt.ko.NC.tsv 
bwtool agg 5000:5000 H3K27ac.allPeaks.bed3 /rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDhetoMF-lib1610.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDhetoMF-lib1612.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDkoMF-lib1611.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDkoMF-lib1613.bw peaks.NC.agg.txt
bwtool agg 5000:5000 ../H3K27ac-EEDko.intersect.broadPeaks.bed /rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDhetoMF-lib1610.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDhetoMF-lib1612.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDkoMF-lib1611.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDkoMF-lib1613.bw H3K27ac-EEDko.agg.txt
bwtool agg 5000:5000 ../H3K27ac-WT.intersect.broadPeaks.bed /rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDhetoMF-lib1610.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDhetoMF-lib1612.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDkoMF-lib1611.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDkoMF-lib1613.bw H3K27ac-WT.agg.txt

awk '$4>0.5' peaks.ko.wt.signalFC.NTR.diff.tsv|cut -f1-3 >peaks.ko.wt.signalFC.gt0.5.bed3
bwtool agg 5000:5000 peaks.ko.wt.signalFC.gt0.5.bed3 /rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDhetoMF-lib1610.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDhetoMF-lib1612.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDkoMF-lib1611.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDkoMF-lib1613.bw peaks.ko.wt.signalFC.gt0.5.NCagg.txt
bedtools intersect -f 1.0 -r -a peaks.wt.ko.NC.tsv -b peaks.ko.wt.signalFC.gt0.5.bed3 >peaks.ko.wt.signalFC.gt0.5.NC.tsv 
awk '$4>1.5' peaks.ko.wt.signalFC.NTR.diff.tsv|cut -f1-3 >peaks.ko.wt.signalFC.bed3
bwtool agg 5000:5000 peaks.ko.wt.signalFC.gt1.5.bed3 /rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDhetoMF-lib1610.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDhetoMF-lib1612.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDkoMF-lib1611.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDkoMF-lib1613.bw peaks.ko.wt.signalFC.gt1.5.NCagg.txt

####2.4.5.3 H3K27me3
ln -s /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/result/H3K27me3-EEDko.intersect.broadPeaks.bed H3K27me3-EEDko.peak.bed
ln -s /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/result/H3K27me3-WT.intersect.broadPeaks.bed H3K27me3-WT.peak.bed
#####Peak category
bedtools intersect -wao -a H3K27me3-WT.peak.bed -b H3K27me3-EEDko.peak.bed |awk '$5<0' |cut -f1-3 >EEDko.loss.peak.bed3
bedtools intersect -wao -a H3K27me3-WT.peak.bed -b H3K27me3-EEDko.peak.bed |awk '$5>0' |awk -v OFS="" -v ORS="" '{print $1"\t";if($2<=$5){print $5"\t"}else{print $2"\t"}if($3<=$6){print $3"\n"}else{print $6"\n"}}' >EEDko.common.peak.bed3
 ls /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/1_H3K27me3*.bw|while read file;do prefix=$(basename $file|sed 's/.adaF.subtract.bw//');bwtool summary -skip-median EEDko.common.peak.bed3 $file EEDko.common.${prefix}.signal.tsv;done;
paste EEDko.common.*.tsv|cut -f1-3,8,16,24,32,40|awk -v OFS="\t" '{print $1,$2,$3,($4+$5)/2,($6+$7+$8)/3}'|awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$4/$5}' >EEDko.common.signals.ko.wt.FC.tsv 
cat EEDko.loss.peak.bed3 <(awk '$6<0.5' EEDko.common.signals.ko.wt.FC.tsv|cut -f1-3) >EEDko.loss.final.bed3
awk '$6>=0.5 && $6<=2' EEDko.common.signals.ko.wt.FC.tsv |cut -f1-3 >EEDko.unchange.final.bed3
cat <(awk '$6>2' EEDko.common.signals.ko.wt.FC.tsv|cut -f1-3) <(bedtools intersect -wo -v -a H3K27me3-EEDko.p     eak.bed -b H3K27me3-WT.peak.bed) >EEDko.gain.final.bed3 
#####NTR
ls ../*.bw |while read file;do prefix=$(basename $file|sed 's/.bw//');bwtool summary -skip-median EEDko.unchange.final.bed3 $file EEDko.unchange.${prefix}.tsv;bwtool summary -skip-median EEDko.loss.final.bed3 $file EEDko.loss.final.${prefix}.tsv;done;
paste EEDko.unchange.*H2BGFP*.tsv |cut -f1-3,8,16,24,32,40,48,56|awk -v OFS="\t" '{print $1,$2,$3,($4+$5)/2-($7+$8)/2,$6-($9+$10)/2}' >EEDko.unchange.wt.ko.NTR.bed3+
paste EEDko.loss.*H2BGFP*.tsv |cut -f1-3,8,16,24,32,40,48,56|awk -v OFS="\t" '{print $1,$2,$3,($4+$5)/2-($7+$8)/2,$6-($9+$10)/2}' >EEDko.loss.wt.ko.NTR.bed3+
paste EEDko.gain.*H2BGFP*.tsv |cut -f1-3,8,16,24,32,40,48,56|awk -v OFS="\t" '{print $1,$2,$3,($4+$5)/2-($7+$8)/2,$6-($9+$10)/2}' |awk -v OFS="\t" '{print $1,$2,$3,$5-$4}'>EEDko.gain.wt.ko.NTRdiff.bed3+
#####ChIP-seq signals
paste *gain*signal.tsv|cut -f1-3,8,16,24,32,40|awk -v OFS="\t" '{if($6+$7+$8!=0){print $1,$2,$3,(($4+$5)/2)/(($6+$7+$8)/3)}else{print $1,$2,$3,"inf"}}' >EEDko.gain.ko.wt.FC.tsv
paste *loss*signal.tsv|cut -f1-3,8,16,24,32,40|awk -v OFS="\t" '{if($6+$7+$8!=0){print $1,$2,$3,(($4+$5)/2)/(($6+$7+$8)/3)}else{print $1,$2,$3,"inf"}}' >EEDko.loss.ko.wt.FC.tsv
paste *unchange*signal.tsv|cut -f1-3,8,16,24,32,40|awk -v OFS="\t" '{if($6+$7+$8!=0){print $1,$2,$3,(($4+$5)/2)/(($6+$7+$8)/3)}else{print $1,$2,$3,"inf"}}' >EEDko.unchange.ko.wt.FC.tsv
#####Combined all peaks together
cat <(paste *unchange*signal.tsv|cut -f1-3,8,16,24,32,40) <(paste *gain*signal.tsv|cut -f1-3,8,16,24,32,40) <(paste *loss*signal.tsv|cut -f1-3,8,16,24,32,40) >all.peaks.ko.wt.repSignal.tsv 
bedtools intersect -wo -f 1.0 -r -a all.peaks.ko.wt.repSignal.tsv -b <(cat EEDko.gain.wt.ko.NTRdiff.bed3+ <(awk -v OFS="\t" '{print $1,$2,$3,$5-$4}' EEDko.loss.wt.ko.NTR.bed3+) <(awk -v OFS="\t" '{print $1,$2,$3,$5-$4}' EEDko.unchange.wt.ko.NTR.bed3+))|cut -f1-8,12 >peaks.ko.wt.repSignal.NTR.diff.tsv
ls /rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/*.bw |while read file;do prefix=$(basename $file
|sed 's/.bw//');cut -f1-3 peaks.ko.wt.signalFC.NTR.diff.tsv|bwtool summary -skip-median stdin $file peaks.${prefix}.NC.tsv;done;
paste *NC.tsv|cut -f1-3,8,16,24,32 >peaks.wt.ko.NC.tsv
cut -f1-3 peaks.ko.wt.signalFC.NTR.diff.tsv >peaks.All.bed3
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
#####Remove regions with extreme high values 
bedtools intersect -v -a <(awk '$4<0.5' peaks.ko.wt.SignalsFCminusToZero.NTR.diff.tsv) -b <(awk '$4>5 && $5>5 && $6>5 && $7>5' peaks.wt.ko.NC.tsv) |cut -f1-3 >peaks.ko.wt.signalsFC.lt0.5.bed
bwtool agg 5000:5000 peaks.ko.wt.signalsFC.lt0.5.bed /rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDhetoMF-lib1610.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDhetoMF-lib1612.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDkoMF-lib1611.bw,/rd1/user/liym/nucleosomeTurnover/MNase/danpos_rst/pooled/1_MNase-EEDkoMF-lib1613.bw peaks.ko.wt.signalsFC.lt0.5.NCagg.txt
#####WT NTR binned by signals
mkdir WT-signals && cd WT-signals
ls /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/1_H3K27me3-WT*.bw|while read file;do prefix=$(basename $file|sed 's/.adaF.subtract.bw//');bwtool summary -skip-median ../H3K27me3-WT.peak.bed $file WT.${prefix}.signal.tsv;done;
paste *.tsv|cut -f1-3,8,16,24 >WT.peaks.signals.tsv
sh /rd1/user/liym/nucleosomeTurnover/scripts/FileForNTR.sh ../H3K27me3-WT.peak.bed ../../../fig4/week0.bw ../../../fig4/week1.bw ../../../fig4/week2.bw ../../../fig4/week4.bw ../../../fig4/week6.bw ../../../fig4/week8.bw WT.peaks.H2BGFP.tsv 
bedtools intersect -f 1.0 -r -wo -a WT.peaks.signals.tsv -b ../../../fig2/H3K27me3.peaks.NTR.tsv|awk -v OFS="\t" '{print $1,$2,$3,($4+$5+$6)/3,$16}' >WT.peaks.signals.NTR.tsv
#Peaks from merged bam files
ls /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/1_H3K27me3-WT*.bw|while read file;do prefix=$(basename $file|sed 's/.adaF.subtract.bw//');cut -f1-3 /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/1_H3K27me3-WT-AdultCM_merged_broad/1_H3K27me3-WT-AdultCM_peaks.broadPeak|bwtool summary -skip-median stdin $file WT.merged.peaks.${prefix}.signal.tsv;done;
sh /rd1/user/liym/nucleosomeTurnover/scripts/FileForNTR.sh /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/1_H3K27me3-WT-AdultCM_merged_broad/1_H3K27me3-WT-AdultCM_peaks.broadPeak ../../../fig4/week0.bw ../../../fig4/week1.bw ../../../fig4/week2.bw ../../../fig4/week4.bw ../../../fig4/week6.bw ../../../fig4/week8.bw WT.merged.peaks.H2BGFP.tsv
Rscript /rd1/user/liym/nucleosomeTurnover/scripts/NTR.R -i=WT.merged.peaks.H2BGFP.tsv -p=0.0001 -s=4 -e=9 -o=WT.merged.peaks.NTR.tsv
paste WT.merged*signal.tsv |cut -f1-3,8,16,24 >WT.merged.peaks.signals.tsv
bedtools intersect -f 1.0 -r -wo -a WT.merged.peaks.signals.tsv -b WT.merged.peaks.NTR.tsv|awk -v OFS="\t" '{print $1,$2,$3,($4+$5+$6)/3,$16}' >WT.merged.peaks.signals.NTR.tsv
bedtools intersect -wo -f 1.0 -r -a WT.merged.peaks.signals.NTR.tsv -b /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/EEDproject/1_H3K27me3-WT-AdultCM_merged_broad/1_H3K27me3-WT-AdultCM_peaks.broadPeak|sort -k12,12n |cut -f1-5,12 >WT.merged.peaks.signals.NTR.FCincrease.tsv

##2.5 Figure5 banding and sham
cd fig6
###2.5.1 prepare data
ls /rd1/user/liym/nucleosomeTurnover/H2BGFP_All/5/danposRst/result/pooled/*.wig |grep -v "Input"|while read file;do prefix=$(basename $file|sed 's/.final.bgsub.Fnor.smooth.wig//');wigToBigWig -clip $file /mnt/share/liym/data/chr.size/mm10.chrom.sizes ${prefix}.bw 2>wigTobw.log & done;
cat /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/banding-sham/result/H3K27acNarrow/banding_up.H3K27ac.bed /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/banding-sham/result/H3K27acNarrow/banding_specific.H3K27ac.bed |cut -f1-3 >banding.unique.peaks.bed3
Rscript /mnt/share/liym/bin/PeakAnnoGO.R banding.unique.peaks.bed3 banding.unique.peaks.anno.tsv banding.unique.peaks.anno.Pie.pdf banding.unique.peaks.anno.dotplot.pdf
###2.5.2 Aggregate lines for different types of peaks
ls /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/banding-sham/result/*.bed |grep -v "intersect"|while read file;do prefix=$(basename $file|sed 's/.bed//');bwtool agg 2000:2000 ${file} 5_H2BGFP-CM-week2.5-banding_Rep1.bw,5_H2BGFP-CM-week2.5-banding_Rep2.bw,5_H2BGFP-CM-week2.5-sham_Rep1.bw,5_H2BGFP-CM-week2.5-sham_Rep2.bw ${prefix}.H2BGFP.profile.txt & done;
ls *H2BGFP.profile.txt|while read file;do prefix=$(echo $file|sed 's/.txt//');Rscript /mnt/share/liym/bin/lines.R -i=$file -x="Peak relatively position" -y="Normalized H2BGFP occupancy" -c="blue,darkblue,red,brown" -l="banding_rep1,banding_rep2,sham_rep1,sham_rep2" -o=${prefix}.pdf;done;
ls *.H2BGFP.profile.txt|while read file;do 
    prefix=$(echo $file|sed 's/.txt//');
    awk -v OFS="\t" '{print $1,($2+$3)/2,($4+$5)/2}' $file|Rscript /mnt/share/liym/bin/lines.R -y1=0.15 -y2=0.28 -x="Peak relatively position" -y="Normalized H2BGFP occupancy" -c="blue,red" -l="banding,sham" -o=${prefix}.repMean.pdf
    done;
###2.5.3 NTR for different types of peaks
mkdir NTR
ls /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/banding-sham/result/*.bed |grep -v "intersect"|while read file;do 
    prefix=$(basename $file|sed 's/.bed//');
    ls *.bw |sed 's/.bw//'|while read bw;do 
        cut -f1-3 $file >tmp.bed;bwtool summary tmp.bed ${bw}.bw NTR/${prefix}.${bw}.tsv -skip-median
    done;
done;
cd NTR
ls *5_H2BGFP-CM-week2.5-banding_Rep1.tsv |sed 's/.5_H2BGFP-CM-week2.5-banding_Rep1.tsv//'|while read file;do paste <(cut -f1-3,8 ${file}.5_H2BGFP-CM-week2.5-banding_Rep1.tsv) <(cut -f8 ${file}.5_H2BGFP-CM-week2.5-banding_Rep2.tsv) <(cut -f8 ${file}.5_H2BGFP-CM-week2.5-sham_Rep1.tsv) <(cut -f8 ${file}.5_H2BGFP-CM-week2.5-sham_Rep2.tsv) |awk -v OFS="\t" '{print $1":"$2"-"$3,$4,$5,$6,$7}' >${file}.summary.tsv ;done;
Rscript run.R .H3K27me3.summary.tsv H3K27me3.boxplot.pdf 
Rscript run.R .H3K4me3.summary.tsv H3K4me3.boxplot.pdf 
Rscript run.R .H3K27ac.summary.tsv H3K27me3.boxplot.pdf 
Rscript run.R .H3K4me1.summary.tsv H3K4me3.boxplot.pdf 
mkdir H3K27acNarrow && cd H3K27acNarrow
ls /rd1/user/liym/nucleosomeTurnover/ChIP-seq/histone/banding-sham/result/H3K27acNarrow/*.bed|grep -v "intersect"|while read file;do
 prefix=$(basename $file|sed 's/.bed//');
    ls *.bw |sed 's/.bw//'|while read bw;do 
        cut -f1-3 $file >tmp.bed;bwtool summary tmp.bed ${bw}.bw NTR/${prefix}.${bw}.tsv -skip-median
    done;
done;
ls *5_H2BGFP-CM-week2.5-banding_Rep1.tsv |sed 's/.5_H2BGFP-CM-week2.5-banding_Rep1.tsv//'|while read file;do paste <(cut -f1-3,8 ${file}.5_H2BGFP-CM-week2.5-banding_Rep1.tsv) <(cut -f8 ${file}.5_H2BGFP-CM-week2.5-banding_Rep2.tsv) <(cut -f8 ${file}.5_H2BGFP-CM-week2.5-sham_Rep1.tsv) <(cut -f8 ${file}.5_H2BGFP-CM-week2.5-sham_Rep2.tsv) |awk -v OFS="\t" '{print $1":"$2"-"$3,$4,$5,$6,$7}' >${file}.summary.tsv ;done;
####Start in R
library("reshape2")
library("ggplot2")
args <- commandArgs(TRUE)
files=list.files(pattern=".H3K27ac.summary.tsv")
class=sub(".H3K27ac.summary.tsv","",files)
for(i in 1:length(files)){
    data=read.delim(file=files[i],header=F)
    data<-cbind(rowMeans(data[,c(2,3)]),rowMeans(data[,c(4,5)]))
    colnames(data)<-c("banding","sham")
    data<-melt(data)[,c(2,3)]
    data<-cbind(replicate(nrow(data),class[i]),data)
    colnames(data)<-c("class","condition","value")
    if(exists("final")){
        final=rbind(final,data)
    }else{
        final=data
    }
}
final[,3]<-final[,3]*(-1)
pdf(file="H3K27ac.narrowPeak.boxplot.pdf")
final[final$class=="banding_down",1]="banding_specific"
final[final$class=="banding_up",1]="sham_specific"
pvalue=0
pvalue[1]<-wilcox.test(final[final$class=="banding_sham.unchange" & final$condition=="banding",3],final[final$class=="banding_sham.unchange" & final$condition=="sham",3])$p.value
pvalue[2]<-wilcox.test(final[final$class=="banding_specific" & final$condition=="banding",3],final[final$class=="banding_specific" & final$condition=="sham",3])$p.value
pvalue[3]<-wilcox.test(final[final$class=="sham_specific" & final$condition=="banding",3],final[final$class=="sham_specific" & final$condition=="sham",3])$p.value
pvalue=round(pvalue,digits=3)
ggplot(final, aes(x = class, y = value, fill = condition)) + geom_boxplot(outlier.shape = NA,notch = TRUE) + ylim(-0.5,0) +labs(x="",y="Relative nucleosome turnover")+theme(axis.text.x = element_text(angle = 90))+annotate(geom="text",label=pvalue[1],x=1,y=0)+annotate(geom="text",label=pvalue[2],x=2,y=0)+annotate(geom="text",label=pvalue[3],x=3,y=0)
dev.off()
####End in R
