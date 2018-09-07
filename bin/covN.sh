#!/bin/bash

usage="$0 vcf1_name vcf2_name outputname Nbam refGenome gatkJar"

if [[ $# -ne 6 ]] || [[ ! -f $1 ]] || [[ ! -f $2 ]] || [[ ! -f $4 ]] || [[ ! -f $5 ]] || [[ ! -f $6 ]]
then
    echo -e $usage
    exit 1
fi

vcf1=$1
vcf2=$2
out=$3
bam1=$4
GENOME=$5
GATKJAR=$6

name_vcf1=$(echo $vcf1 | sed "s/.vcf//g")
name_vcf2=$(echo $vcf2 | sed "s/.vcf//g")
name_out=$(echo $out | sed "s/.tsv//g")
vcf2bed --deletions < $vcf2 > ${name_vcf2}_deletions.bed
vcf2bed --snvs < $vcf2 > ${name_vcf2}_snvs.bed
vcf2bed --deletions < $vcf1 > ${name_vcf1}_deletions.bed
vcf2bed --snvs < $vcf1 > ${name_vcf1}_snvs.bed
bedops --everything {${name_vcf2},${name_vcf1}}_{deletions,snvs}.bed | awk 'BEGIN{OFS="\t"}{print($1,$2,$3)}' > ${name_out}.bed
java -Xms512m -Xmx6G -jar $GATKJAR HaplotypeCaller -R $GENOME -I $bam1 -O "$name_out.vcf" --intervals ${name_out}.bed --output-mode EMIT_ALL_SITES > "$name_out.log" 2>&1
cat "$name_out.vcf" | sed "/^#/d" | perl -lane '$F[9]=~s/^[^:]*:([^:]*).*/$1/;@reads=split(",",$F[9]);$reads[1]=="" and $reads[1]=0;if($reads[0] eq "./."){$readsref=0;$readsout=0}else{$readsref=splice(@reads,0,1);$readsout=join(",",@reads)};print join("\t",@F[0,1,3,4],$readsref,$readsout)' > $out
