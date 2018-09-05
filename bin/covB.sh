#!/bin/bash

usage="$0 vcf1 vcf2 outputvcf_datavcf1bedvcf2 outputvcf_datavcf2bedvcf1 bam1 bam2 genome gatkjar"

if [[ $# -ne 8 ]] || [[ ! -f $1 ]] || [[ ! -f $2 ]] || [[ ! -f $5 ]] || [[ ! -f $6 ]]
then
    echo -e $usage
    exit 1
fi

vcf1=$1
vcf2=$2
out1=$3
out2=$4
bam1=$5
bam2=$6
GENOME=$7
GATKJAR=$8

dothething ()
{
    bam1=$1
    vcf2=$2
    out=$3

    name_vcf2=$(echo $2 | sed "s/.vcf//g")
    name_out=$(echo $out | sed "s/.tsv//g")
    vcf2bed --deletions < $vcf2 > ${name_vcf2}_deletions.bed
    vcf2bed --snvs < $vcf2 > ${name_vcf2}_snvs.bed
    bedops --everything ${name_vcf2}_{deletions,snvs}.bed | awk 'BEGIN{OFS="\t"}{print($1,$2,$3)}' > ${name_vcf2}.bed
    java -Xms512m -Xmx6G -jar $GATKJAR HaplotypeCaller -R $GENOME -I $bam1 -O "$name_out.vcf" --intervals ${name_vcf2}.bed --output-mode EMIT_ALL_SITES > "$name_out.log" 2>&1
    cat "$name_out.vcf" | sed "/^#/d" | perl -lane '$F[9]=~s/^[^:]*:([^:]*).*/$1/;@reads=split(",",$F[9]);$reads[1]=="" and $reads[1]=0;if($reads[0] eq "./."){$readsref=0;$readsout=0}else{$readsref=splice(@reads,0,1);$readsout=join(",",@reads)};print join("\t",@F[0,1,3,4],$readsref,$readsout)' > $out
}

dothething $bam1 $vcf2 $out1
dothething $bam2 $vcf1 $out2
