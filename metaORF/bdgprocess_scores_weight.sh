#!/bin/sh

#Author:Mohamood Adhil, Date:10/10/2017
function usage()
{
    echo "This script average bedgraph files accross the all features for a given base length"
    echo ""
    echo "-h --help"
    echo "bdg input file paths comma seperated (Mandatory) --bdgpaths=$bdgpaths"
    echo "sample name comma seperated (Mandatory) --samplename=$samplename"
    echo "Strand information to be considered (Mandatory yes or no) --strand=$strand"
    echo "Feature bed file (Tab seperated) with no header. If --strand=yes, then the file should contain strand information in the fourth column (chrom,start,end,strand) (Mandatory) --bedfile=$bedfile"
    echo "Number of bases for avaeraging all the feature (Mandatory eg 1000) --avglength=$avglength"
    echo "Number of bases to add as a flanking region (Mandatory eg 500) --flanking=$flanking"
    echo "Threshold for genecount (Mandatory eg 1.5) --gthr=$gthr"
    echo "Output Path (Mandatory) --outpath=$outpath"
}

while [ "$1" != "" ]; do
    PARAM=`echo $1 | awk -F= '{print $1}'`
    VALUE=`echo $1 | awk -F= '{print $2}'`
    case $PARAM in
        -h | --help)
            usage
            exit
            ;;
        --bdgpaths)
            bdgpaths=$VALUE
            ;;
        --samplename)
            samplename=$VALUE
            ;;
        --strand)
            strand=$VALUE
            ;;
        --bedfile)
            bedfile=$VALUE
            ;;
	--avglength)
            avglength=$VALUE
            ;;
        --flanking)
            flanking=$VALUE
            ;;
	--gthr)
            gthr=$VALUE
            ;;
        --outpath)
            outpath=$VALUE
            ;;
        *)
            echo "ERROR: unknown parameter \"$PARAM\""
            usage
            exit 1
            ;;
    esac
    shift
done

#bed graph file into array
arr2=(${bdgpaths//,/ })
#sample name into array
arr3=(${samplename//,/ })

echo "$strand"
if [ "$strand" == "yes" ]; then
#negstrandfile
echo "$bedfile"
awk -v flanking="${flanking}" -F "\t" '{print $1"\t"$2-flanking"\t"$3+flanking"\t"$2"\t"$3"\t"$4}' $bedfile | awk -F "\t" '{if ($2 < 0) print $1"\t""0""\t"$3"\t"$4"\t"$5"\t"$6; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' |  awk -F "\t" '{if ($6 ~ /-1/) print $0}' > $outpath/negstrand_features.txt
#posstrandfile
awk -v flanking="${flanking}" -F "\t" '{print $1"\t"$2-flanking"\t"$3+flanking"\t"$2"\t"$3"\t"$4}' $bedfile | awk -F "\t" '{if ($2 < 0) print $1"\t""0""\t"$3"\t"$4"\t"$5"\t"$6; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' | awk -F "\t" '{if ($6 == 1) print $0}' > $outpath/posstrand_features.txt
for index in ${!arr2[@]}
do
#avg score
bedtools intersect -a ${arr2[$index]} -b $outpath/posstrand_features.txt -wb | awk -F "\t" '{printf("%s\t%s\t%s\t%.6f\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9)}' | awk -F "\t" '{OFS="\t"}{for (i=1;i<=$3-$2;i++) print $2+i,$4,$8,$9,$5$8$9}' | awk -v avglength="${avglength}" -F "\t" '{if ($3 >= $1) print $1-$3"\t"$2"\t"$5; else if (($3 < $1) && ($4 > $1)) print (($1-$3)/($4-$3))*avglength"\t"$2"\t"$5; else if ($4 <= $1) print ($1-$4)+avglength"\t"$2"\t"$5}' > $outpath/pos\_${arr3[$index]}_scores.txt
bedtools intersect -a ${arr2[$index]} -b $outpath/negstrand_features.txt -wb | awk -F "\t" '{printf("%s\t%s\t%s\t%.6f\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9)}' | awk -F "\t" '{OFS="\t"}{for (i=1;i<=$3-$2;i++) print $2+i,$4,$8,$9,$5$8$9}' | awk -v avglength="${avglength}" -F "\t" '{if ($3 >= $1) print $1-$3"\t"$2"\t"$5; else if (($3 < $1) && ($4 > $1)) print (($1-$3)/($4-$3))*avglength"\t"$2"\t"$5; else if ($4 <= $1) print ($1-$4)+avglength"\t"$2"\t"$5}' > $outpath/neg\_${arr3[$index]}_scores.txt
#avg count
bedtools intersect -a ${arr2[$index]} -b $outpath/posstrand_features.txt -wb | awk -F "\t" '{printf("%s\t%s\t%s\t%.6f\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9)}' | awk -v gthr="${gthr}" -F "\t" '{if ($4 >= gthr) print $1"\t"$2"\t"$3"\t""1""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9; else print $1"\t"$2"\t"$3"\t""0""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | awk -F "\t" '{OFS="\t"}{for (i=1;i<=$3-$2;i++) print $2+i,$4,$8,$9,$5$8$9}' | awk -v avglength="${avglength}" -F "\t" '{if ($3 >= $1) print $1-$3"\t"$2"\t"$5; else if (($3 < $1) && ($4 > $1)) print (($1-$3)/($4-$3))*avglength"\t"$2"\t"$5; else if ($4 <= $1) print ($1-$4)+avglength"\t"$2"\t"$5}' > $outpath/pos\_${arr3[$index]}_weight.txt
bedtools intersect -a ${arr2[$index]} -b $outpath/negstrand_features.txt -wb | awk -F "\t" '{printf("%s\t%s\t%s\t%.6f\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9)}' | awk -v gthr="${gthr}" -F "\t" '{if ($4 >= gthr) print $1"\t"$2"\t"$3"\t""1""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9; else print $1"\t"$2"\t"$3"\t""0""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | awk -F "\t" '{OFS="\t"}{for (i=1;i<=$3-$2;i++) print $2+i,$4,$8,$9,$5$8$9}' | awk -v avglength="${avglength}" -F "\t" '{if ($3 >= $1) print $1-$3"\t"$2"\t"$5; else if (($3 < $1) && ($4 > $1)) print (($1-$3)/($4-$3))*avglength"\t"$2"\t"$5; else if ($4 <= $1) print ($1-$4)+avglength"\t"$2"\t"$5}' > $outpath/neg\_${arr3[$index]}_weight.txt
done
rm $outpath/posstrand_features.txt $outpath/negstrand_features.txt
elif [ "$strand" == "no" ]; then
awk -v flanking="${flanking}" -F "\t" '{print $1"\t"$2-flanking"\t"$3+flanking"\t"$2"\t"$3"\t"$4}' $bedfile | awk -F "\t" '{if ($2 < 0) print $1"\t""0""\t"$3"\t"$4"\t"$5"\t"$6; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > $outpath/features.txt
for index in ${!arr2[@]}
do
#avg score
bedtools intersect -a ${arr2[$index]} -b $outpath/features.txt -wb | awk -F "\t" '{printf("%s\t%s\t%s\t%.6f\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9)}' | awk -F "\t" '{OFS="\t"}{for (i=1;i<=$3-$2;i++) print $2+i,$4,$8,$9,$5$8$9}' | awk -v avglength="${avglength}" -F "\t" '{if ($3 >= $1) print $1-$3"\t"$2"\t"$5; else if (($3 < $1) && ($4 > $1)) print (($1-$3)/($4-$3))*avglength"\t"$2"\t"$5; else if ($4 <= $1) print ($1-$4)+avglength"\t"$2"\t"$5}' > $outpath/${arr3[$index]}_scores.txt
#avg count
bedtools intersect -a ${arr2[$index]} -b $outpath/features.txt -wb | awk -F "\t" '{printf("%s\t%s\t%s\t%.6f\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9)}' | awk -v gthr="${gthr}" -F "\t" '{if ($4 >= gthr) print $1"\t"$2"\t"$3"\t""1""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9; else print $1"\t"$2"\t"$3"\t""0""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | awk -F "\t" '{OFS="\t"}{for (i=1;i<=$3-$2;i++) print $2+i,$4,$8,$9,$5$8$9}' | awk -v avglength="${avglength}" -F "\t" '{if ($3 >= $1) print $1-$3"\t"$2"\t"$5; else if (($3 < $1) && ($4 > $1)) print (($1-$3)/($4-$3))*avglength"\t"$2"\t"$5; else if ($4 <= $1) print ($1-$4)+avglength"\t"$2"\t"$5}' > $outpath/${arr3[$index]}_weight.txt
done
rm $outpath/features.txt
elif [ "$strand" != "yes" ] || [ "$strand" = "no" ]; then
echo "Error: The program didn't run, please provide yes or no for --strand"
fi
