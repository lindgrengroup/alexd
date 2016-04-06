source 0_config.sh

$PLINK --file  $INPUT --out $INPUT --make-bed
$PLINK --file  $INPUT --out ${INPUT}_R --recode AD


# for chr in $(seq 1 22); do
#      plink --bfile ~/TRIOS/TRIO_ID_SWAP \
#            --chr $chr \
#            --make-bed \
#            --out ~/TRIOS/IMPUTE/INPUT/TRIO_chr$chr ;
# done

for chr in $(seq 1 22); do
    $PLINK --file  $INPUT \
           --chr $chr \
           --recode AD \
           --out ${INPUT}_R_${chr}
done

for chr in $(seq 1 22); do
    $PLINK --file  $INPUT \
           --chr $chr \
           --make-just-bim \
           --out ${INPUT}_R_${chr}
done



$GTOOL -P --ped $INPUT.ped --map $INPUT.map --og  $INPUT.gen --os  $INPUT.sample

# awk '{print $0 >> "TRIO_chr"$1".gen"}' TRIO.gen