##Run bash deduplicate.sh ${gff3file} ${TSVfile} ${faafile} ${cdsfile} to get deduplicate file after processing interproscan.sh step
##But only a little duplicate genes predicted by BRAKER2, so this step is optional

gff3=$1
TSV=$2
faa=$3
cds=$4
cat $1|awk '{if($3 == "gene")print $0}'| gawk 'BEGIN {OFS="\t"} {split($1, a, "[:\\-+]"); print a[1], (a[2]+$4), (a[2]+$5),  $9,$8,$7 }'|sort -k1,1 -k2,2n|bedtools merge -s -c 4 -o collapse -delim "|"   >tmp_merged
duplicate=`cat tmp_merged|grep '|'|awk '{print $4}'`
touch ${gff3%.gff3}_duplicate.list
for line in $duplicate
do
#Get the gene list for overlap genes
        gene=`echo "$line"|sed -e 's/ID=//g'|sed -e 's/;|/ /g'|sed -e 's/;//g'`
        IFS=' ' read -r -a gene_list <<< "$gene"
#Compare the cds of each gene
                        for gene1 in "${gene_list[@]}"
                        do
                                        grep -E `echo $(paste <(echo "$gene1") <(echo  "\\.") -d '')` $2  >>tmp.tsv
                        done
                        if [ -s tmp.tsv ]; then
                                        awk -F '\t' ' {print $3, $7, $8, $12, $1}' tmp.tsv|sed -e 's/.t.*//'|grep -v '-'|sort -k5,5 -k4,4 -k1,1|sort -k5,5 -k4,4 -u|sort -k1,1n >tmp.tsv1
                                        tmp_ID1=`cat tmp.tsv1|grep ${gene_list[0]}|awk '{print $4}'`
                                        tmp_IPR1=`echo $tmp_ID1|wc -w`
                                        tmp_ID2=`cat tmp.tsv1|grep ${gene_list[1]}|awk '{print $4}'`
                                        tmp_IPR2=`echo $tmp_ID2|wc -w`
                                                if [ $tmp_IPR1 == 0 && $tmp_IPR2 >0 ]; then
                                                                echo ${gene_list[0]} >>${gff3%.gff3}_duplicate.list
                                                elif [ $tmp_IPR2 == 0 && $tmp_IPR1 >0]; then
                                                                echo ${gene_list[1]} >>${gff3%.gff3}_duplicate.list
                                                elif [[( $tmp_IPR1 < $tmp_IPR2 && $tmp_IPR1 > 0 && $tmp_IPR2 > 0)]]; then
                                                        for i in $tmp_ID1; do
                                                                if ! echo "$tmp_ID2" |grep -q "$i" ; then
                                                                        changed=1
                                                        fi
                                                        done
                                                        if [ "$changed" ]; then
                                                                echo "Condtion1, keep both"
                                                                echo $tmp_ID1
                                                                echo $tmp_ID2
                                                        else
                                                                echo "Condition1 overlap, keep ${gene_list[1]}"
                                                                echo ${gene_list[0]} $tmp_ID1
                                                                echo ${gene_list[1]} $tmp_ID2
                                                                echo ${gene_list[0]}>>${gff3%.gff3}_duplicate.list
                                                        fi
                                                elif [[( $tmp_IPR1 > $tmp_IPR2 && $tmp_IPR1 > 0 && $tmp_IPR2 > 0)]]; then
                                                        for i in $tmp_ID2; do
                                                                if ! echo "$tmp_ID1" |grep -q "$i" ; then
                                                                        changed=1
                                                        fi
                                                        done
                                                        if [ "$changed" ]; then
                                                                echo "Condition 2, keep both"
                                                                echo ${gene_list[0]} $tmp_ID1
                                                                echo ${gene_list[1]} $tmp_ID2
                                                        else
                                                                echo "Condition 2 overlap, keep ${gene_list[0]}"
                                                                echo ${gene_list[0]} $tmp_ID1
                                                                echo ${gene_list[1]} $tmp_ID2
                                                                echo ${gene_list[1]}>>${gff3%.gff3}_duplicate.list
                                                        fi
                                                else
                                                        echo $tmp_ID1
                                                        echo $tmp_ID2
                                                fi

                                 else
                                        echo "no protein domain identified in both gene model"
                                        cat tmp.tsv ${gff3%.gff3}_duplicate.list >>${gff3%.gff3}_duplicate.list
                        fi
                        rm tmp.tsv
done
rm tmp_merged
for i in `cat "${gff3%.gff3}_deduplicate.list"`; do grep -v "${i};\|${i}.t" $2 >>${gff3%.gff3}_deduplicate.gff3.tmp ; done
sort -k1,1 -k4,4n -k9,9 -u ${gff3%.gff3}_deduplicate.gff3.tmp >${gff3%.gff3}_deduplicate.gff3
rm ${gff3%.gff3}_deduplicate.gff3.tmp
awk '{if($3 == "mRNA") print $9}' ${gff3%.gff3}_deduplicate.gff3 |sed -e 's/ID=//g' -e 's/;Parent.*//g'|sort -k1,1 -u >${gff3%.gff3}_protein_deduplicate.list
seqtk subseq $3 ${gff3%.gff3}_protein_deduplicate.list >${gff3%.gff3}_deduplicate_aa.fasta
for i in `cat "${gff3%.gff3}_deduplicate.list"`; do grep -v "${i}.t" $2 >>${gff3%.gff3}_deduplicate.tsv.tmp ; done
seqtk subseq $4 ${gff3%.gff3}_protein_deduplicate.list >${gff3%.gff3}_deduplicate.cds
sort -k1,9 -u ${gff3%.gff3}_deduplicate.tsv.tmp >${gff3%.gff3}_deduplicate.tsv
rm ${gff3%.gff3}_deduplicate.tsv.tmp
