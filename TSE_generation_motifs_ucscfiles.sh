#################generate tooth specific enhancers
cd /home/FCAM/ewentworth/cotney/rawdata/chromatinsegmentations/mouse_18state/tissue_specific_enhancer_analysis
module load bedtools
scp -r /home/FCAM/ahardy/ANALYSIS/ChIP/ChromHMM_all_tissues/CHROMHMM_OUTPUT_FILES/18State/*dense*bed .

ls *.bed | sed 's/_rep.*//g; s/_Rep.*//g' | awk '!visited[$0]++' > samples.txt
ls *.bed | sed 's/_rep.*//g; s/_Rep.*//g' | awk '!visited[$0]++' | egrep "fac|Tooth|forebrain|NCC|tube" | egrep -v "postnatal_0_day_forebrain" > CF_samples.txt
egrep -v -f CF_samples.txt samples.txt | egrep "postnatal|embryonic_16.5|unknown" > non-CF_samples.txt


for sample in *.bed
do
export NAME=`echo $sample | sed 's/.bed//g'`
egrep "9_EnhA1|10_EnhA2|8_EnhG2" $sample | cut -f1,2,3 > "$NAME"_enhancers.bed
done


for ROOT in $( cat samples.txt)
do
ls "$ROOT"*enhancers.bed
cat "$ROOT"*enhancers.bed | bedtools sort -i stdin | bedtools merge -i stdin > "$ROOT"_merged_sorted_all_strong_enhancers.bed
done

for ROOT in $( cat non-CF_samples.txt)
do
echo $ROOT
ls "$ROOT"*_merged_sorted_all_strong_enhancers.bed
cat "$ROOT"*_merged_sorted_all_strong_enhancers.bed | bedtools sort -i stdin | bedtools merge -i stdin | wc -l
cat "$ROOT"*_merged_sorted_all_strong_enhancers.bed | bedtools sort -i stdin | bedtools merge -i stdin >> all_strong_enhancer_segments_non-CF_samples.bed
done

bedtools sort -i all_strong_enhancer_segments_non-CF_samples.bed | bedtools merge -i stdin > merged_sorted_all_strong_enhancer_segments_non-CF_samples.bed


for ROOT in $( cat CF_samples.txt)
do
echo $ROOT
echo "$ROOT"*_merged_sorted_all_strong_enhancers.bed
cat "$ROOT"*_merged_sorted_all_strong_enhancers.bed | bedtools sort -i stdin | bedtools merge -i stdin >> all_strong_enhancer_segments_CF_samples.bed
done
cat all_strong_enhancer_segments_CF_samples.bed | bedtools sort -i stdin | bedtools merge -i stdin > merged_sorted_all_strong_enhancer_segments_CF_samples.bed

bedtools intersect -c -a merged_sorted_all_strong_enhancer_segments_CF_samples.bed -b all_strong_enhancer_segments_CF_samples.bed -b merged_sorted_all_strong_enhancer_segments_CF_samples.bed > count_CF_reproducibility.txt

cat count_CF_reproducibility.txt | awk '{if($4>=2) print $0}' > reproducible_merged_sorted_strong_CF_enhancers_2.bed



bedtools intersect -v -a embryonic_13.5_day_Tooth_merged_sorted_all_strong_enhancers.bed -b merged_sorted_all_strong_enhancer_segments_non-CF_samples.bed > strong_novel_embryonic_tooth_enhancers.bed

bedtools intersect -v -a reproducible_merged_sorted_strong_CF_enhancers_2.bed -b merged_sorted_all_strong_enhancer_segments_non-CF_samples.bed | bedtools intersect -v -a stdin -b strong_novel_embryonic_tooth_enhancers.bed > strong_novel_embryonic_CF_enhancers.bed



awk '{print $0, $1"_"$2"_"$3}' strong_novel_embryonic_tooth_enhancers.bed | sed 's/ /\t/g' > temp.txt
mv temp.txt strong_novel_embryonic_tooth_enhancers.bed


cat *merged_sorted_all_strong_enhancers.bed | bedtools sort -i stdin | bedtools merge -i stdin | cut -f1,2,3 | bedtools sort -i stdin > all_strong_enhancers_background.txt


##############################lift these over, then back, for conserved region identification

liftOver -minMatch=0.25 strong_novel_embryonic_tooth_enhancers.bed ~/cotney/genome/mm10/mm10ToHg19.over.chain.gz temp.txt unmapped.bed
sort -k1,1 -k2,2n temp.txt > hg19-25pct-strong_novel_embryonic_tooth_enhancers.bed 

liftOver -minMatch=0.25 hg19-25pct-strong_novel_embryonic_tooth_enhancers.bed ~/cotney/genome/hg19/hg19ToMm10.over.chain.gz temp.txt unmapped.bed
sort -k1,1 -k2,2n temp.txt > conserved_trong_novel_embryonic_tooth_enhancers.bed 


liftOver -minMatch=0.25 strong_novel_embryonic_CF_enhancers.bed ~/cotney/genome/mm10/mm10ToHg19.over.chain.gz temp.txt unmapped.bed
sort -k1,1 -k2,2n temp.txt > hg19-25pct-strong_novel_embryonic_CF_enhancers.bed 

liftOver -minMatch=0.25 hg19-25pct-strong_novel_embryonic_CF_enhancers.bed ~/cotney/genome/hg19/hg19ToMm10.over.chain.gz temp.txt unmapped.bed
sort -k1,1 -k2,2n temp.txt > conserved_strong_novel_embryonic_CF_enhancers.bed 


liftOver -minMatch=0.25 all_strong_enhancers_background.txt ~/cotney/genome/mm10/mm10ToHg19.over.chain.gz temp.txt unmapped.bed
sort -k1,1 -k2,2n temp.txt > hg19-25pct-conserved_all_strong_enhancers_background.txt

liftOver -minMatch=0.25 hg19-25pct-conserved_all_strong_enhancers_background.txt ~/cotney/genome/hg19/hg19ToMm10.over.chain.gz temp.txt unmapped.bed
sort -k1,1 -k2,2n temp.txt > conserved_all_strong_enhancers_background.txt 

cut -f1,2,3 conserved_trong_novel_embryonic_tooth_enhancers.bed | cat - conserved_all_strong_enhancers_background.txt | bedtools sort -i stdin > conserved_all_strong_enhancers_background.bed

scp -r conserved* ~/cotney/analysis/bulkRNA/tooth/WGCNA/.
scp -r all_*.bed all_*.txt ~/cotney/analysis/bulkRNA/tooth/WGCNA/.

for sample in *merged_sorted_all_strong_enhancers.bed
do
liftOver -minMatch=0.25 $sample ~/cotney/genome/mm10/mm10ToHg19.over.chain.gz temp.txt unmapped.bed
sort -k1,1 -k2,2n temp.txt > hg19-25pct-$sample
done
scp -r hg19*.bed /home/FCAM/ewentworth/cotney/rawdata/chromatinsegmentations/mouse_18state/.


scp -r *ooth* *CF*.bed conserved_all_strong_enhancers_background.bed ~/cotney/analysis/bulkRNA/tooth/WGCNA/.

#################################look for overlaps of the TSEs with vista enhancers 
cd /home/FCAM/ewentworth/cotney/analysis/chipseq/mousetooth/VISTA
scp -r /home/FCAM/ewentworth/cotney/rawdata/chromatinsegmentations/mouse_18state/tissue_specific_enhancer_analysis/conserved_trong_novel_embryonic_tooth_enhancers.bed .
scp -r /home/FCAM/ewentworth/cotney/rawdata/chromatinsegmentations/mouse_18state/tissue_specific_enhancer_analysis/embryonic_13.5_day_Tooth_merged_sorted_all_strong_enhancers.bed .

cat all_VISTA_8_2020.txt | grep ">"> VISTA_enhancer_coord_8_2020.txt
cat VISTA_enhancer_coord_8_2020.txt | grep Human > Human_VISTA_enhancer_coord_8_2020.txt
cat VISTA_enhancer_coord_8_2020.txt | grep Mouse > Mouse_VISTA_enhancer_coord_8_2020.txt
cat Human_VISTA_enhancer_coord_8_2020.txt | grep -v negative > Human_positive_VISTA.txt
cat Human_VISTA_enhancer_coord_8_2020.txt | grep negative > Human_negative_VISTA.txt
cat Mouse_VISTA_enhancer_coord_8_2020.txt | grep -v negative > Mouse_positive_VISTA.txt
cat Mouse_VISTA_enhancer_coord_8_2020.txt | grep negative > Mouse_negative_VISTA.txt

sed 's/|/\t/g; s/:/\t/g; s/-/\t/g' Human_positive_VISTA.txt | sort -k1,1 -k2,2n | awk '!visited[$0]++' > temp.txt
mv temp.txt Human_positive_VISTA.txt 
sed 's/|/\t/g; s/:/\t/g; s/-/\t/g' Mouse_positive_VISTA.txt | sort -k1,1 -k2,2n | awk '!visited[$0]++' > temp.txt
mv temp.txt Mouse_positive_VISTA.txt 

cat Human_positive_VISTA.txt Mouse_positive_VISTA.txt | sort -k1,1 -k2,2n > Combined_positive_VISTA.txt


sed 's/|/\t/g; s/:/\t/g; s/-/\t/g' Human_positive_VISTA.txt | egrep "branchial|snout|tooth|face|facial" > Human_cf_positive_VISTA.txt 
sed 's/|/\t/g; s/:/\t/g; s/-/\t/g' Mouse_positive_VISTA.txt | egrep "branchial|snout|tooth|face|facial" > Mouse_cf_positive_VISTA.txt
cat Human_cf_positive_VISTA.txt Mouse_cf_positive_VISTA.txt | awk '!visited[$0]++' > Combined_cf_positive_VISTA.txt

sed -i "s/>//g; s/positive.*//g" Combined_cf_positive_VISTA.txt
sed -i "s/>//g; s/positive.*//g" Combined_positive_VISTA.txt

awk '{print $2, $3, $4, $5"_"$6}' Human_positive_VISTA.txt | sed 's/ /\t/g' | bedtools sort -i stdin > temp2.txt
liftOver -minMatch=0.25 -bedPlus=4 -tab temp2.txt ~/cotney/genome/hg19/hg19ToMm10.over.chain.gz temp.txt unmapped.bed
sort -k1,1 -k2,2n temp.txt | awk '{print $0, "Human"}' > mm10_Human_positive_VISTA.txt

awk '{print $2, $3, $4, $5"_"$6}' Mouse_positive_VISTA.txt | sed 's/ /\t/g' | bedtools sort -i stdin > temp2.txt
liftOver -minMatch=0.25 -bedPlus=4 -tab temp2.txt ~/cotney/genome/mm9/mm9ToMm10.over.chain.gz temp.txt unmapped.bed
sort -k1,1 -k2,2n temp.txt | awk '{print $0, "Mouse"}' > mm10_Mouse_positive_VISTA.txt

#now all the vista coordinates are in mm10. put them together in one file
cat mm10_Mouse_positive_VISTA.txt mm10_Human_positive_VISTA.txt | awk '{print $1, $2, $3, $5"_"$4}' | sed 's/ /\t/g' | bedtools sort -i stdin > mm10_combined_positive_VISTA.txt

# now intersect with the tooth specific enhancers, but pad the distance to either side of the enhancer by 1kbp
awk '{print $1, $2-2000, $3+2000, $4}' mm10_combined_positive_VISTA.txt | sed 's/ /\t/g' > Combined_positive_VISTA_1kb_paddedcoords.txt
#1604 total enhancers that lift over to mm10

bedtools intersect -wb -a conserved_trong_novel_embryonic_tooth_enhancers.bed -b Combined_positive_VISTA_1kb_paddedcoords.txt > strong_novel_embryonic_tooth_enhancers_VISTA_overlap.bed
bedtools intersect -wb -a embryonic_13.5_day_Tooth_merged_sorted_all_strong_enhancers.bed -b Combined_positive_VISTA_1kb_paddedcoords.txt > embryonic_tooth_enhancers_VISTA_overlap.bed



##########################################motif analysis from homer

/home/FCAM/ewentworth/cotney/rawdata/chromatinsegmentations/superenhancers/mouse/E13.5_tooth_se_notindbSuper.bed

/home/FCAM/ewentworth/cotney/rawdata/chromatinsegmentations/superenhancers/mouse/embryonic_13_5_day_Tooth_Rep1_Superenhancer.bed

findMotifsGenome.pl strong_novel_embryonic_tooth_enhancers.bed mm10 tooth_motifs -bg conserved_all_strong_enhancers_background.bed -size 200

findMotifsGenome.pl conserved_trong_novel_embryonic_tooth_enhancers.bed mm10 conserved_tooth_motifs -bg conserved_all_strong_enhancers_background.bed -size 200

findMotifsGenome.pl conserved_strong_novel_embryonic_CF_enhancers.bed mm10 conserved_CF_motifs -bg conserved_all_strong_enhancers_background.bed -size 200

rm Tooth.superenhancers.formotifs.bed

cat embryonic_13_5_day_Tooth_Rep1_Superenhancer.bed | while read line
do
echo $line | sed 's/ /\t/g' > temp.bed
bedtools intersect -a temp.bed -b embryonic_13.5_day_Tooth_merged_sorted_all_strong_enhancers.bed -wb | cut -f1,2,3 >> Tooth.superenhancers.formotifs.bed
done

findMotifsGenome.pl Tooth.superenhancers.formotifs.bed mm10 tooth_SE_motifs -bg conserved_all_strong_enhancers_background.bed -size 200

rm Tooth.superenhancers.formotifs.bed

cat E13.5_tooth_se_notindbSuper.bed | while read line
do
echo $line | sed 's/ /\t/g' > temp.bed
bedtools intersect -a temp.bed -b embryonic_13.5_day_Tooth_merged_sorted_all_strong_enhancers.bed -wb | cut -f1,2,3 >> Tooth.superenhancers.formotifs.bed
done

findMotifsGenome.pl Tooth.superenhancers.formotifs.bed mm10 tooth_novel_SE_motifs -bg conserved_all_strong_enhancers_background.bed -size 200



###################make loop files of novel enhancers:target genes bed file from great
export GENCODE=/home/FCAM/ewentworth/cotney/genome/mm10/mm10_refGene_TSS_coord.bed
awk '{if ($3 == "transcript") {print $0}}' /home/FCAM/ewentworth/cotney/genome/mm10/mm10.refGene.gtf | awk '{if ($7 == "-") {print $1, $5, $5+1, $7, $10} else if ($7 == "+") {print $1, $5, $5-1, $7, $10}}' | sed 's/\"//g; s/\;//g; s/ /\t/g' | sort -k1,1 -k2,2n | awk '!visited[$0]++' > /home/FCAM/ewentworth/cotney/genome/mm10/mm10_refGene_TSS_coord.bed


awk '{print $2, $3, $4, $5, $6, $7, $8}' strong_novel_tooth_enh_genetarget_enhancerloops.txt > strong_novel_tooth_enh_genetarget_enhancerloops.bed

#### in R
gencode<-read.table(file='/home/FCAM/ewentworth/cotney/genome/mm10/mm10_refGene_TSS_coord.bed')
targets<-read.table(file='tooth_target_gene_coords.txt')
enhancers1<-read.table(file='strong_novel_tooth_enh_genetarget_enhancerloops.txt')
enhancers<-enhancers1[,c(2,3,4,5,6,7,8)]
rm<-setdiff(enhancers$V6, gencode$V7)
enhancers<-enhancers[which(enhancers$V7 %in% gencode$V5),]
enhancers<-enhancers[grep("Rik", enhancers$V7, invert=TRUE),]
enhancers$V8<-as.numeric(enhancers$V8)
enhancers$V8<-abs(enhancers$V8)
enhancers=enhancers[order(enhancers[,'V2'], enhancers[,'V3'], enhancers[,'V4'], enhancers[,'V8']),]
enhancers$name<-paste(enhancers$V2, enhancers$V3, enhancers$V4, sep="_")
enhancers<-enhancers[!duplicated(enhancers$name),]
gencode<-gencode[which(gencode$V5 %in% enhancers $V7),]
gencode = gencode[!duplicated(gencode$V5),]
colnames(gencode)<-c('V1', 'V2', 'V3', 'V4', 'genes')
colnames(enhancers)<-c('V1', 'V2', 'V3', 'V4', 'V5', 'genes', 'V6')
combined<-merge(enhancers,gencode, by='genes', all=TRUE)
combined<-combined[order(combined$V1.x, combined$V2.x),]
combined$sourcename<-paste('tooth', rownames(combined), sep='_')
combined$score<-c('1000')
combined$value<-c('1000')
combined$color<-c('0,0,0')
combined$sourcestrand<-c('.')
combined$expname<-c('incisor')
combined$name<-paste(rownames(combined), combined$genes, sep="_")
#combined<-combined[,c(2,3,4,18,13,14,17,15,8,9,10,1,11,2,3,4,12,16)]
combined<-combined[,c(2,3,4,8,14,15,18,16,9,10,11,1,17,2,3,4,13,17)]


write.table(gencode, file="tooth_target_gene_coords.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(enhancers, file="strong_novel_tooth_enh_genetarget_enhancerloops.bed", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(combined, file="strong_novel_tooth_enh_genetarget_enhancerloops.interact", col.names=FALSE, row.names=FALSE, quote=FALSE)

####not in R

sed -i 's/ /\t/g' strong_novel_tooth_enh_genetarget_enhancerloops.interact 
bedToBigBed -type=bed5+13 -as=interact.as  -extraIndex=name -unc strong_novel_tooth_enh_genetarget_enhancerloops.interact /home/FCAM/ewentworth/cotney/genome/mm10/mm10.chrom.sizes /tgc/TGCore_User_Data/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/strong_novel_tooth_enh_genetarget_enhancerloops.bigInteract

chmod -R 755 /tgc/TGCore_User_Data/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/strong_novel_tooth_enh_genetarget_enhancerloops.bigInteract




liftOver -minMatch=0.25 -bedPlus=3 -tab strong_novel_tooth_enh_genetarget_enhancerloops.interact  ~/cotney/genome/mm10/mm10ToHg19.over.chain.gz temp.txt unmapped.bed

awk '{if ($11 > $10) {print $9, $10, $11, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $1, $2, $3, $17, $18} else {print $9, $11, $10, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $1, $2, $3, $17, $18}}' temp.txt | sort -k1,1 -k2,2n | sed 's/ /\t/g' > temp2.txt

liftOver -minMatch=0.25 -bedPlus=3 -tab temp2.txt ~/cotney/genome/mm10/mm10ToHg19.over.chain.gz temp.txt unmapped.bed

awk '{if ($4 = $1) {print $4, $5, $6, $7, $8, $9, $10, $11, $1, $2, $3, $15, $16, $17, $18, $19, $20, $21}}' temp.txt | sort -k1,1 -k2,2n | sed 's/ /\t/g' > hg19-25pct-strong_novel_tooth_enh_genetarget_enhancerloops.interact 

bedToBigBed -type=bed5+13 -as=interact.as  -extraIndex=name -unc hg19-25pct-strong_novel_tooth_enh_genetarget_enhancerloops.interact  /home/FCAM/ewentworth/cotney/genome/hg19/hg19.chrom.sizes /tgc/TGCore_User_Data/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/hg19-25pct-strong_novel_tooth_enh_genetarget_enhancerloops.bigInteract 
chmod -R 755 /tgc/TGCore_User_Data/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/hg19-25pct-strong_novel_tooth_enh_genetarget_enhancerloops.bigInteract 






















cut -f1,2,3 hg19-25pct-embryonic_13.5_day_Tooth_merged_sorted_all_strong_enhancers.bed > hg19-25pct-embryonic_13.5_day_Tooth_merged_sorted_all_strong_enhancers.bed.prep
bedToBigBed -type=bed3 hg19-25pct-embryonic_13.5_day_Tooth_merged_sorted_all_strong_enhancers.bed.prep /home/FCAM/ewentworth/cotney/genome/hg19/hg19.chrom.sizes /tgc/TGCore_User_Data/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/hg19-25pct-embryonic_13.5_day_Tooth_merged_sorted_all_strong_enhancers.bigBed

chmod -R 755 /tgc/TGCore_User_Data/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/hg19-25pct-embryonic_13.5_day_Tooth_merged_sorted_all_strong_enhancers.bigBed



cat hg19-25pct-embryonic_13.5_day_Tooth_rep1_state*.bed | bedtools sort -i stdin > hg19-25pct-embryonic_13.5_day_Tooth_rep1.dense.bed
bedToBigBed -type=bed9 hg19-25pct-embryonic_13.5_day_Tooth_rep1.dense.bed /home/FCAM/ewentworth/cotney/genome/hg19/hg19.chrom.sizes /tgc/TGCore_User_Data/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/hg19-25pct-embryonic_13.5_day_Tooth_rep1.dense.bigBed
chmod -R 755 /tgc/TGCore_User_Data/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/hg19-25pct-embryonic_13.5_day_Tooth_rep1.dense.bigBed

bedToBigBed -type=bed3+1 strong_novel_embryonic_tooth_enhancers.bed /home/FCAM/ewentworth/cotney/genome/mm10/mm10.chrom.sizes /tgc/TGCore_User_Data/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/strong_novel_embryonic_tooth_enhancers.bigBed
chmod -R 755 /tgc/TGCore_User_Data/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/strong_novel_embryonic_tooth_enhancers.bigBed

bedToBigBed -type=bed3+1 hg19-25pct-strong_novel_embryonic_tooth_enhancers.bed /home/FCAM/ewentworth/cotney/genome/hg19/hg19.chrom.sizes /tgc/TGCore_User_Data/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/hg19-25pct-strong_novel_embryonic_tooth_enhancers.bigBed
chmod -R 755 /tgc/TGCore_User_Data/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/hg19-25pct-strong_novel_embryonic_tooth_enhancers.bigBed

for sample in $(cat samples.txt | egrep "fac|ooth|palate")
do
cat "$sample"* | bedtools sort -i stdin > $sample.dense.bed
bedToBigBed -type=bed9 $sample.dense.bed /home/FCAM/ewentworth/cotney/genome/hg19/hg19.chrom.sizes /tgc/TGCore_User_Data/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/$sample.bigBed
chmod -R 755 /tgc/TGCore_User_Data/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/$sample.bigBed
done


for sample in /tgc/TGCore_User_Data/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/hg19-25pct-embryonic*bigBed
do
export NAME=`echo $sample | sed 's/\/tgc\/TGCore_User_Data\/WebData\/cotney\/hubs\/ChIP\/ewentworth\/mousetooth\/hg19-25pct-//g; s/.bigBed//g'`
echo -e "track type=bigBed name='$NAME' description='$NAME' bigDataUrl=http://graveleylab.cam.uchc.edu/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/hg19-25pct-$NAME.bigBed itemRgb='On'"
done 




for sample in /tgc/TGCore_User_Data/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/*strong*ooth*bigBed 
do
export NAME=`echo $sample | sed 's/\/tgc\/TGCore_User_Data\/WebData\/cotney\/hubs\/ChIP\/ewentworth\/mousetooth\///g; s/.bigBed//g'`
echo -e "track type=bigBed name='$NAME' description='$NAME' bigDataUrl=http://graveleylab.cam.uchc.edu/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/$NAME.bigBed itemRgb='On'"
done 
