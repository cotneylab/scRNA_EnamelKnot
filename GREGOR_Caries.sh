#!/bin/bash
#SBATCH --job-name=cariesgregor
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mail-type=END
#SBATCH --mem=64G
#SBATCH --mail-user=wentworth@uchc.edu
#SBATCH -o %j.out
#SBATCH -e %j.err
source ~/.bash_profile
export DIR=/home/FCAM/ewentworth/cotney/analysis/gwas/gregor/caries

cd $DIR


export DIR1=/home/FCAM/ewentworth/cotney/rawdata/chromatinsegmentations/superenhancers/human
export DIR2=/home/FCAM/ewentworth/cotney/rawdata/chromatinsegmentations/human_25state
ls $DIR2/gregor-*.bed > index_human_25state.txt
#ls $DIR1/*.bed >> index_human_25state.txt


egrep -hi "gregor-13|gregor-14|gregor-18|gregor-15|enhancers" index_human_25state.txt > index_human_activeenhancers_25state.txt

export DIR3=/home/FCAM/ewentworth/cotney/rawdata/chromatinsegmentations/human_atac

ls $DIR3/*enhancers.bed > index_human_atac_enhanceroverlap.txt
ls $DIR3/*hg19* | egrep -v 'enhancers' > index_human_atac.txt

export DIR4=/home/FCAM/ewentworth/cotney/rawdata/chromatinsegmentations/human_18state
ls $DIR4/*.bed > index_human_18state.txt
egrep 'state9|state10|state8|merge|activeenhancers.bed' index_human_18state.txt > index_human_activeenhancers_18state.txt


export DIR5=/home/FCAM/ewentworth/cotney/rawdata/chromatinsegmentations/mouse_18state
export DIR6=/home/FCAM/ewentworth/cotney/rawdata/chromatinsegmentations/mouse_18state/tissue_specific_enhancer_analysis
ls $DIR5/hg19*.bed | egrep -v "consensus" > index_mouse.txt
egrep -v "-e11.5-wt_pas_" index_mouse.txt > temp.txt
mv temp.txt index_mouse.txt
egrep 'state9|state10|state8|merge|strong_novel_embryonic_toot' index_mouse.txt > index_mouse_activeenhancers.txt


############dental phenotypeaddgregateSNPS.txt is a cat list of dentaldisease SNPs from GWAS catalog (mostly caries) and odontogenesis SNPs from GWAS catalog

for sample in index_*.txt
do
export NAME=`echo $sample | sed 's/index_//g' | sed 's/.txt//g'`
echo $sample
for filename in *SNPs.txt
do
export filename2=`echo $filename | sed 's/SNPs.txt//g'`
mkdir $DIR/$NAME-$filename2-out
echo -e "INDEX_SNP_FILE = $DIR/$filename\nBED_FILE_INDEX = $DIR/$sample\nREF_DIR = /home/FCAM/ewentworth/cotney/tools/gregor/refdir\nR2THRESHOLD = 0.8\nLDWINDOWSIZE = 1000000\nOUT_DIR = $DIR/$NAME-$filename2-out\nMIN_NEIGHBOR_NUM = 500\nBEDFILE_IS_SORTED = True\nPOPULATION = EUR\nTOPNBEDFILES = 2\nJOBNUMBER = 10\n###############################################################################\nBATCHTYPE = slurm\nBATCHOPTS = --partition=general --qos=general -n 2 -c 16 --mem=150G " > $DIR/$NAME-$filename2-conf.file
echo -e "#/bin/bash\n#SBATCH --job-name=caries\n#SBATCH -N 1\n#SBATCH -n 2\n#SBATCH -c 16\n#SBATCH --partition=general\n#SBATCH --qos=general\n#SBATCH --mail-type=END\n#SBATCH --mem=150G\n#SBATCH --mail-user=wentworth@uchc.edu\n#SBATCH -o %j.out\n#SBATCH -e %j.err\ncd $DIR\nperl /home/FCAM/ewentworth/cotney/tools/gregor/GREGOR/script/GREGOR.pl --conf $DIR/$NAME-$filename2-conf.file" | sed 's/\/bin/!\/bin/g' > $DIR/gregor-$NAME-$filename2.sh
#sbatch $DIR/gregor-$NAME-$filename2.sh
done
done
####make sure to use the asian population setting for the yu clefting analysis
for sample in *yu*conf.file
do
sed -i 's/POPULATION = EUR/POPULATION = ASN/g' $sample
done

for sample in *active*18*sh *mouse*active*sh
do
sbatch $sample
done


for sample in */
do
export NAME=`echo $sample | sed 's/-out\///g'`
export STAT=`echo $sample | sed 's/$/StatisticSummaryFile.txt/g'`
export DISEASE=`echo $sample | sed 's/mouse//g; s/human//g; s/_activeenhancers//g; s/_atac-//g; s/_enhanceroverlap-//g; s/_18state-//g; s/_25state-//g; s/-//g; s/out//g; s/\///g; s/_clefts/Clefting/g; s/yu/Yu /g; s/ludwig_2017/Ludwig/g; s/craniofacial_measures/Craniofacial Variation/g; s/_atacc/C/g; s/caries/Caries/g; s/crohns/Crohns Disease/g'`
export ASSAY=`echo $NAME | awk '{if ($1 ~ /human_atac_enhanceroverlap/) {print "Human Enhancer Overlap ATAC Peaks"} else if ($1 ~ /human_atac/) {print "Human ATAC Peaks"} else if ($1 ~ /human_activeenh/) {print "Human Active Enhancers"} else if ($1 ~ /mouse_active/) {print "Mouse 18 State Active Enhancers"} else if ($1 ~ /mouse/) {print "Mouse 18 State"} else if ($1 ~ /human_18state/) {print "Human 18 State"} else {print "Human 25 State"}}'`
awk 'NR==1{print;next}{if ($3=="NA") {print $0, 0} else {print $0, (log($2/$3)/log(2))}}' $STAT | tail -n +2 | sed 's/dense.sorted.bed//g; s/_25_imputed12marks_dense.enhancer_states.bed//g; s/reproducible_enhancer_states.25.bed/All_CS_Enhancers/g; s/reproducible_craniofacial_specific_enhancer_states.bed/CS_Specific_Enhancers/g; s/hg19-25pct-//g; s/-sorted//g; s/gregor-//g; s/.bed//g; s/_18_core_K27ac_segments//g; s/.enhancers//g; s/.hg19.narrowPeak//g; s/impute_//g' | awk '{if ($1 ~ /ENCFF/) {print $0, "misc"} else if ($1 ~ /CS|F2|face|cranio|lower|upper|ooth|NCC|snout|faci|CF|cs|Late|Early|Middle|Cells|Mesen|derm|Neural|Sox|Muscle|Endo|embryonic_9.5_day_embryonic|cfse/) {print $0, "craniofacial"} else if ($1 ~ /ENCFF/) {print $0, "misc"} else if ($1 ~ /eart|cardi/) {print $0, "heart"} else {print $0, "misc"}}' | awk '{if ($1 ~ /E080|E081|E082|E083|E084|E085|E086|E087|E088|E089|fetal/) {print $0, "fetal"} else if ($1 ~ /CS|F2|E001|E002|E003|E004|E005|E006|E007|E008|E009|E010|E011|E012|E013|E014|E015|E016|E024|E018|E019|E020|E021|E022|eart|embryonic/) {print $0, "embryonic"} else if ($1 ~ /ENCFF/) {print $0, "adult"} else {print $0, "adult"}}' | awk '{if ($1 ~ /merged_sorted_all_strong/) {print $0, "mouse"} else {print $0, "human"}}' | sed 's/ /\t/g' | egrep -v "e11.5-wt_pas_mesenchyme_markerpeaks|SKNSH|K562|forebrain|diff|A549|Blood|Sox|snout|palate" | sed 's/cs/CS/g; s/Cells/Cells_Pseudobulk/g' | egrep -v 'ncc|unknown|combined', > $DIR/$NAME-out.txt
done

#Combine output files for human and mouse analysis

cat human_activeenhancers_18state-caries-out.txt mouse_activeenhancers-caries-out.txt > combined_activeenhancers-caries-out.txt
cat human_activeenhancers_18state-dentaldisease-out.txt mouse_activeenhancers-dentaldisease-out.txt > combined_activeenhancers-dentaldisease-out.txt
cat human_activeenhancers_18state-odontogenesis-out.txt mouse_activeenhancers-odontogenesis-out.txt > combined_activeenhancers-odontogenesis-out.txt
cat human_activeenhancers_18state-crohns-out.txt mouse_activeenhancers-crohns-out.txt > combined_activeenhancers-crohns-out.txt
cat human_activeenhancers_18state-combined-out.txt mouse_activeenhancers-combined-out.txt > combined_activeenhancers-combined-out.txt
cat human_activeenhancers_18state-lupus-out.txt mouse_activeenhancers-lupus-out.txt > combined_activeenhancers-lupus-out.txt
cat human_activeenhancers_18state-yu_cleft-out.txt mouse_activeenhancers-yu_cleft-out.txt > combined_activeenhancers-yu_cleft-out.txt
cat human_activeenhancers_18state-ruiz_facialmorphology-out.txt mouse_activeenhancers-ruiz_facialmorphology-out.txt > combined_activeenhancers-ruiz_facialmorphology-out.txt
cat human_activeenhancers_18state-ludwig_cleft-out.txt mouse_activeenhancers-ludwig_cleft-out.txt > combined_activeenhancers-ludwig_cleft-out.txt
cat human_activeenhancers_18state-dentalphenotypeaggregate-out.txt mouse_activeenhancers-dentalphenotypeaggregate-out.txt > combined_activeenhancers-dentalphenotypeaggregate-out.txt

####Generate some figures and correct those p values using BH!

for sample in combined*out.txt
do
export NAME=`echo $sample | sed 's/-out.txt//g'`
export DISEASE=`echo $sample | awk '{if ($0~ /ludwig/) {print "Ludwig Clefting"} else if ($0~ /yu/) {print "Yu Clefting"} else if ($0~ /craniofacial/) {print "Craniofacial Metrics"} else if ($0~ /caries/) {print "Caries"} else if ($0~ /crohn/) {print "Crohns Disease"}}'`
#export DISEASE=`echo $sample | sed 's/-out.txt//g; s/combined_activenhancers_18state-//g; s/_/ /g;'`
echo -e "library(ggplot2)\nlibrary(ggrepel)\nlibrary(stringr)\nlibrary(writexl)\nlibrary(patchwork)\nx<-read.table(\""$DIR/$NAME-out.txt"\")\nname<-gsub(\"-out.txt\", \"\", \""$DIR/$NAME-out.txt"\")\nx[is.na(x)]<-1\nx\$V9<-gsub(\"state\", \"state \", gsub(\"activ\", \"merged\", gsub(\"_|-\", \" \", gsub(\"hg19-25pct-|_rep1|_rep2|_rep3|_rep4|.bed|_sorted_all_strong|-sorted|_18_core_K27ac_segments\", \"\", gsub(\"strong_novel_embryonic_tooth\", \"TSEs\", x\$V1)))))\nw<-cbind(x, (p.adjust(x\$V4, \"BH\")), (p.adjust(x\$V4, \"bonferroni\")), -log10(p.adjust(x\$V4, \"BH\")), -log10(p.adjust(x\$V4, \"bonferroni\")))\nx<-w\ncolnames(x)<-c(\"Sample\", \"ObservedSNPs\", \"ExpectedSNPs\", \"PValue\", \"Log2SNPFoldEnrichment\", \"Lineage\", \"DevelopmentalStage\", \"Organism\", \"Label\", \"BHCorrectedPValue\", \"BonferroniCorrectedPValue\", \"BHCorrectedPValueLog\", \"BonferroniCorrectedPValueLog\")\nx<-x[!duplicated(x),]\nn<-head(x[order(-x\$Log2SNPFoldEnrichment),], 35)\nq<-head(x[order(-x\$BHCorrectedPValueLog),], 35)\nx\$State<-gsub(\".*merged\", \"merged\", gsub(\".*state\", \"state\", x\$Label))\nminy<-0.1\nmaxy<-range(x\$BHCorrectedPValueLog)[2]+1\nminx<-range(x\$Log2SNPFoldEnrichment)[1]-0.25\nmaxx<-range(x\$Log2SNPFoldEnrichment)[2]+0.25\ns<-ggplot(x, aes(x=Log2SNPFoldEnrichment, y=BHCorrectedPValueLog, colour=State, shape=Lineage, label=Label)) + geom_point(aes(size=ObservedSNPs))+ geom_hline(yintercept=1.3) + ggtitle(\""$DISEASE"\") + scale_color_manual(values=c(\"grey\", \"darkorange\", \"greenyellow\", \"orange\", \"darkgreen\")) + labs(x=\"Fold Enrichment (Log2)\", y=\"BH Corrected P Value (-Log10)\", colour=\"State\", size=\"Number of SNPs per Sample\", shape=\"Lineage\") + xlim(minx, maxx) + ylim(miny,maxy)\nt<-ggplot(x, aes(x=Log2SNPFoldEnrichment, y=BHCorrectedPValueLog, colour=Lineage, shape=Organism, label=Label)) + geom_point(aes(size=ObservedSNPs))+ geom_hline(yintercept=1.3) + ggtitle(\""$DISEASE"\") + scale_color_manual(values=c(\"orange\", \"grey\", \"grey\")) + labs(x=\"Fold Enrichment (Log2)\", y=\"BH Corrected P Value (-Log10)\", colour=\"Lineage\", size=\"Number of SNPs per Sample\", shape=\"Organism\") + xlim(minx, maxx) + ylim(miny, maxy)\npdf(file=paste(\""$NAME.pdf"\"), height=8.5, width=11)\ns + geom_label_repel(aes(Log2SNPFoldEnrichment, BHCorrectedPValueLog, label=Label), xlim=c(minx, maxx), ylim=c(miny,maxy), data=subset(x, (Lineage == \"craniofacial\" & BHCorrectedPValueLog >= q[15,11] & Log2SNPFoldEnrichment >= n[15,5]) | (Log2SNPFoldEnrichment >= n[15,5] & BHCorrectedPValueLog >= q[15,11]) | Label == 'E13.5 Tooth'))\nt + geom_label_repel(aes(Log2SNPFoldEnrichment, BHCorrectedPValueLog, label=Label), xlim=c(minx, maxx), ylim=c(miny,maxy), data=subset(x, (Lineage == \"craniofacial\" & BHCorrectedPValueLog >= q[15,11] & Log2SNPFoldEnrichment >= n[15,5]) | (Log2SNPFoldEnrichment >= n[15,5] & BHCorrectedPValueLog >= q[15,11]) | Label == \"E13.5 Tooth\"))\ndev.off()\n" > $DIR/$NAME-figures.Rscript
done

for sample in *Rscript
do
echo $sample
Rscript $sample
done


