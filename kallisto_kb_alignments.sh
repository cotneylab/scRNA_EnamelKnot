module load kallisto
kb ref --overwrite -i $GENOME/mm10_kallisto_index.idx -g t2g.txt -f1 $GENOME/cDNA.fasta $GENOME/Mus_musculus.GRCm39.dna.primary_assembly.fa $GENOME/kallisto.gencode.vM25.annotation.gtf

##############bulk RNA-seq alignment. renamed files based on the replicates.
for sample in e*molar/
do
export NAME=`echo $sample | sed 's/\///g'`
for i in {1..7}
do
echo -e "kallisto quant -i $GENOME/mus_musculus/transcriptome.idx -o $sample/fastqs/"$NAME"_"$i" --single -l 180 -s 20 --bootstrap-samples=10 $sample/fastqs/"$NAME"_"$i".fastq.gz >> ~/cotney/rawdata/GEO/mouse/embryonictoot/bulkRNAseq/initial_analysis.sh
done
done

cd /home/FCAM/ewentworth/cotney/rawdata/GEO/mouse/embryonictoot/bulkRNAseq/
ls -d e* | egrep -v ".fastq" > temp.txt

###make s2c.txt file for importing output files into DESeq using tximport.
for file in $(cat temp.txt); do 
export NAME=`echo $file | sed 's/_molar_/_molar\t/g' | cut -f1`
echo -e $file $NAME /home/FCAM/ewentworth/cotney/rawdata/GEO/mouse/embryonictoot/bulkRNAseq/$file >> s2c.txt; 
done

for sample in $(awk '{print $3}' s2c.txt)
do
kallisto h5dump -o $sample $sample/abundance.h5
done


#########################single cell RNA-seq alignment
##download files from geo
module load sratoolkit/2.10.8
egrep "SRR" SRA_runs_table.tsv | sed 's/,/\t/g' | cut -f1 | awk 'BEGIN{ORS=" "}1' | awk '{print "fastq-dump --split-files -o fastqs", $0}'

module load kallisto
cd fastqs/
for sample in e*R1*fastq.gz
do
export name=`echo $sample | sed 's/_R1.*//g'`
export R2=`ls $name*R2*`
echo $name
echo $R2
echo $sample
mkdir $name-rna-out
kb count -i $GENOME/transcriptome.idx \
-g $GENOME/transcripts_to_genes.txt -x 10XV3 -o $name-rna-out --overwrite \
$sample $R2
cd /home/CAM/ewentworth/cotney/rawdata/GEO/mouse/embryonictoot/scRNAseq/combined/fastqs
done
