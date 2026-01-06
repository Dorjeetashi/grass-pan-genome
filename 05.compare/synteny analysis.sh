#Preparing CDS(coding DNA sequences) files
ln -s ../01.prepare_data/cds/Sorghum.fasta Sorghum.cds
ln -s ../01.prepare_data/cds/Maize.fasta Maize.cds
ln -s ../01.prepare_data/cds/Rice.fasta Rice.cds

#Preparing GFF(general feature format) files
ln -s ../01.prepare_data/Sorghum.longest_isoform.gff3 Sorghum.longest_isoform.gff3
ln -s ../01.prepare_data/Maize.longest_isoform.gff3 Maize.longest_isoform.gff3
ln -s ../01.prepare_data/Rice.longest_isoform.gff3 Rice.longest_isoform.gff3

#Generate gene location BED File
python3 -m jcvi.formats.gff bed --type=mRNA Sorghum.longest_isoform.gff3 -o Sorghum.bed
python3 -m jcvi.formats.gff bed --type=mRNA Maize.longest_isoform.gff3 -o Maize.bed
python3 -m jcvi.formats.gff bed --type=mRNA Rice.longest_isoform.gff3 -o Rice.bed

#Identification of synteny blocks between species
python3 -m jcvi.compara.catalog ortholog Sorghum Maize --no_strip_names --cscore 0.99
python3 -m jcvi.compara.catalog ortholog Maize Rice --no_strip_names --cscore 0.99

#Filtering and simplification of synteny blocks for visualization
python3 -m jcvi.compara.synteny screen --minspan=30 --simple Sorghum.Maize.anchors Sorghum.Maize.anchors.new
python3 -m jcvi.compara.synteny screen --minspan=30 --simple Maize.Rice.anchors Maize.Rice.anchors.new

#Obtain chromosome information
cat Sorghum.bed|cut -f1|sort -k1,1n|uniq|tr '\n' ','
cat Maize.bed|cut -f1|sort -k1,1n|uniq|tr '\n' ','
cat Rice.bed|cut -f1|sort -k1,1n|uniq|tr '\n' ','

#Preparing seqids file
echo "1,2,3,4,5,6,7,8,9,10
chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10
Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8,Chr9,Chr10,Chr11,Chr12" >grass.seqids

#Preparing layout file
echo "# y, xstart, xend, rotation, color, label, va,  bed
 .6,     .1,    .8,       0,      , Sorghum, top, Sorghum.bed
 .4,     .1,    .8,       0,      , Maize, bottom, Maize.bed
 .2,     .1,    .8,       0,      , Rice, bottom, Rice.bed
# edges
e, 0, 1, Sorghum.Maize.anchors.simple
e, 1, 2, Maize.Rice.anchors.simple" >grass.layout

#Synteny Visualization
python3 -m jcvi.graphics.karyotype  grass.seqids grass.layout -o grass.karyotype.pdf
