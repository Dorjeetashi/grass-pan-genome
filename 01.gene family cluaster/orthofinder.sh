#Reference genome–based gene clustering of sorghum, maize, and rice
orthofinder  -f ref_pep \
 -S diamond \
 -M msa \
 -A muscle \
 -T fasttree \
 -t 40 \
 -a 2 
 #Pangenome-based gene clustering of sorghum, maize, and rice
orthofinder  -f pan_pep \
 -S diamond \
 -M msa \
 -A muscle \
 -T fasttree \
 -t 40 \
 -a 2 