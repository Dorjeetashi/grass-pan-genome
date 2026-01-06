#Sorghum graph pangenome construction
cactus-pangenome ./js Sb_evolverPrimates.txt --outDir Sb_pan --outName Sb_pan --reference Sb_BTx623 --vcf --giraffe --gfa --gbz

#Maize graph pangenome construction
cactus-pangenome ./js Zm_evolverPrimates.txt --outDir Zm_pan --outName Zm_pan --reference Zm_B73 --vcf --giraffe --gfa --gbz

#Rice graph pangenome construction
cactus-pangenome ./js Os_evolverPrimates.txt --outDir Os_pan --outName Os_pan --reference Os_MSU --vcf --giraffe --gfa --gbz