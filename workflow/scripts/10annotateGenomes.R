### Annotate genome assemblies with prokka

bact <- c('BIGb248', 'BIGb306', 'BIGb446', 'BIGb468')

#The --proteins option is recommended when you have good quality reference genomes and want to ensure gene naming is consistent. Some species use specific terminology which will be often lost if you rely on the default Swiss-Prot database included with Prokka.

for (i in 1:length(bact)) {

  #command <- paste0('prokka --force --rfam --protein --outdir output/assembly/hgap4/annotation/prokka/', bact[1], ' input/assembly/HGAP4/fasta/', bact[1], '.fasta')
  command <- paste0('prokka --force --outdir output/assembly/hgap4/annotation/prokka/', bact[i], ' --prefix prokka input/assembly/HGAP4/fasta/', bact[i], '.fasta > prokka.log 2>&1')
  print(command)
  system(command)
}


#seqret -sequence input/assembly/HGAP4/fasta/BIGb248.fasta -feature -fformat gff -fopenfile output/assembly/hgap4/annotation/prokka/BIGb248/PROKKA_02282019.gff -osformat embl -auto -outseq output/assembly/hgap4/annotation/prokka/BIGb248/PROKKA_02282019.embl
#seqret -sequence input/assembly/HGAP4/fasta/BIGb306.fasta -feature -fformat gff -fopenfile output/assembly/hgap4/annotation/prokka/BIGb306/PROKKA_02282019.gff -osformat embl -auto -outseq output/assembly/hgap4/annotation/prokka/BIGb306/PROKKA_02282019.embl
#seqret -sequence input/assembly/HGAP4/fasta/BIGb446.fasta -feature -fformat gff -fopenfile output/assembly/hgap4/annotation/prokka/BIGb446/PROKKA_02282019.gff -osformat embl -auto -outseq output/assembly/hgap4/annotation/prokka/BIGb446/PROKKA_02282019.embl
#seqret -sequence input/assembly/HGAP4/fasta/BIGb468.fasta -feature -fformat gff -fopenfile output/assembly/hgap4/annotation/prokka/BIGb468/PROKKA_02282019.gff -osformat embl -auto -outseq output/assembly/hgap4/annotation/prokka/BIGb468/PROKKA_02282019.embl
