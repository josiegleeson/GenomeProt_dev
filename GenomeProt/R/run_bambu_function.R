
# edit to take in output path so we can supply session ID

run_bambu_function <- function(bam_file_list, gtf, organism, output_directory) { 
  
  # set organism
  if (organism == "human") {
    library(BSgenome.Hsapiens.UCSC.hg38)
    genomedb <- BSgenome.Hsapiens.UCSC.hg38
  } else if (organism == "mouse") {
    library(BSgenome.Mmusculus.UCSC.mm39)
    genomedb <- BSgenome.Mmusculus.UCSC.mm39
  } else if (organism == "celegans") {
    library(BSgenome.Celegans.UCSC.ce11)
    genomedb <- BSgenome.Celegans.UCSC.ce11
  } else if (organism == "drosophila") {
    library(BSgenome.Dmelanogaster.UCSC.dm6)
    genomedb <- BSgenome.Dmelanogaster.UCSC.dm6
  } else if (organism == "rat") {
    library(BSgenome.Rnorvegicus.UCSC.rn7)
    genomedb <- BSgenome.Rnorvegicus.UCSC.rn7
  } else if (organism == "zebrafish") {
    library(BSgenome.Drerio.UCSC.danRer11)
    genomedb <- BSgenome.Drerio.UCSC.danRer11
  }
  
  bambuAnnotations <- prepareAnnotations(gtf)
  se <- bambu(reads = bam_file_list, annotations = bambuAnnotations, genome = genomedb)
  writeBambuOutput(se, path = output_directory)
  
  tx_data <- as.data.frame(mcols(se))
  tx_data <- as.data.frame(apply(tx_data, 2, as.character))
  tx_data <- tx_data %>% dplyr::filter(novelTranscript == "TRUE")
  
  write.csv(tx_data, paste0(output_directory, "/novel_transcript_classes.csv"), row.names=F, quote=F)
  
}
