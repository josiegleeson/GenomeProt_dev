# functions

# import proteomics file (fragpipe, maxquant etc)
import_proteomics_data <- function(proteomics_file) {
  
  # required columns
  # Peptide, Protein, Mapped Proteins
  # optional: Protein Start
  
  # if its metamorpheus:
  # Peptide column is called 'Base Sequence', proteins matched are 'Protein Groups'
  # # Use str_replace to capture groups and replace with a delimiter
  # # delimits based on TA= (any num chars) followed by a |
  # split_data <- str_replace_all(data, "(TA=[A-Za-z0-9.]+)\\|", "\\1_SPLIT_")
  # 
  # # Split based on the new delimiter
  # split_data <- str_split(split_data, "_SPLIT_", simplify = TRUE)
  # 
  # # Convert to a data frame
  # df <- as.data.frame(split_data, stringsAsFactors = FALSE)
  
  # import proteomics
  protfile <- fread(proteomics_file)
  
  # ensure compatibility
  if (c("Peptide") %in% colnames(protfile)) {
    # skip
  } else if (c("Peptide Sequence") %in% colnames(protfile)) {
    protfile$Peptide <- protfile$`Peptide Sequence`
    protfile$`Peptide Sequence` <- NULL
  } else if (c("Stripped.Sequence") %in% colnames(protfile)) {
    protfile$Peptide <- protfile$Stripped.Sequence
    protfile$Stripped.Sequence <- NULL
  }
  
  
  # check if Protein Start is a valid column
  if (c("Protein Start") %in% colnames(protfile)) {
    protfile <- protfile %>% dplyr::select(Peptide, Protein, `Protein Start`, `Mapped Proteins`)
    # add a column of every mapped ORF
    protfile <- protfile %>% dplyr::mutate(
      all_mappings = case_when(
        (is.na(`Mapped Proteins`) == TRUE | nchar(`Mapped Proteins`) == 0) ~ paste0(Protein),
        TRUE ~ paste0(Protein, ", ", `Mapped Proteins`)
      )
    )
  } else if (c("All Mapped Proteins") %in% colnames(protfile)) {
    protfile <- protfile %>% dplyr::select(Peptide, `All Mapped Proteins`)
    protfile$all_mappings <- protfile$`All Mapped Proteins`
    protfile$`All Mapped Proteins` <- NULL
  }  else {
    protfile <- protfile %>% dplyr::select(Peptide, Protein, `Mapped Proteins`)
    # add a column of every mapped ORF
    protfile <- protfile %>% dplyr::mutate(
      all_mappings = case_when(
        (is.na(`Mapped Proteins`) == TRUE | nchar(`Mapped Proteins`) == 0) ~ paste0(Protein),
        TRUE ~ paste0(Protein, ", ", `Mapped Proteins`)
      )
    )
  }
  
  # separate into one row per mapped ORF
  prot_expanded <- separate_rows(protfile, all_mappings, sep = "\\,(?!chr)|\\; |\\;|\\, ")
  
  #prot_expanded <- prot_expanded[!grepl("\\,chr", prot_expanded$all_mappings),]
  
  # remove white space and characters
  prot_expanded <- prot_expanded %>% 
    dplyr::mutate(across(where(is.character), ~ gsub(" ", "", .))) %>%
    dplyr::mutate(PID = all_mappings)
  
  # remove any duplications
  prot_expanded <- prot_expanded[!(base::duplicated(prot_expanded)),]
  
  # select columns to keep
  if (c("Protein Start") %in% colnames(prot_expanded)) {
    prot_expanded <- prot_expanded %>% dplyr::select(PID, Peptide, `Protein Start`)
  } else {
    prot_expanded <- prot_expanded %>% dplyr::select(PID, Peptide)
  }
  
  # remove peptides mapped to reversed sequences
  prot_expanded <- prot_expanded %>% 
    dplyr::filter(!startsWith(PID, "rev"))
  
  # rename column
  prot_expanded$peptide <- prot_expanded$Peptide
  prot_expanded$Peptide <- NULL
  
  return(prot_expanded)
  
}


# import custom database metadata of ORFs
import_orf_metadata <- function(metadata_file) {
  
  orf_data <- fread(metadata_file)
  orf_data$transcript_id <- orf_data$transcript
  orf_data$transcript <- NULL
  orf_data$PID <- paste0(orf_data$accession, "|CO=", orf_data$orf_genomic_coordinates)
  orf_data <- orf_data %>% dplyr::select(PID, accession, orf_genomic_coordinates, transcript_id, transcript_biotype, orf_type, localisation, uniprot_status, openprot_id, 
                                         `molecular_weight(kDA)`, isoelectric_point, hydrophobicity, aliphatic_index, longest_orf_in_transcript)
  
  return(orf_data)
  
}

# import custom database FASTA of ORFs
import_fasta <- function(fasta_file, proteomics_data, gtf_file) {
  
  # fasta_file <- fasta_import_file
  # gtf_file <- gtf_import_file
  # proteomics_data <- pd
  
  # get transcript lengths
  gtf_txdb <- makeTxDbFromGFF(gtf_file)
  tx_lengths <- transcriptLengths(gtf_txdb)
  tx_lengths$transcript_id <- tx_lengths$tx_name
  tx_lengths <- tx_lengths %>% dplyr::select(transcript_id, tx_len)
  
  # get gene names
  gtf_import <- rtracklayer::import(gtf_file, format="gtf") %>% 
    as_tibble() %>% 
    dplyr::filter(type == "transcript") %>% 
    dplyr::select(transcript_id, gene_name, strand)
  
  # merge tx info
  tx_data <- merge(tx_lengths, gtf_import, by="transcript_id", all.x=T, all.y=F)
  
  # get transcripts for mapping
  txs <- exonsBy(gtf_txdb, by=c("tx"), use.names=T)
  
  # import fasta
  db <- readAAStringSet(fasta_file, format="fasta", use.names=TRUE)
  
  # convert to df
  convert_AA_to_df <- function(AAstring) {
    data.frame(name=names(AAstring),
               protein_sequence=as.character(AAstring),
               protein_length=as.numeric(nchar(as.character(AAstring))),
               row.names=NULL)
  }
  
  # convert
  df <- convert_AA_to_df(db)
  
  # get nt length of AA protein length
  df$nt_length <- df$protein_length * 3
  
  # format headers into columns
  suppressWarnings({
    df <- df %>%
      separate(name, into = c("PID"), sep = "\\ ", remove = FALSE)
    
    # filter for proteins with mapped peptides
    df <- df %>% dplyr::filter(PID %in% proteomics_data$PID)
    
    # filter to remove ORFs with multiple genomic loci
    # area we need to change
    df <- df[!grepl("\\,chr", df$PID),]
    
    df <- df %>% 
      separate(name, into = c("header", "gene_id", "gene_name", "transcript_id"), sep = "\\ ", remove = FALSE) %>%
      separate(gene_id, into = c("gene_id"), sep = "\\,", remove = FALSE) %>%
      separate(PID, into = c("protein_name", "location"), sep = "\\|CO=", remove = FALSE) %>%
      separate(location, into = c("chromosome", "start", "end"), sep = "\\:|-", remove = FALSE)
  })
  
  df$gene_name <- NULL
  
  df <- df %>%
    dplyr::mutate(across(where(is.character), ~ gsub("GA=", "", .))) %>%
    dplyr::mutate(across(where(is.character), ~ gsub("GN=", "", .))) %>%
    dplyr::mutate(across(where(is.character), ~ gsub("TA=", "", .)))
  
  # add one row per transcript
  df_expanded <- separate_rows(df, transcript_id, sep = "\\,|\\, ")
  df_expanded <- df_expanded[!(base::duplicated(df_expanded)),]
  
  df_expanded <- merge(df_expanded, tx_data, by="transcript_id", all.x=T, all.y=F)
  
  df_expanded$start <- as.numeric(df_expanded$start)
  df_expanded$end <- as.numeric(df_expanded$end)
  
  # get df of orf genomic coords
  orf_genomic_coords_df <- df_expanded
  
  # re orient start and end to match strand for mapping
  orf_genomic_coords_df <- orf_genomic_coords_df %>% 
    rowwise() %>% 
    mutate(stranded_start = case_when(
      strand == "+" ~ min(start,end),
      strand == "-" ~ max(start,end),
    ))
  
  # filter txdb for only those in df
  txs_filtered <- txs[names(txs) %in% orf_genomic_coords_df$transcript_id]
  txs_unlisted <- unlist(txs_filtered)
  
  # function to convert genomic start position to transcript position
  convert_gene_pos_to_transcript_pos <- function(txdb, input_df) {
    
    # get all relevant exons
    all_exons <- txdb[names(txdb) %in% input_df$transcript_id]
    
    # all exons info df
    exon_data <- data.frame(
      transcript_id = names(all_exons),
      exon_start = start(all_exons),
      exon_end = end(all_exons),
      exon_rank = mcols(all_exons)$exon_rank
    ) %>%
      group_by(transcript_id) %>%
      arrange(exon_rank) %>%
      mutate(
        exon_length = exon_end - exon_start + 1,
        total_length = cumsum(exon_length)
      )
    
    # join input data with exon data
    result <- input_df %>%
      left_join(exon_data, by = "transcript_id", relationship = "many-to-many") %>%
      dplyr::filter(stranded_start >= exon_start & stranded_start <= exon_end) %>%
      group_by(transcript_id, stranded_start) %>%
      slice(1) %>% 
      ungroup() %>%
      mutate(
        txstart = case_when(
          strand == "+" & exon_rank == 1 ~ stranded_start - exon_start,
          strand == "+" ~ (stranded_start - exon_start) + (total_length - exon_length),
          strand == "-" & exon_rank == 1 ~ exon_end - stranded_start,
          strand == "-" ~ (exon_end - stranded_start) + (total_length - exon_length)
        )
      )
    
    return(result)
  }
  
  # apply to orf genomic coords
  orf_transcriptomic_coords <- convert_gene_pos_to_transcript_pos(txs_unlisted, orf_genomic_coords_df)
  
  # ORFs that could not be mapped to transcripts are returned with -1 starts
  orf_transcriptomic_coords <- orf_transcriptomic_coords %>% dplyr::filter(txstart >= 0)
  
  # get transcript end coord
  orf_transcriptomic_coords$txend <- orf_transcriptomic_coords$txstart + orf_transcriptomic_coords$nt_length - 1
  
  # remove ORF coords outside transcript ends
  orf_transcriptomic_coords <- orf_transcriptomic_coords %>% dplyr::filter(txend < tx_len)
  
  # combine with proteomics data
  metadata <- full_join(orf_transcriptomic_coords, proteomics_data, by = "PID")
  metadata <- metadata[!(base::duplicated(metadata)),]
  metadata <- metadata %>% dplyr::filter(!is.na(transcript_id) & !is.na(txstart))
  
  find_peptide_position <- function(peptides, protein_sequences) {
    
    # exact match
    exact_matches <- stri_locate_first_fixed(protein_sequences, peptides)
    
    # non-exact matches use fuzzy matching
    non_exact_mask <- is.na(exact_matches[,1])
    
    if (any(non_exact_mask)) {
      fuzzy_matches <- vapply(which(non_exact_mask), function(i) {
        peptide <- peptides[i]
        protein <- protein_sequences[i]
        
        # generate regex pattern allowing up to 3 mismatches
        pattern <- paste0(strsplit(peptide, "")[[1]], collapse = "(.{0,1})")
        
        match <- stri_locate_first_regex(protein, pattern)
        if (is.na(match[1])) -1L else as.integer(match[1])
      }, integer(1))
      
      exact_matches[non_exact_mask,1] <- fuzzy_matches
    }
    return(as.integer(exact_matches[,1]))
  }
  
  # apply to dataframe
  setDT(metadata)
  metadata[, mapped_pep_start := find_peptide_position(peptide, protein_sequence)]
  
  # if exact match or one mismatch isn't found, use proteomics defined start site
  if (c("Protein Start") %in% colnames(metadata)) {
    metadata <- metadata %>% mutate(pep_start = case_when(
      (mapped_pep_start != `Protein Start` & mapped_pep_start != -1) ~ mapped_pep_start,
      (mapped_pep_start != `Protein Start` & mapped_pep_start == -1) ~ `Protein Start`,
      TRUE ~ mapped_pep_start
    ))
  } else {
    metadata <- metadata %>% dplyr::filter(mapped_pep_start != -1 & !is.na(mapped_pep_start))
    metadata$pep_start <- metadata$mapped_pep_start
  }
  
  # get peptide end location within every ORF
  metadata$pep_end <- metadata$pep_start + nchar(metadata$peptide) - 1
  
  metadata_output <- metadata %>% dplyr::filter(pep_start < protein_length & pep_end <= protein_length)
  
  metadata_output$mapped_pep_start <- NULL
  metadata_output$name <- NULL
  metadata_output$header <- NULL
  
  return(metadata_output)
  
}

# get peptide coords from transcript coords
extract_peptide_coords <- function(metadata_df) {
  
  # subset metadata for conversion of peptide coords
  txcoordsdf_subset <- metadata_df %>% dplyr::select(orf_tx_id, transcript_id, PID, gene_id, pep_start, pep_end, peptide, txstart, txend)
  
  # get peptide locations within every transcript
  txcoordsdf_subset$txstart <- as.numeric(txcoordsdf_subset$txstart)
  txcoordsdf_subset$txend <- as.numeric(txcoordsdf_subset$txend)
  
  txcoordsdf_subset$start <- txcoordsdf_subset$txstart + ((txcoordsdf_subset$pep_start-1) * 3)
  txcoordsdf_subset$end <- txcoordsdf_subset$start + (nchar(txcoordsdf_subset$peptide) * 3)
  #txcoordsdf_subset$width <- txcoordsdf_subset$end - txcoordsdf_subset$start
  txcoordsdf_subset$peptide <- as.factor(txcoordsdf_subset$peptide)
  
  txcoordsdf_subset <- txcoordsdf_subset %>% 
    dplyr::filter(start >= txstart & end <= txend)
  
  # GRanges from df of transcriptomic coords of peptides
  coords_peptides <- makeGRangesFromDataFrame(txcoordsdf_subset,
                                              keep.extra.columns=TRUE,
                                              ignore.strand=FALSE,
                                              seqinfo=NULL,
                                              seqnames.field="transcript_id",
                                              start.field="start",
                                              end.field="end",
                                              strand.field="strand",
                                              starts.in.df.are.0based=TRUE,
                                              na.rm=TRUE)
  # set names
  names(coords_peptides) <- c(coords_peptides$transcript)
  
  # filter peptides for bad mappings
  coords_peptides <- subset(coords_peptides, start(coords_peptides) != 0)
  coords_peptides <- subset(coords_peptides, end(coords_peptides) <= coords_peptides$txend)
  
  return(coords_peptides)
  
}


