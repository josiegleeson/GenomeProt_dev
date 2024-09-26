
# inputs: MQ/FragPipe peptides.tsv output, GTF, metadata
# outputs: BED12/GTF files of peptides, ORFs and transcripts, database of peptides with info on locations etc, summary file of peptides

library(stringi)
source("global.R")
source("R/integration_functions.R")

option_list = list(
  make_option(c("-p", "--proteomics"), type="character", default=NULL,
              help="Proteomics data file", metavar="character"),
  make_option(c("-f", "--fasta"), type="character", default=NULL,
              help="Custom FASTA used for proteomics", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="Custom metadata used for proteomics", metavar="character"),
  make_option(c("-g", "--gtf"), type="character", default=NULL,
              help="GTF used to generate custom FASTA", metavar="character"),
  make_option(c("-s", "--savepath"), type="character", default=NULL,
              help="Output directory", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

options(scipen=999)

proteomics_import_file <- opt$proteomics
fasta_import_file <- opt$fasta
metadata_import_file <- opt$metadata
gtf_import_file <- opt$gtf
output_directory <- opt$savepath
 
# source("~/Documents/GenomeProt/GenomeProt/R/integration_functions.R")
# proteomics_import_file <- "~/Documents/miguel_data/august/diann/report.pr_matrix.tsv"
# fasta_import_file <- "~/Documents/miguel_data/august/database/proteome_database_filt.fasta"
# metadata_import_file <- "~/Documents/miguel_data/august/database/proteome_database_metadata_filt.txt"
# gtf_import_file <- "~/Documents/miguel_data/august/database/proteome_database_transcripts.gtf"
# output_directory <- "~/Documents"

# ------------- import files ------------- #

pd <- suppressWarnings(import_proteomics_data(proteomics_import_file))

gtf <- makeTxDbFromGFF(gtf_import_file) # make txdb of gtf

orf_df <- import_orf_metadata(metadata_import_file)

md <- import_fasta(fasta_import_file, pd, gtf_import_file)

# ---------------------------------------- #


# ------------- run analysis ------------- #

# extract ORF transcript coordinates to df
md$orf_tx_id <- paste0(md$protein_name, "_", md$transcript_id)

# get just unique orf and transcript for mapping
orf_transcript_coords_df <- md %>% dplyr::select(orf_tx_id, txstart, txend, transcript_id, gene_id, strand)
orf_transcript_coords_df <- orf_transcript_coords_df[!(base::duplicated(orf_transcript_coords_df)),]

# make GRanges from df of ORF transcript coordinates
orf_transcript_coords <- makeGRangesFromDataFrame(orf_transcript_coords_df,
                                                  keep.extra.columns=TRUE, ignore.strand=FALSE, seqinfo=NULL,
                                                  seqnames.field="transcript_id", start.field="txstart", end.field="txend", strand.field="strand",
                                                  starts.in.df.are.0based=FALSE, na.rm=TRUE)

names(orf_transcript_coords) <- c(orf_transcript_coords$orf_tx_id) # set names
mcols(orf_transcript_coords)$gene_id <- c(orf_transcript_coords_df$gene_id)

# get exons for mapping
exons <- exonsBy(gtf, "tx", use.names=TRUE) # get exon data per transcript
exons_filt <- exons[names(exons) %in% orf_transcript_coords_df$transcript_id] # filter for only transcripts with peptides

orf_tx_names <- as.character(seqnames(orf_transcript_coords)) # get tx names

# match names of transcripts, return index of match
names(orf_transcript_coords) <- match(orf_tx_names, names(exons_filt)) 

# original peptides in the GRanges
orf_ids <- orf_transcript_coords$orf_tx_id

# ORFik map to genome coordinates
orf_in_genomic <- ORFik::pmapFromTranscriptF(orf_transcript_coords, exons_filt, removeEmpty = T)

# map back to GRangesList, with group information
orf_in_genomic@unlistData$PID <- orf_ids[groupings(orf_in_genomic)]

# add exon_number for GTF export
orf_in_genomic_gr <- unlist(orf_in_genomic, use.names=F) # convert to GRanges

# create vector of exon number per peptide and transcript
exon_number_vec <- ave(seq_along(orf_in_genomic_gr), mcols(orf_in_genomic_gr)$PID, FUN = seq_along)

# add to GRanges
mcols(orf_in_genomic_gr)$exon_number <- exon_number_vec
# re-list
orf_in_genomic <- split(orf_in_genomic_gr, ~ mcols(orf_in_genomic_gr)$PID)

# peptides
# use ORF transcript coords to determine peptide transcript coords
peptide_transcript_coords <- extract_peptide_coords(md)

# map to genomic coords
peptide_tx_names <- as.character(seqnames(peptide_transcript_coords)) # get tx names

# match names of transcripts, return index of match
names(peptide_transcript_coords) <- match(peptide_tx_names, names(exons_filt)) 

# original peptides in the GRanges
pep_ids <- peptide_transcript_coords$peptide
pep_PID_ids <- peptide_transcript_coords$PID
pep_gene_ids <- peptide_transcript_coords$gene_id

# ORFik map to genome coordinates

# causes script to fail
pep_in_genomic <- ORFik::pmapFromTranscriptF(peptide_transcript_coords, exons_filt, removeEmpty = F)

# map back to GRangesList, with group information
pep_in_genomic@unlistData$peptide <- pep_ids[groupings(pep_in_genomic)]
pep_in_genomic@unlistData$PID <- pep_PID_ids[groupings(pep_in_genomic)]
pep_in_genomic@unlistData$gene <- pep_gene_ids[groupings(pep_in_genomic)]

# remove 0 ranges
pep_in_genomic_gr <- unlist(pep_in_genomic, use.names=F) # convert to GRanges
# rename with transcript and peptide
tx_pep_names <- c(paste0(names(pep_in_genomic_gr), "_", pep_in_genomic_gr$peptide))
names(pep_in_genomic_gr) <- tx_pep_names # set names
pep_in_genomic_gr$tx_pid_grouping <- paste0(pep_in_genomic_gr$PID, "_", names(pep_in_genomic_gr))
# remove 0 ranges
pep_in_genomic_gr <- subset(pep_in_genomic_gr, (start(pep_in_genomic_gr) != 0 & end(pep_in_genomic_gr) != 0)  )

# create vector of exon number per peptide and transcript
exon_number_vec <- ave(seq_along(pep_in_genomic_gr), pep_in_genomic_gr$tx_pid_grouping, FUN = seq_along)
# add to GRanges
mcols(pep_in_genomic_gr)$exon_number <- exon_number_vec
mcols(pep_in_genomic_gr)$tx_pid_grouping <- NULL

# re-list
pep_in_genomic <- split(pep_in_genomic_gr, ~ names(pep_in_genomic_gr))

# ---------------------------------------- #


# ------------- export ------------- #

# export bed12 of peptides
ORFik::export.bed12(pep_in_genomic, paste0(output_directory, "/peptides.bed12"), rgb = 0)

# export bed12 of ORFs
ORFik::export.bed12(orf_in_genomic, paste0(output_directory, "/ORFs.bed12"), rgb = 0)

# format GTF of ORFs
orf_in_genomic_gr$source <- c("custom")
orf_in_genomic_gr$type <- c("CDS")
orf_in_genomic_gr$phase <- 0
orf_in_genomic_gr$ORF_id <- names(orf_in_genomic_gr)
orf_in_genomic_gr$transcript_id <- names(orf_in_genomic_gr)
names(orf_in_genomic_gr) <- NULL
orf_in_genomic_gr$group_id <- "ORFs"

# format GTF of all transcripts that had mapped peptides
gtf_for_exporting <- import(gtf_import_file, format="gtf")
gtf_filtered <- gtf_for_exporting[mcols(gtf_for_exporting)$transcript_id %in% md$transcript]
gtf_filtered$group_id <- "transcripts"

# reformat exons for bed12
gtf_as_bed12 <- gtf_filtered[mcols(gtf_filtered)$type == "exon"]

names(gtf_as_bed12) <- paste0(gtf_as_bed12$transcript_id, "_", gtf_as_bed12$gene_id)

# convert to grl
tx_in_genomic <- split(gtf_as_bed12, ~ names(gtf_as_bed12))

# export bed12 of transcripts
ORFik::export.bed12(tx_in_genomic, paste0(output_directory, "/transcripts.bed12"), rgb = 0)

# ---------------------------------------- #


# ------- summary file -------- #

# convert to df
mcols(pep_in_genomic_gr)$txname <- names(pep_in_genomic_gr)
results_pept_df <- pep_in_genomic_gr %>% as_tibble()
results_pept_df <- separate(results_pept_df, txname, into = c("transcript_id", "peptide"), sep = "_", remove = TRUE)
results_pept_df$gene_id <- results_pept_df$gene
results_pept_df$gene <- NULL

# group by peptide and transcript to summarise based on how many exons peptide spans
results_pept_df_unique <- results_pept_df %>% 
  dplyr::group_by(peptide, transcript_id) %>% 
  dplyr::slice_max(exon_number) %>% 
  dplyr::mutate(number_exons = exon_number) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-start, -end, -width, -exon_number)

# merge results with metadata
metadata_to_merge <- md %>% 
  dplyr::select(PID, peptide, transcript_id, gene_id, gene_name, protein_length, tx_len)

peptide_result <- left_join(results_pept_df_unique, metadata_to_merge, by=c("PID", "peptide", "gene_id", "transcript_id"))

peptide_result <- peptide_result[!(base::duplicated(peptide_result)),]

peptide_result <- left_join(peptide_result, orf_df, by=c("PID", "transcript_id"))

peptide_result <- peptide_result[!(base::duplicated(peptide_result)),]

peptide_result$seqnames <- NULL

peptide_result <- peptide_result %>% 
  dplyr::mutate(longest_orf_in_transcript = case_when(
    longest_orf_in_transcript == "Y" ~ TRUE,
    longest_orf_in_transcript == "N" ~ FALSE
  ))

setDT(peptide_result)

peptide_result[, c("peptide_ids_gene", "peptide_ids_orf", "peptide_ids_transcript", "shared_novel_protein_peptide") := 
                 .(length(unique(gene_id)) == 1,
                   length(unique(gene_id)) == 1 & length(unique(PID)) == 1,
                   length(unique(gene_id)) == 1 & length(unique(PID)) == 1 & length(unique(transcript_id)) == 1,
                   length(unique(PID)) > 1 & all(startsWith(PID, "ORF"))),
               by = peptide]

peptide_result[, orf_identified := any(peptide_ids_orf == TRUE), by = PID]
peptide_result[, gene_identified := any(peptide_ids_gene == TRUE), by = gene_id]
peptide_result[, transcript_identified := any(peptide_ids_transcript == TRUE), by = transcript_id]

# get missing peptides back in the output
# missing_peptides <- pd %>% dplyr::filter(!(PID %in% peptide_result$PID))
missing_peptides <- pd %>% dplyr::filter(!(peptide %in% peptide_result$peptide))

# ensure same col names
columns_to_add <- setdiff(names(peptide_result), names(missing_peptides))

# add missing columns to missing df with NA values
missing_peptides[columns_to_add] <- NA

# ensure the column order is same before rbind
missing_peptides <- missing_peptides[, names(peptide_result)]

combined_peptide_result <- rbind(peptide_result, missing_peptides)

combined_peptide_result$PID <- gsub(",", ".", combined_peptide_result$PID)
combined_peptide_result <- combined_peptide_result[!(base::duplicated(combined_peptide_result)),]

# export summary data
write.csv(combined_peptide_result, paste0(output_directory, "/peptide_info.csv"), row.names=F, quote=F)

# include orf_status and peptide_status in GTF mcols
results_to_merge_with_granges <- left_join(results_pept_df, peptide_result, by=c("transcript_id", "peptide", "strand", "PID", "gene_id"))
results_to_merge_with_granges <- results_to_merge_with_granges[!(duplicated(results_to_merge_with_granges)),]
results_to_merge_with_granges <- results_to_merge_with_granges %>% 
  dplyr::select(-openprot_id, -`molecular_weight(kDA)`, -isoelectric_point, -hydrophobicity, -aliphatic_index)
results_to_merge_with_granges$naming <- paste0(results_to_merge_with_granges$transcript_id, "_", results_to_merge_with_granges$peptide)

# make GRanges from df of ORF transcript coordinates
pep_in_genomic_gr_export <- makeGRangesFromDataFrame(results_to_merge_with_granges,
                                                     keep.extra.columns=TRUE, ignore.strand=FALSE, seqinfo=NULL,
                                                     seqnames.field="seqnames", start.field="start", end.field="end", strand.field="strand",
                                                     starts.in.df.are.0based=FALSE, na.rm=TRUE)

names(pep_in_genomic_gr_export) <- c(pep_in_genomic_gr_export$naming) # set names

# add mcols
pep_in_genomic_gr_export$source <- "custom"
pep_in_genomic_gr_export$type <- "exon"
pep_in_genomic_gr_export$group_id <- "peptides"

# export annotations for vis
combined <- c(pep_in_genomic_gr_export, orf_in_genomic_gr, gtf_filtered)
rtracklayer::export(combined, paste0(output_directory, "/combined_annotations.gtf"), format="gtf")

# ---------------------------------------- #



