cat("\n🚀 Loading genome annotation and core dependency packages...\n")
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

target_files <- list(
    Down_Regulated = "H2A.Z_SpikeIn_DESeq2_NC_vs_siCebpb_Down_Peaks.bed",
    Up_Regulated   = "H2A.Z_SpikeIn_DESeq2_NC_vs_siCebpb_Up_Peaks.bed"
)

bed_files <- list()
for(name in names(target_files)) {
    file_path <- target_files[[name]]
    if(file.exists(file_path) && file.info(file_path)$size > 0) {
        bed_files[[name]] <- file_path
    } else {
        cat(paste0("⚠️  Warning: Valid BED file not found for ", name, " (", file_path, "). Skipping...\n"))
    }
}

if(length(bed_files) == 0) {
    stop("Error: No valid differential BED files detected in the current directory.")
}

cat("📂 Reading BED files and converting to GenomicRanges objects...\n")
peak_list <- lapply(bed_files, readPeakFile)

cat("⚡ Running ChIPseeker core annotation engine with unified TSS ± 5kb window...\n")
anno_list <- lapply(peak_list, annotatePeak, 
                    TxDb = txdb, 
                    tssRegion = c(-5000, 5000), 
                    annoDb = "org.Mm.eg.db")

cat("💾 Applying strict 5kb absolute distance filtering and exporting tables...\n")
for(name in names(anno_list)) {
    anno_df <- as.data.frame(anno_list[[name]])
    
    anno_df_filtered <- anno_df[!is.na(anno_df$distanceToTSS) & abs(anno_df$distanceToTSS) <= 5000, ]
    
    core_columns <- c("seqnames", "start", "end", "V4", "V5", "annotation", "geneId", "SYMBOL", "GENENAME", "distanceToTSS")
    existing_cols <- core_columns[core_columns %in% colnames(anno_df_filtered)]
    remaining_cols <- setdiff(colnames(anno_df_filtered), existing_cols)
    
    anno_df_final <- anno_df_filtered[, c(existing_cols, remaining_cols)]
    anno_df_final <- anno_df_final[order(abs(anno_df_final$distanceToTSS)), ]
    
    output_txt <- paste0("H2A.Z_SpikeIn_DESeq2_NC_vs_siCebpb_", name, "_Peaks_Within_5kb.txt")
    write.table(anno_df_final, output_txt, sep = "\t", quote = FALSE, row.names = FALSE)
    
    unique_genes <- unique(anno_df_final$SYMBOL)
    unique_genes <- unique_genes[!is.na(unique_genes) & unique_genes != ""]
    
    output_genes <- paste0("H2A.Z_SpikeIn_DESeq2_NC_vs_siCebpb_", name, "_Within_5kb_Unique_Genes.txt")
    writeLines(unique_genes, output_genes)
    
    cat(paste0("  🎉 Success [", name, "]: Captured ", nrow(anno_df_final), " peaks within TSS ± 5kb, mapping to ", length(unique_genes), " unique target genes.\n"))
}

save.image("h2az_cutandtag_final_analysis_session_normalized.RData")
cat("\n==================================================================\n")
cat("🎉 ChIPseeker annotation pipeline completed successfully!\n")
cat("==================================================================\n")
