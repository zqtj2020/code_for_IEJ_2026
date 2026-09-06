load("h2az_safe_spikein_session.RData")
library(DiffBind)
library(DESeq2)

cat("\n5-1. Extracting raw count matrix and metadata...\n")
counts_matrix <- dba.peakset(h2az_counts, bRetrieve = TRUE)
counts_df     <- as.data.frame(mcols(counts_matrix))

sample_names  <- h2az_counts$samples$SampleID
raw_counts    <- as.matrix(counts_df[, sample_names])
rownames(raw_counts) <- paste0(seqnames(counts_matrix), ":", start(counts_matrix), "-", end(counts_matrix))

col_data <- data.frame(
    Condition = factor(h2az_counts$samples$Condition, levels = c("Control", "Knockdown")), 
    row.names = sample_names
)

cat("5-2. Injecting external Spike-in size factors and building DESeq2 object...\n")
dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = col_data, design = ~ Condition)

if(!exists("size_factors")){
    stop("Error: External normalisation variable 'size_factors' not found in RData!")
}
if(!all(names(size_factors) %in% colnames(dds))) {
    if(length(size_factors) == ncol(dds)){
        names(size_factors) <- colnames(dds)
    } else {
        stop("Error: The number of size_factors does not match the sample count!")
    }
}
sizeFactors(dds) <- size_factors

cat("5-3. Running DESeq2 statistical testing and apeglm LFC shrinkage...\n")
dds <- DESeq(dds)
res_deseq2 <- lfcShrink(dds, coef = "Condition_Knockdown_vs_Control", type = "apeglm")

cat("5-4. Merging results and parsing genomic coordinates...\n")
res_df <- as.data.frame(res_deseq2)
coords <- do.call(rbind, strsplit(rownames(res_df), "[:-]"))
res_df$seqnames <- as.character(coords[, 1])
res_df$start    <- as.integer(coords[, 2])
res_df$end      <- as.integer(coords[, 3])
res_df$FDR      <- res_df$padj  
res_df$Log2FoldChange <- res_df$log2FoldChange

cat("5-5. Performing dual threshold filtering (FDR < 0.05, Fold-Change > 1.2)...\n")
valid_res  <- res_df[!is.na(res_df$FDR), ]
down_peaks <- valid_res[valid_res$FDR < 0.05 & valid_res$Log2FoldChange < -0.263, ]
up_peaks   <- valid_res[valid_res$FDR < 0.05 & valid_res$Log2FoldChange > 0.263, ]

write_bed <- function(df, prefix, filename) {
    if(nrow(df) > 0) {
        bed <- data.frame(
            chrom  = as.character(df$seqnames),
            start  = as.integer(df$start) - 1, 
            end    = as.integer(df$end),
            name   = paste0(prefix, "_", seq_len(nrow(df))),
            score  = pmin(round(abs(df$Log2FoldChange) * 100), 1000), 
            strand = rep(".", nrow(df))
        )
        write.table(bed, filename, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        cat(paste0("  --> Successfully written: ", filename, " (Total regions: ", nrow(df), ")\n"))
    } else {
        cat(paste0("  ⚠️ Warning: ", prefix, " filtered result is empty, file not generated.\n"))
    }
}

cat("5-6. Generating standardized differential Peaks files...\n")
write_bed(down_peaks, "CutTag_Down_Peak", "H2A.Z_SpikeIn_DESeq2_NC_vs_siCebpb_Down_Peaks.bed")
write_bed(up_peaks, "CutTag_Up_Peak", "H2A.Z_SpikeIn_DESeq2_NC_vs_siCebpb_Up_Peaks.bed")

cat("\n==================================================================")
cat("\n🎉 DESeq2 calculation pipeline executed successfully!")
cat("\n📂 Output files deposited in the current directory:")
cat("\n  - Down-regulated: H2A.Z_SpikeIn_DESeq2_NC_vs_siCebpb_Down_Peaks.bed")
cat("\n  - Up-regulated: H2A.Z_SpikeIn_DESeq2_NC_vs_siCebpb_Up_Peaks.bed")
cat("\n==================================================================\n")
