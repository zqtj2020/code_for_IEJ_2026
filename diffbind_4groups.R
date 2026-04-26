library(DiffBind)

dbObj <- dba(sampleSheet="diffbind_full.csv")
dbObj <- dba.count(dbObj, minOverlap=1, bParallel=TRUE)

run_logic_full <- function(mask, label) {
    # Subset object based on mask
    sub <- dba(dbObj, mask=mask)
    d3_idx <- which(sub$masks$D3)
    d0_idx <- which(sub$masks$D0)
    
    sub <- dba.contrast(sub, group1=d3_idx, group2=d0_idx, name1="D3", name2="D0")
    sub <- dba.normalize(sub, method=DBA_DESEQ2, normalize=DBA_NORM_TMM)
    sub <- dba.analyze(sub)
    
    res <- as.data.frame(dba.report(sub, th=1))
    d3_enriched <- res[res$p.value < 0.05 & res$Fold > 0, ]
    d0_enriched <- res[res$p.value < 0.05 & res$Fold < 0, ]
    
    file_label <- if(label == "siC") "siCebpb" else label
    
    write.table(d3_enriched[,1:3], paste0("D3_", file_label, "_DARs.bed"), sep="\t", quote=F, row.names=F, col.names=F)
    write.table(d0_enriched[,1:3], paste0("D0_", file_label, "_DARs.bed"), sep="\t", quote=F, row.names=F, col.names=F)
    
    return(c(up_in_D3=nrow(d3_enriched), up_in_D0=nrow(d0_enriched)))
}

nc_res <- run_logic_full(dbObj$masks$NC, "NC")
sic_res <- run_logic_full(dbObj$masks$siC, "siC")

