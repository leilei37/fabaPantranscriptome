library(data.table)

load('~/Downloads/legume/isoseq/filter/stablegene_preferisoform.Rdata')
acc_cols <- setdiff(colnames(gene_ex), "gene_id")
#normalized gene_ex
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = gene_ex[, -1],
                              colData = data.frame(row.names = acc_cols),
                              design = ~ 1)
dds <- estimateSizeFactors(dds)
norm_gene_ex <- counts(dds, normalized=TRUE)
norm_gene_ex<-cbind(gene_ex[,1,drop=FALSE],norm_gene_ex)
# ---------- gene: wide -> long ----------
gene_long <- melt(
  norm_gene_ex,
  id.vars = "gene_id",
  measure.vars = acc_cols,
  variable.name = "accession",
  value.name = "expr"
)

# ---------- isoform: wide -> long ----------
acc_cols_iso <- setdiff(colnames(newex), c("gene_id", "transcript_id"))

iso_long <- melt(
  newex,
  id.vars = c("gene_id", "transcript_id"),
  measure.vars = acc_cols_iso,
  variable.name = "accession",
  value.name = "expr"
)

setDT(gene_long)
setDT(iso_long)

gene_long[, expr := as.numeric(expr)]
iso_long[,  expr := as.numeric(expr)]

gene_dt <- gene_long
iso_dt  <- iso_long

check <- iso_long[, .(iso_sum = sum(expr, na.rm=TRUE)),
                  by=.(gene_id, accession)]
iso_long$gene_id<-as.character(iso_long$gene_id)
check <- merge(check, gene_long, by=c("gene_id","accession"))
check[, ratio := iso_sum / expr]

summary(check$ratio)

# Expression filters (adjust as needed)
min_gene_mean_expr <- 1      # TPM>=1 typically; if counts, raise depending library size
min_iso_mean_expr  <- 1    # TPM>=0.1 (or counts >= e.g. 5)
min_accessions_expressed <- 5  # require gene/isoform measurable in enough accessions

# Gene stability thresholds
gene_cv_max <- 0.20          # stable total expression across accessions
gene_iqr_over_median_max <- 0.30  # alternative robust stability threshold

# Isoform usage-change thresholds
min_if_mean <- 0.05          # keep isoforms contributing at least 5% on average
min_if_range <- 0.30         # max(IF)-min(IF) >= 0.30 indicates strong usage change
min_if_sd <- 0.10            # sd(IF) >= 0.10
min_dom_switch_frac <- 0.20  # dominant isoform differs in >=20% accessions
min_entropy_sd <- 0.15       # optional: entropy sd threshold

# Handling zeros in entropy
eps <- 1e-12
# --------------------------------

# ---------- helpers ----------
cv <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (is.na(m) || m == 0) return(NA_real_)
  s / m
}

iqr_over_median <- function(x) {
  med <- median(x, na.rm = TRUE)
  if (is.na(med) || med == 0) return(NA_real_)
  IQR(x, na.rm = TRUE) / med
}

shannon_entropy <- function(p) {
  # p: numeric vector, should sum to 1 (or close)
  p <- p[!is.na(p) & p > 0]
  if (length(p) == 0) return(NA_real_)
  -sum(p * log(p))
}

# normalize column names (flexible)
setnames(gene_dt, tolower(names(gene_dt)))
setnames(iso_dt,  tolower(names(iso_dt)))

# Expected columns:
# gene_dt: gene_id, accession, expr
# iso_dt:  gene_id, transcript_id, accession, expr
stopifnot(all(c("gene_id","accession","expr") %in% names(gene_dt)))
stopifnot(all(c("gene_id","transcript_id","accession","expr") %in% names(iso_dt)))

# Ensure types
gene_dt[, expr := as.numeric(expr)]
iso_dt[,  expr := as.numeric(expr)]

# ---------- basic filtering ----------
# keep genes expressed in enough accessions
gene_expr_n <- gene_dt[expr > 0, .(n_acc = uniqueN(accession), mean_expr = mean(expr, na.rm=TRUE)), by=gene_id]
gene_keep <- gene_expr_n[n_acc >= min_accessions_expressed & mean_expr >= min_gene_mean_expr, gene_id]
gene_dt <- gene_dt[gene_id %in% gene_keep]

# keep isoforms for those genes
iso_dt <- iso_dt[gene_id %in% gene_keep]

# keep isoforms expressed in enough accessions and above mean
iso_expr_n <- iso_dt[expr > 0, .(n_acc = uniqueN(accession), mean_expr = mean(expr, na.rm=TRUE)), by=.(gene_id, transcript_id)]
iso_keep <- iso_expr_n[n_acc >= min_accessions_expressed & mean_expr >= min_iso_mean_expr]
iso_dt <- iso_dt[iso_keep, on=.(gene_id, transcript_id), nomatch=0]

# ---------- compute Isoform Fraction (IF) ----------
# Merge gene expression onto isoform table by gene_id+accession
setkey(gene_dt, gene_id, accession)
setkey(iso_dt,  gene_id, accession)

iso_m <- iso_dt[gene_dt, nomatch=0]
# iso_m has: gene_id, transcript_id, accession, expr (iso), i.expr (gene)
setnames(iso_m, c("expr","i.expr"), c("iso_expr","gene_expr"))

# IF: iso / gene; protect divide-by-zero
iso_m[, IF := fifelse(is.na(gene_expr) | gene_expr <= 0, NA_real_, iso_expr / gene_expr)]

# optional: within each gene+accession, renormalize IF to sum to 1
# (helps when gene_expr not exactly sum iso_expr due to quant method)
iso_m[, IF_norm := {
  s <- sum(IF, na.rm=TRUE)
  if (is.na(s) || s <= 0) rep(NA_real_, .N) else IF / s
}, by=.(gene_id, accession)]

# choose IF used for stats
iso_m[, IF_use := IF_norm]

# ---------- gene stability metrics ----------
gene_stats <- gene_dt[, .(
  n_acc = uniqueN(accession),
  mean_expr = mean(expr, na.rm=TRUE),
  median_expr = median(expr, na.rm=TRUE),
  sd_expr = sd(expr, na.rm=TRUE),
  cv_expr = cv(expr),
  iqr_med = iqr_over_median(expr)
), by=gene_id]

quantile(gene_stats$cv_expr, probs = seq(0, 1, 0.05), na.rm=TRUE)
cv_cutoff  <- quantile(gene_stats$cv_expr, 0.20, na.rm=TRUE)
iqr_cutoff <- quantile(gene_stats$iqr_med, 0.20, na.rm=TRUE)

stable_genes <- gene_stats[
  cv_expr <= cv_cutoff, 
  gene_id
]
#stable_genes <- gene_stats[
#  cv_expr <= gene_cv_max & iqr_med <= gene_iqr_over_median_max, gene_id
#]

# ---------- isoform usage variability metrics ----------
# keep only stable genes
iso_s <- iso_m[gene_id %in% stable_genes & !is.na(IF_use)]
#iso_s <- iso_stats[gene_id %in% stable_genes & !is.na(IF_use)]

# drop tiny-contribution isoforms (mean IF)
iso_mean_if <- iso_s[, .(mean_if = mean(IF_use, na.rm=TRUE)), by=.(gene_id, transcript_id)]
iso_s <- iso_s[iso_mean_if[mean_if >= min_if_mean], on=.(gene_id, transcript_id), nomatch=0]

# isoform-level IF variability
iso_if_stats <- iso_s[, .(
  mean_if = mean(IF_use, na.rm=TRUE),
  sd_if   = sd(IF_use, na.rm=TRUE),
  min_if  = min(IF_use, na.rm=TRUE),
  max_if  = max(IF_use, na.rm=TRUE),
  range_if = max(IF_use, na.rm=TRUE) - min(IF_use, na.rm=TRUE),
  n_acc = uniqueN(accession)
), by=.(gene_id, transcript_id)]

# dominant isoform per accession (within gene)
dom <- iso_s[order(gene_id, accession, -IF_use)]
dom <- dom[, .SD[1], by=.(gene_id, accession)]
setnames(dom, "transcript_id", "dom_transcript")

# dominant switching summary per gene
dom_switch <- dom[, .(
  n_acc = uniqueN(accession),
  n_dom = uniqueN(dom_transcript),
  top_dom = names(sort(table(dom_transcript), decreasing=TRUE))[1],
  top_dom_frac = max(table(dom_transcript)) / .N
), by=gene_id]
dom_switch[, dom_switch_frac := 1 - top_dom_frac]  # fraction NOT using top dominant isoform

# entropy per gene+accession
# build vector of IF within gene+accession
entropy_dt <- iso_s[, .(
  entropy = shannon_entropy(IF_use + eps)
), by=.(gene_id, accession)]

entropy_gene <- entropy_dt[, .(
  entropy_mean = mean(entropy, na.rm=TRUE),
  entropy_sd   = sd(entropy, na.rm=TRUE),
  entropy_min  = min(entropy, na.rm=TRUE),
  entropy_max  = max(entropy, na.rm=TRUE)
), by=gene_id]

# gene-level usage-change summary: based on isoform stats + dominant switching + entropy
# 1) at least 2 isoforms kept
n_iso <- iso_if_stats[, .(n_iso = uniqueN(transcript_id)), by=gene_id]

# 2) gene-level max range/sd among isoforms
gene_usage <- iso_if_stats[, .(
  max_range_if = max(range_if, na.rm=TRUE),
  max_sd_if    = max(sd_if, na.rm=TRUE)
), by=gene_id]

# merge everything
res <- merge(gene_stats, n_iso, by="gene_id", all.x=FALSE, all.y=FALSE)
res <- merge(res, gene_usage, by="gene_id", all.x=TRUE)
res <- merge(res, dom_switch, by="gene_id", all.x=TRUE)
res <- merge(res, entropy_gene, by="gene_id", all.x=TRUE)

# apply final criteria for "stable gene but isoform usage changed"
hit <- res[
  n_iso >= 2 &
    #max_range_if >= min_if_range &
    max_sd_if >= min_if_sd 
    #dom_switch_frac >= min_dom_switch_frac 
    #entropy_sd >= min_entropy_sd
]
hist(res$max_sd_if,breaks=100)
hist(res$max_range_if,breaks=100)

# rank: stronger usage change first, then more dominant switching
setorder(hit, -max_range_if, -dom_switch_frac, -entropy_sd)
#which(hit$gene_id %in% ssp.gene$HedinV2)

# ---------- write outputs ----------
fwrite(hit, "stable_gene_isoform_change.tsv", sep="\t")

# also output long table for plotting
plot_dt <- iso_s[gene_id %in% hit$gene_id, .(gene_id, transcript_id, accession, gene_expr, iso_expr, IF=IF_use,IF_norm)]
fwrite(plot_dt, "stable_gene_isoform_change_long_if.tsv", sep="\t")

message("Done.\nOutputs:\n  stable_gene_isoform_change.tsv\n  stable_gene_isoform_change_long_if.tsv")

stable_set <- unique(stable_genes)
hit_set    <- unique(hit$gene_id)

n_stable <- length(stable_set)
n_iso_diff <- length(intersect(stable_set, hit_set))
prop <- n_iso_diff / n_stable

cat("Stable genes:", n_stable, "\n")
cat("Isoform-variable among stable:", n_iso_diff, "\n")
cat("Proportion:", round(prop, 4), "\n")

######################################plot
library(viridis)
p3_hex <- ggplot(plot_df, aes(x = cv_expr, y = dom_sd_usage)) +
  geom_hex(bins = 100) +
  scale_fill_gradient(low = "lightblue", high = "darkorange") +  # 自定义渐变
  geom_vline(xintercept = cv_cutoff, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = min_if_sd, linetype = "dashed", color = "grey40") +
  labs(
    x = "Gene-level expression variation (CV)",
    y = "Isoform usage variation (SD of dominant isoform proportion)",
    fill = "Count"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  coord_cartesian(
    xlim = c(0, quantile(plot_df$cv_expr, 0.99, na.rm = TRUE))
  )
pdf("~/Downloads/fabaisoseqfig_pdf/stablegene_preferisoform_summary_hex.pdf", width = 5, height = 7)
print(p3_hex)
dev.off()


