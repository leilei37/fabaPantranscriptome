library(data.table)

save(newex,gene_ex,file='~/Downloads/legume/isoseq/filter/stablegene_preferisoform.Rdata')
load('~/Downloads/legume/isoseq/filter/stablegene_preferisoform.Rdata')

names(geno)
dt <- iso_m[gene_id == 'PB.23384']
dt <- iso_m[gene_id == 'Vfaba.Hedin2.R2.3g018667']

###########################################################################################
dt[, highlight := ifelse(
  transcript_id == "PB.23384.240",
  "PB.23384.240",
  "Other isoforms"
)]

dt[, highlight := factor(
  highlight,
  levels = c("PB.23384.240", "Other isoforms")
)]
ord_dt <- dt[transcript_id == "PB.23384.240",
             .(IF_240 = IF_norm),
             by = .(geno, accession)]
ord_dt[, acc_order := {
  accession[order(IF_240)]
}, by = geno]
acc_levels <- ord_dt[
  order(geno, IF_240),
  accession
]
dt[, accession := factor(accession, levels = unique(acc_levels))]
dt[, transcript_id := factor(
  transcript_id,
  levels = c(
    "PB.23384.240",
    setdiff(unique(transcript_id), "PB.23384.240")
  )
)]
other_iso <- setdiff(levels(dt$transcript_id), "PB.23384.240")
grey_cols <- gray.colors(
  length(other_iso),
  start = 0.85, 
  end   = 0.4   
)

fill_cols <- c(
  "PB.23384.240" = "pink",
  setNames(grey_cols, other_iso)
)
p_iso <- ggplot(
  dt,
  aes(x = accession, y = IF_norm, fill = transcript_id)
) +
  geom_bar(stat = "identity", width = 1) +
  facet_wrap(~ geno, scales = "free_x", nrow = 1) +
  scale_fill_manual(
    values = fill_cols,
    name = "Isoform"
  ) +
  theme_bw() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(
    title = "Isoform usage (IF)",
    y = "Isoform Fraction"
  )
print(p_iso)
######################################
dt <- dt[which(!is.na(geno)), ]

dt$gene_id <- newex$gene_id[
  match(dt$transcript_id, newex$transcript_id)
]

dt[, highlight := ifelse(
  transcript_id == "PB.23384.240",
  "PB.23384.240",
  "Other isoforms"
)]

dt[, highlight := factor(
  highlight,
  levels = c("PB.23384.240", "Other isoforms")
)]

# ------------------------------
# accession order（based on PB.23384.240 IF）
# ------------------------------

ord_dt <- dt[
  gene_id == "Vfaba.Hedin2.R2.3g018667",
  .(IF_667 = sum(IF_norm, na.rm = TRUE)),
  by = .(geno, accession)
]

ord_dt[, acc_order := accession[order(IF_667)], by = geno]

acc_levels <- ord_dt[
  order(geno, IF_667),
  accession
]

dt[, accession := factor(
  accession,
  levels = unique(acc_levels)
)]
############################################
ord_dt <- dt[
  transcript_id == "PB.23384.240",
  .(IF_240 = IF_norm),
  by = .(geno, accession)
]

ord_dt[, acc_order := {
  accession[order(IF_240)]
}, by = geno]

acc_levels <- ord_dt[
  order(geno, IF_240),
  accession
]

dt[, accession := factor(
  accession,
  levels = unique(acc_levels)
)]
# ------------------------------
# accession order（based on gene 667 first，then based on PB.23384.240）
# ------------------------------

ord_667 <- dt[
  gene_id == "Vfaba.Hedin2.R2.3g018667",
  .(IF_667 = sum(IF_norm, na.rm = TRUE)),
  by = .(geno, accession)
]

ord_240 <- dt[
  transcript_id == "PB.23384.240",
  .(IF_240 = IF_norm),
  by = .(geno, accession)
]

ord_dt <- merge(
  ord_667, ord_240,
  by = c("geno", "accession"),
  all = TRUE
)

ord_dt[is.na(IF_667), IF_667 := 0]
ord_dt[is.na(IF_240), IF_240 := 0]

ord_dt[, acc_order := accession[
  order(IF_667, IF_240)
], by = geno]

acc_levels <- ord_dt[
  order(geno, IF_667, IF_240),
  accession
]

dt[, accession := factor(
  accession,
  levels = unique(acc_levels)
)]

# ------------------------------
# transcript order（240 on the top）
# ------------------------------

tx_240 <- "PB.23384.240"

tx_667_other <- unique(dt$transcript_id[
  dt$gene_id == "Vfaba.Hedin2.R2.3g018667" &
    dt$transcript_id != tx_240
])

tx_666 <- unique(dt$transcript_id[
  dt$gene_id == "Vfaba.Hedin2.R2.3g018666"
])

dt[, transcript_id := factor(
  transcript_id,
  levels = c(
    tx_240,         
    tx_667_other,  
    tx_666          
  )
)]

# ------------------------------
# color structure
# ------------------------------

tx_240 <- "PB.23384.240"

tx_667 <- unique(dt$transcript_id[
  dt$gene_id == "Vfaba.Hedin2.R2.3g018667"
])

tx_666 <- unique(dt$transcript_id[
  dt$gene_id == "Vfaba.Hedin2.R2.3g018666"
])

tx_667_other <- setdiff(tx_667, tx_240)

#FDE0DD
red_cols <- colorRampPalette(c("#CCCCCC", "#CCCCCC"))(
  length(tx_667_other)
)
names(red_cols) <- tx_667_other

grey_cols <- gray.colors(
  length(tx_666),
  start = 0.8,
  end   = 0.8
)
names(grey_cols) <- tx_666

pink_cols <- c("PB.23384.240" = "pink")

fill_cols <- c(pink_cols, red_cols, grey_cols)
fill_cols[2]<-"lightblue"
fill_cols[18]<-"lightgreen"
fill_cols[23]<-'lightgoldenrod'
fill_cols[25]<-'lightcoral'
fill_cols[13]<-'lightcyan'


save(dt, file = '~/Downloads/legume/isoseq/filter/240_plotdata.Rdata')
load('~/Downloads/legume/isoseq/filter/240_plotdata.Rdata')
# ------------------------------
# plot
# ------------------------------
dt.67<-dt%>% filter(gene_id == "Vfaba.Hedin2.R2.3g018667")
p_iso <- ggplot(
  dt,
  aes(x = accession, y = IF_norm, fill = transcript_id)
) +
  geom_bar(stat = "identity", width = 1) +
  facet_wrap(~ geno, scales = "free_x", nrow = 1) +
  scale_fill_manual(
    values = fill_cols,
    name = "Isoform"
  ) +
  theme_bw() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid   = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(
    title = "Isoform usage (IF)",
    y = "Isoform Fraction",
    x = NULL
  )

print(p_iso)

pdf("~/Downloads/fabaisoseqfig_pdf/chr3_385770429isoform_usage_by_genotype_allisoform_rmlowisoform_highlightorfgene_order.pdf", width = 15, height = 5)
print(p_iso)
dev.off()
