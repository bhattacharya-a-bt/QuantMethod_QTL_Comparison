library(data.table)
library(ggplot2)
library(RColorBrewer)
library(ggVennDiagram)
library(dplyr)
library(tidyr)
library(rtracklayer)
library(patchwork)
library(GenomicRanges)
library(dplyr)
library(purrr)
library(cowplot)
library(UpSetR)

### read in data

dat <- {}
for(chr in 1:22){
  load(paste0("/Users/sthead/OneDrive - Inside MD Anderson/Bhattacharya,Arjun's files - GTEx GENCODE Comp/5_simulations/results/pass2/eqtl_res_with_true_beta_chr",chr,".RData"))
  out <- data.frame(out)
  out$is_qtl <- lengths(out$transcript_id) > 0
  out$gene_trim <-sub("\\..*$", "", out$phe_id)
  out <- out[out$adj_beta_pval < 0.05,]
  out$max_ld_with_true_eqtls[out$is_qtl==T] <- 1
  out$shared_eqtl <- (lengths(out$transcript_id) > 1)
  dat <- rbind(dat,out)
  rm(out)
}

dat$setting <- paste(dat$quant,dat$annot,dat$measure,dat$norm,sep="/")
settings <- unique(dat$setting)




safe_summary <- function(x, prefix){
  out <- c(
    mean = NA,
    sd   = NA,
    min  = NA,
    max  = NA,
    q1=NA,
    q2=NA,
    med=NA
  )
  
  if(length(x) > 0 && any(!is.na(x))){
    out["mean"] <- mean(x, na.rm=TRUE)
    out["sd"]   <- sd(x, na.rm=TRUE)
    out["min"]  <- min(x, na.rm=TRUE)
    out["max"]  <- max(x, na.rm=TRUE)
    out["q1"] <- quantile(x, 0.25, na.rm=T)
    out["q2"] <- quantile(x, 0.75, na.rm=T)
    out["med"] <- quantile(x, 0.5, na.rm=T)
  }
  
  names(out) <- paste0(prefix, names(out))
  return(out)
}

extract_tx_features <- function(x_list){
  x <- unlist(x_list)
  
  if(length(x) == 0 || all(is.na(x))){
    return(c(mean=NA, max=NA, sd=NA))
  }
  
  return(c(
    mean = mean(x, na.rm=TRUE),
    max  = max(x, na.rm=TRUE),
    sd   = sd(x, na.rm=TRUE)
  ))
}

out <- list()
counter <- 1

for(i in 1:length(settings)){
  for(j in 1:length(settings)){
    
    if(j > i){
      
      s1 <- settings[i]
      s2 <- settings[j]
      
      df <- data.frame(s1=s1, s2=s2)
      
      d1 <- dat[dat$setting == s1, ]
      d2 <- dat[dat$setting == s2, ]
      
      # -----------------------------
      # eGene counts
      # -----------------------------
      
      df$n_egene1 <- length(unique(d1$gene_trim))
      df$n_egene2 <- length(unique(d2$gene_trim))
      
      df$n_true_egene1 <- sum(d1$true_egene)
      df$n_true_egene2 <- sum(d2$true_egene)
      
      # concordance
      egene_concord <- intersect(d1$gene_trim, d2$gene_trim)
      df$n_egene_concord <- length(egene_concord)
      
      true_egene_concord <- intersect(
        d1$gene_trim[d1$true_egene],
        d2$gene_trim[d2$true_egene]
      )
      df$n_true_egene_concord <- length(true_egene_concord)
      
      # -----------------------------
      # TRUE eQTL recovery
      # -----------------------------
      d1_true <- d1[d1$true_egene,]
      d2_true <- d2[d2$true_egene,]
      
      df$prop_true_qtl1 <- if(nrow(d1_true)>0) mean(d1_true$is_qtl) else NA
      df$prop_true_qtl2 <- if(nrow(d2_true)>0) mean(d2_true$is_qtl) else NA
      
      # ---------------------------------
      # max LD of lead SNP with true QTLs
      # ---------------------------------
      d1_true <- d1[d1$true_egene,]
      d2_true <- d2[d2$true_egene,]
      
      df$max_ld_qtl1_mean <- if(nrow(d1_true)>0) mean(d1_true$max_ld_with_true_eqtls,na.rm=T) else NA
      df$max_ld_qtl2_mean <- if(nrow(d2_true)>0) mean(d2_true$max_ld_with_true_eqtls,na.rm=T) else NA
      
      df$max_ld_qtl1_med <- if(nrow(d1_true)>0) median(d1_true$max_ld_with_true_eqtls,na.rm=T) else NA
      df$max_ld_qtl2_med <- if(nrow(d2_true)>0) median(d2_true$max_ld_with_true_eqtls,na.rm=T) else NA
      
      df$max_ld_qtl1_25 <- if(nrow(d1_true)>0) quantile(d1_true$max_ld_with_true_eqtls,.25,na.rm=T) else NA
      df$max_ld_qtl2_25 <- if(nrow(d2_true)>0) quantile(d2_true$max_ld_with_true_eqtls,.25,na.rm=T) else NA
      
      df$max_ld_qtl1_75 <- if(nrow(d1_true)>0) quantile(d1_true$max_ld_with_true_eqtls,.75,na.rm=T) else NA
      df$max_ld_qtl2_75 <- if(nrow(d2_true)>0) quantile(d2_true$max_ld_with_true_eqtls,.75,na.rm=T) else NA
      
      # -----------------------------
      # Align concordant genes
      # -----------------------------
      #rownames(d1) <- d1$gene_trim
      #rownames(d2) <- d2$gene_trim
      
      tmp1 <- d1[d1$gene_trim %in% true_egene_concord, ]
      tmp2 <- d2[d2$gene_trim %in% true_egene_concord, ]
      
      # -----------------------------
      # Top variant concordance
      # -----------------------------
      if(nrow(tmp1) > 0){
        same_var <- tmp1$var_id == tmp2$var_id
        df$prop_topvar_concord <- mean(same_var, na.rm=TRUE)
        df$n_topvar_concord <- sum(same_var, na.rm=TRUE)
      } else {
        df$prop_topvar_concord <- NA
        df$n_topvar_concord <- NA
      }
      
      # -----------------------------
      # Concordant gene features
      # -----------------------------
      df <- cbind(df,
        as.data.frame(t(safe_summary(tmp1$n_transcripts, "concord_ntx_"))),
        as.data.frame(t(safe_summary(tmp1$gene_h2, "concord_h2_"))),
        as.data.frame(t(safe_summary(tmp1$n_qtl_shared, "concord_nqtl_shared")))
      )
      
      # among concordant true egenes, proportion of top QTLs that are shared QTLs
      df$concord_prop_sharedqtl1 <- mean(tmp1$shared_eqtl, na.rm=TRUE)
      df$concord_prop_sharedqtl2 <- mean(tmp2$shared_eqtl, na.rm=TRUE)
      
      # -----------------------------
      # Transcript-level features
      # -----------------------------
      if(nrow(tmp1) > 0){
        tx_stats <- t(sapply(tmp1$transcript_h2, extract_tx_features))
        
        df$tx_h2_mean_concord <- mean(tx_stats[,1], na.rm=TRUE)
        df$tx_h2_max_concord  <- mean(tx_stats[,2], na.rm=TRUE)
        df$tx_h2_sd_concord   <- mean(tx_stats[,3], na.rm=TRUE)

              } else {
        df$tx_h2_mean_concord <- NA
        df$tx_h2_max_concord  <- NA
        df$tx_h2_sd_concord   <- NA
      }
      
      # -----------------------------
      # Discordant genes (detection)
      # -----------------------------
      disc1 <- d1[!(d1$gene_trim %in% d2$gene_trim) & d1$true_egene, ]
      disc2 <- d2[!(d2$gene_trim %in% d1$gene_trim) & d2$true_egene, ]
      
      df <- cbind(df,
        as.data.frame(t(safe_summary(disc1$n_transcripts, "disc_ntx_"))),
        as.data.frame(t(safe_summary(disc1$gene_h2, "disc_h2_"))),
        as.data.frame(t(safe_summary(disc1$n_qtl_shared, "disc_nqtl_")))
      )
      
      # among discordant true egenes, proportion of top QTLs that are shared QTLs
      df$discord_prop_sharedqtl1 <- mean(disc1$shared_eqtl, na.rm=TRUE)
      df$discord_prop_sharedqtl2 <- mean(disc2$shared_eqtl, na.rm=TRUE)
      
      # -----------------------------
      # Variant discordance (same gene, diff SNP)
      # -----------------------------
      if(nrow(tmp1) > 0){
        keep <- which(tmp1$var_id != tmp2$var_id)
        
        vdisc1 <- tmp1[keep, , drop=FALSE]
        
        df <- cbind(df,
          as.data.frame(t(safe_summary(vdisc1$r_squared, "vdisc1_r2_")))
        )
      }
      
      if(nrow(tmp2) > 0){
        keep <- which(tmp2$var_id != tmp1$var_id)
        
        vdisc2 <- tmp2[keep, , drop=FALSE]
        
        df <- cbind(df,
          as.data.frame(t(safe_summary(vdisc2$r_squared, "vdisc2_r2_")))
        )
      }
      
      # -----------------------------
      # Save
      # -----------------------------
      out[[counter]] <- df
      counter <- counter + 1
    }
  }
}

out_df <- do.call(rbind, out)

out_df <- out_df[-grep("normalized",out_df$s1),]
out_df <- out_df[-grep("normalized",out_df$s2),]

out_df <- out_df[-grep("abundance",out_df$s1),]
out_df <- out_df[-grep("abundance",out_df$s2),]

##########################################################
# begin simulation main figure 

# plot comparings heritability between discordant and concordant genes
parse_setting <- function(x) {
  parts <- strsplit(x, "/")
  data.frame(
    method     = sapply(parts, `[`, 1),
    annot      = sapply(parts, `[`, 2),
    norm_method = sapply(parts, `[`, 3),
    inv_rnorm  = sapply(parts, `[`, 4)
  )
}

s1_parsed <- parse_setting(out_df$s1)
s2_parsed <- parse_setting(out_df$s2)

plot_df <- cbind(out_df, 
                 s1_method = s1_parsed$method,
                 s1_annot  = s1_parsed$annot,
                 s1_norm = s1_parsed$inv_rnorm,
                 s2_method = s2_parsed$method,
                 s2_annot  = s2_parsed$annot,
                 s2_norm = s2_parsed$inv_rnorm)

plot_df$label_new <- paste0(plot_df$s1_method,"/",plot_df$s1_annot," vs. ",plot_df$s2_method,"/",plot_df$s2_annot)

# holding method constant
box_df1 <- map_dfr(c("salmon", "featureCounts"), function(method) {
  tmp <- plot_df[plot_df$s1_norm == "normY" & plot_df$s2_norm == "normY", ]
  tmp <- tmp[tmp$s1_method == method & tmp$s2_method == method, ]
  tmp <- tmp[-grep("sub", tmp$label_new), ]
  tmp$label_new <- paste0(tmp$s1_annot, " vs. ", tmp$s2_annot)
  
  tmp %>%
    select(label_new,
           concord_h2_min, concord_h2_q1, concord_h2_q2, concord_h2_max, concord_h2_mean, concord_h2_med,
           disc_h2_min,    disc_h2_q1,    disc_h2_q2,    disc_h2_max,    disc_h2_mean,    disc_h2_med) %>%
    pivot_longer(-label_new,
                 names_to = c("type", ".value"),
                 names_pattern = "(concord|disc)_h2_(.*)") %>%
    mutate(type   = recode(type, concord = "Concordant", disc = "Discordant"),
           method = method)
})
box_df1$variable <- "Gene heritability"
names(box_df1)[names(box_df1)=="method"] <- "constant"


# holding annot constant
box_df2 <- map_dfr(c("gencode_v27", "gencode_v38","gencode_v45"), function(annot) {
  tmp <- plot_df[plot_df$s1_norm == "normY" & plot_df$s2_norm == "normY", ]
  tmp <- tmp[tmp$s1_annot == annot & tmp$s2_annot == annot, ]
  tmp$label_new <- paste0(tmp$s1_method, " vs. ", tmp$s2_method)
  
  tmp %>%
    select(label_new,
           concord_h2_min, concord_h2_q1, concord_h2_q2, concord_h2_max, concord_h2_mean, concord_h2_med,
           disc_h2_min,    disc_h2_q1,    disc_h2_q2,    disc_h2_max,    disc_h2_mean,    disc_h2_med) %>%
    pivot_longer(-label_new,
                 names_to = c("type", ".value"),
                 names_pattern = "(concord|disc)_h2_(.*)") %>%
    mutate(type   = recode(type, concord = "Concordant", disc = "Discordant"),
           annot = annot)
})
box_df2$variable <- "Gene heritability"
names(box_df2)[names(box_df2)=="annot"] <- "constant"

# comparing n tx between discordant and concordant genes

# holding method constant
box_df3 <- map_dfr(c("salmon", "featureCounts"), function(method) {
  tmp <- plot_df[plot_df$s1_norm == "normY" & plot_df$s2_norm == "normY", ]
  tmp <- tmp[tmp$s1_method == method & tmp$s2_method == method, ]
  tmp <- tmp[-grep("sub", tmp$label_new), ]
  tmp$label_new <- paste0(tmp$s1_annot, " vs. ", tmp$s2_annot)
  
  tmp %>%
    select(label_new,
        concord_ntx_min, concord_ntx_q1, concord_ntx_q2, concord_ntx_max, concord_ntx_mean,concord_ntx_med,
         disc_ntx_min,    disc_ntx_q1,    disc_ntx_q2,    disc_ntx_max,    disc_ntx_mean, disc_ntx_med) %>%
    pivot_longer(-label_new,
                 names_to = c("type", ".value"),
                 names_pattern = "(concord|disc)_ntx_(.*)") %>%
    mutate(type   = recode(type, concord = "Concordant", disc = "Discordant"),
           method = method)
})
box_df3$variable <- "Gene isoform count"
names(box_df3)[names(box_df3)=="method"] <- "constant"

# holding annot constant
box_df4 <- map_dfr(c("gencode_v27", "gencode_v38","gencode_v45"), function(annot) {
  tmp <- plot_df[plot_df$s1_norm == "normY" & plot_df$s2_norm == "normY", ]
  tmp <- tmp[tmp$s1_annot == annot & tmp$s2_annot == annot, ]
  tmp$label_new <- paste0(tmp$s1_method, " vs. ", tmp$s2_method)
  
  tmp %>%
    select(label_new,
         concord_ntx_min, concord_ntx_q1, concord_ntx_q2, concord_ntx_max, concord_ntx_mean,concord_ntx_med,
         disc_ntx_min,    disc_ntx_q1,    disc_ntx_q2,    disc_ntx_max,    disc_ntx_mean, disc_ntx_med) %>%
    pivot_longer(-label_new,
                 names_to = c("type", ".value"),
                 names_pattern = "(concord|disc)_ntx_(.*)") %>%
    mutate(type   = recode(type, concord = "Concordant", disc = "Discordant"),
           annot = annot)
})
box_df4$variable <- "Gene isoform count"
names(box_df4)[names(box_df4)=="annot"] <- "constant"

box_df <- rbind(box_df1,box_df2,box_df3,box_df4)
box_df$comp <- "Between annotations"
box_df$comp[box_df$label_new=="featureCounts vs. salmon"] <- "Between methods"

#############################
# main simulations figure
#############################

annot_labels <- c(
  "gencode_v27" = "v27",
  "gencode_v38" = "v38",
  "gencode_v45" = "v45",
  "salmon" = "Salmon"
)

box_df <- box_df %>%
  mutate(label_new = stringr::str_replace_all(label_new, annot_labels))

box_df <- box_df %>%
mutate(constant = stringr::str_replace_all(constant, annot_labels))

box_df$constant <- factor(box_df$constant,levels=c("featureCounts","Salmon","v27","v38","v45"))



# first boxplot
tmp <- box_df[box_df$label_new=="featureCounts vs. Salmon",]

gg1 <- ggplot(tmp, aes(x = label_new, y = med, color = type, shape = constant, group = interaction(constant, type))) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = q1, ymax = q2), width = 0.2, position = position_dodge(width = 0.6)) +
  facet_wrap(~ variable, scales = "free_y") +
  # scale_shape_manual(values = c("Salmon" = 16, "featureCounts" = 17, "v27" = 1, "v38" = 2, "v45" = 0)) +
  scale_shape_manual(values = c("Salmon" = 16, "featureCounts" = 17, "v27" = 1, "v38" = 23, "v45" = 25)) +
  scale_color_manual(values = c("Concordant" = "#0072B2", "Discordant" = "#E69F00")) +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(x="Methods",title = "eGene comparison between methods", y = "Value", color = "eGene type", shape = "GENCODE")

tmp <- box_df[!box_df$label_new=="featureCounts vs. Salmon",]

gg2 <- ggplot(tmp, aes(x = label_new, y = med, color = type, shape = constant, group = interaction(constant, type))) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = q1, ymax = q2), width = 0.2, position = position_dodge(width = 0.6)) +
  facet_wrap(~ variable, scales = "free_y") +
  scale_shape_manual(values = c("Salmon" = 16, "featureCounts" = 17, "GENCODEv27" = 1, "GENCODEv38" = 2, "GENCODEv45" = 0)) +
  scale_color_manual(values = c("Concordant" = "#0072B2", "Discordant" = "#E69F00")) +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(x="Annotations",title = "eGene comparison between annotations", y = "Value", color = "eGene type", shape = "Method")


fin_fig <- cowplot::plot_grid(
  gg1, 
  gg2,
  nrow = 2,
  ncol = 1,
  rel_widths = c(1, 1),
  labels = c("A", "B"),
  label_x = 0.0,
  label_y = .995,
  hjust = 0,
  vjust = 1,
  label_size = 16,
  label_fontface = "bold")

ggsave(plot = fin_fig,
       filename = "/Users/sthead/OneDrive - Inside MD Anderson/Bhattacharya,Arjun's files - GTEx GENCODE Comp/Manuscript/Figures/r1/r1_fig3_sims.pdf",
       height = 6,
       width = 8)

